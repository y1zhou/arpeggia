use crate::interactions::{InteractionComplex, Interactions, ResultEntry};
use crate::residues::ResidueId;
use crate::utils::{load_model, write_df_to_file, DataFrameFileType};
use clap::Parser;
use pdbtbx::*;
use polars::prelude::*;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use tracing::{debug, error, info, trace, warn};

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub(crate) struct Args {
    /// Path to the PDB or mmCIF file to be analyzed
    #[arg(short, long)]
    input: PathBuf,

    /// Output directory
    #[arg(short, long)]
    output: PathBuf,

    /// Group chains for interactions:
    /// e.g. A,B/C,D
    /// where chains A and B are the "ligand" and C and D are the "receptor".
    /// Chains can exist on both sides, in which case intra-chain interactions will be calculated.
    /// If only one group is provided, all remaining chains will be considered as the other group.
    /// If no groups are provided ('/'), all inter- and intra-chain interactions will be calculated.
    #[arg(short, long, default_value_t = String::from("/"))]
    groups: String,

    /// Name of the output file
    #[arg(short, long, default_value_t = String::from("contacts"))]
    name: String,

    /// Output file type
    #[arg(short = 't', long, default_value_t = DataFrameFileType::Csv)]
    output_format: DataFrameFileType,

    /// Compensation factor for VdW radii dependent interaction types
    #[arg(short = 'c', long = "vdw-comp", default_value_t = 0.1)]
    vdw_comp: f64,

    /// Distance cutoff when searching for neighboring atoms
    #[arg(short, long, default_value_t = 6.5)]
    dist_cutoff: f64,

    /// Number of threads to use for parallel processing. One thread should be sufficient unless the system is very large
    #[arg(short = 'j', long = "num-threads", default_value_t = 1)]
    num_threads: usize,
}

pub(crate) fn run(args: &Args) {
    trace!("{args:?}");

    // Create Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.num_threads)
        .build_global()
        .unwrap();
    debug!("Using {} thread(s)", rayon::current_num_threads());

    // Make sure `input` exists
    let input_path = match Path::new(&args.input).canonicalize() {
        Ok(path) => path,
        Err(e) => {
            error!("Failed to retrieve input file: {}", e);
            return;
        }
    };
    let output_path = match std::path::absolute(&args.output) {
        Ok(path) => path,
        Err(e) => {
            error!("Failed to resolve the output directory: {}", e);
            return;
        }
    };
    let input_file: String = input_path.to_str().unwrap().parse().unwrap();

    // Load file as complex structure
    let (pdb, pdb_warnings) = load_model(&input_file);
    if !pdb_warnings.is_empty() {
        pdb_warnings.iter().for_each(|e| match e.level() {
            pdbtbx::ErrorLevel::BreakingError => error!("{e}"),
            pdbtbx::ErrorLevel::InvalidatingError => error!("{e}"),
            _ => warn!("{e}"),
        });
    }

    let mut df_contacts = get_contacts(&pdb, args.groups.as_str(), args.vdw_comp, args.dist_cutoff);

    // Prepare output directory
    let _ = std::fs::create_dir_all(output_path.clone());
    let output_file = output_path
        .join(args.name.clone())
        .with_extension(args.output_format.to_string());

    // Save results and log the identified interactions
    let df_clash = df_contacts
        .clone()
        .lazy()
        .filter(col("interaction").eq(lit("StericClash")))
        .collect()
        .unwrap();
    if df_clash.height() > 0 {
        warn!(
            "Found {} steric {}\n{}",
            df_clash.height(),
            match df_clash.height() {
                1 => "clash",
                _ => "clashes",
            },
            df_clash
        );
    }

    // Save res to CSV files
    write_df_to_file(&mut df_contacts, &output_file, args.output_format);
    let output_file_str = output_file.to_str().unwrap();
    info!("Results saved to {output_file_str}");
}

pub fn get_contacts<'a>(
    pdb: &'a PDB,
    groups: &'a str,
    vdw_comp: f64,
    dist_cutoff: f64,
) -> DataFrame {
    let (i_complex, build_ring_err) =
        InteractionComplex::new(pdb, groups, vdw_comp, dist_cutoff).unwrap();

    if !build_ring_err.is_empty() {
        build_ring_err.iter().for_each(|e| warn!("{e}"));
    }

    // Information on the sequence of the chains in the model
    debug!(
        "Parsed ligand chains {lig:?}; receptor chains {receptor:?}",
        lig = i_complex.ligand,
        receptor = i_complex.receptor
    );

    // Find interactions
    let atomic_contacts = i_complex.get_atomic_contacts();
    let df_atomic = results_to_df(&atomic_contacts);
    debug!(
        "Found {} atom-atom contacts\n{}",
        df_atomic.height(),
        df_atomic
    );

    let mut ring_contacts: Vec<ResultEntry> = Vec::new();
    ring_contacts.extend(i_complex.get_ring_atom_contacts());
    ring_contacts.extend(i_complex.get_ring_ring_contacts());
    let df_ring = results_to_df(&ring_contacts);
    debug!("Found {} ring contacts\n{}", df_ring.height(), df_ring);

    // Annotate sidechain centroid distances and dihedrals
    let mut sc_dist_dihedrals = HashMap::new();
    sc_dist_dihedrals.extend(i_complex.collect_sc_stats(&atomic_contacts));
    sc_dist_dihedrals.extend(i_complex.collect_sc_stats(&ring_contacts));

    let df_sc_stats = sc_results_to_df(&sc_dist_dihedrals);

    // Concate dataframes
    let contacts_sc_join_cols = [
        col("model"),
        col("from_chain"),
        col("from_resi"),
        col("from_insertion"),
        col("from_altloc"),
        col("to_chain"),
        col("to_resi"),
        col("to_insertion"),
        col("to_altloc"),
    ];
    concat([df_atomic.lazy(), df_ring.lazy()], UnionArgs::default())
        .unwrap()
        // Annotate sidechain stats
        .join(
            df_sc_stats.lazy(),
            &contacts_sc_join_cols,
            &contacts_sc_join_cols,
            JoinArgs::new(JoinType::Left),
        )
        .sort(
            [
                "model",
                "from_chain",
                "to_chain",
                "from_resi",
                "from_altloc",
                "from_atomi",
                "to_resi",
                "to_altloc",
                "to_atomi",
                "interaction",
            ],
            Default::default(),
        )
        .collect()
        .unwrap()
}

fn results_to_df(res: &[ResultEntry]) -> DataFrame {
    df!(
        "model" => res.iter().map(|x| x.model as u32).collect::<Vec<u32>>(),
        "interaction" => res.iter().map(|x| x.interaction.to_string()).collect::<Vec<String>>(),
        "distance" => res.iter().map(|x| x.distance as f32).collect::<Vec<f32>>(),
        "from_chain" => res.iter().map(|x| x.ligand.chain.to_owned()).collect::<Vec<String>>(),
        "from_resn" => res.iter().map(|x| x.ligand.resn.to_owned()).collect::<Vec<String>>(),
        "from_resi" => res.iter().map(|x| x.ligand.resi as i32).collect::<Vec<i32>>(),
        "from_insertion" => res.iter().map(|x| x.ligand.insertion.to_owned()).collect::<Vec<String>>(),
        "from_altloc" => res.iter().map(|x| x.ligand.altloc.to_owned()).collect::<Vec<String>>(),
        "from_atomn" => res.iter().map(|x| x.ligand.atomn.to_owned()).collect::<Vec<String>>(),
        "from_atomi" => res.iter().map(|x| x.ligand.atomi as i32).collect::<Vec<i32>>(),
        "to_chain" => res.iter().map(|x| x.receptor.chain.to_owned()).collect::<Vec<String>>(),
        "to_resn" => res.iter().map(|x| x.receptor.resn.to_owned()).collect::<Vec<String>>(),
        "to_resi" => res.iter().map(|x| x.receptor.resi as i32).collect::<Vec<i32>>(),
        "to_insertion" => res.iter().map(|x| x.receptor.insertion.to_owned()).collect::<Vec<String>>(),
        "to_altloc" => res.iter().map(|x| x.receptor.altloc.to_owned()).collect::<Vec<String>>(),
        "to_atomn" => res.iter().map(|x| x.receptor.atomn.to_owned()).collect::<Vec<String>>(),
        "to_atomi" => res.iter().map(|x| x.receptor.atomi as i32).collect::<Vec<i32>>(),
    )
    .unwrap()
}

fn sc_results_to_df(res: &HashMap<(ResidueId, ResidueId), (f64, f64, f64)>) -> DataFrame {
    df!(
        "model" => res.iter().map(|(k, _)| k.0.model as u32).collect::<Vec<u32>>(),
        "from_chain" => res.iter().map(|(k, _)| k.0.chain.to_owned()).collect::<Vec<String>>(),
        "from_resi" => res.iter().map(|(k, _)| k.0.resi as i32).collect::<Vec<i32>>(),
        "from_insertion" => res.iter().map(|(k, _)| k.0.insertion.to_owned()).collect::<Vec<String>>(),
        "from_altloc" => res.iter().map(|(k, _)| k.0.altloc.to_owned()).collect::<Vec<String>>(),
        "to_chain" => res.iter().map(|(k, _)| k.1.chain.to_owned()).collect::<Vec<String>>(),
        "to_resi" => res.iter().map(|(k, _)| k.1.resi as i32).collect::<Vec<i32>>(),
        "to_insertion" => res.iter().map(|(k, _)| k.1.insertion.to_owned()).collect::<Vec<String>>(),
        "to_altloc" => res.iter().map(|(k, _)| k.1.altloc.to_owned()).collect::<Vec<String>>(),
        "sc_centroid_dist" => res.iter().map(|(_, v)| v.0 as f32).collect::<Vec<f32>>(),
        "sc_dihedral" => res.iter().map(|(_, v)| v.1 as f32).collect::<Vec<f32>>(),
        "sc_centroid_angle" => res.iter().map(|(_, v)| v.2 as f32).collect::<Vec<f32>>(),
    )
    .unwrap()
}
