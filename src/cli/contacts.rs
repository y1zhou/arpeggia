use crate::chains::ChainExt;
use crate::interactions::{InteractionComplex, Interactions, ResultEntry};
use crate::utils::{load_model, write_df_to_csv};
use clap::Parser;
use pdbtbx::*;
use polars::prelude::*;
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
    /// where chains A and B are the "ligand" and C and D are the "receptor"
    #[arg(short, long)]
    groups: String,

    /// Compensation factor for VdW radii dependent interaction types
    #[arg(short = 'c', long = "vdw-comp", default_value_t = 0.1)]
    vdw_comp: f64,

    /// Distance cutoff when searching for neighboring atoms
    #[arg(short, long, default_value_t = 4.5)]
    dist_cutoff: f64,

    /// Number of threads to use for parallel processing
    #[arg(short = 'j', long = "num-threads", default_value_t = 0)]
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
    let input_path = Path::new(&args.input).canonicalize().unwrap();
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

    let (df_atomic, df_ring, i_complex) =
        get_contacts(pdb, args.groups.as_str(), args.vdw_comp, args.dist_cutoff);

    // Information on the sequence of the chains in the model
    info!(
        "Loaded {} {}",
        i_complex.model.chain_count(),
        match i_complex.model.chain_count() {
            1 => "chain",
            _ => "chains",
        }
    );
    for chain in i_complex.model.chains() {
        debug!(">{}: {}", chain.id(), chain.pdb_seq().join(""));
    }

    debug!(
        "Parsed ligand chains {lig:?}; receptor chains {receptor:?}",
        lig = i_complex.ligand,
        receptor = i_complex.receptor
    );

    // Prepare output directory
    // let file_id = input_path.file_stem().unwrap().to_str().unwrap();
    let output_path = Path::new(&args.output).canonicalize().unwrap();
    let _ = std::fs::create_dir_all(output_path.clone());
    let output_dir = output_path.to_str().unwrap();

    debug!("Results will be saved to {output_dir}/contacts.csv");

    // Save results and log the identified interactions

    info!(
        "Found {} atom-atom contacts\n{}",
        df_atomic.shape().0,
        df_atomic
    );
    let df_clash = df_atomic
        .clone()
        .lazy()
        .filter(col("interaction").eq(lit("StericClash")))
        .collect()
        .unwrap();
    if df_clash.height() > 0 {
        warn!("Found {} steric clashes\n{}", df_clash.shape().0, df_clash);
    }
    info!("Found {} ring contacts\n{}", df_ring.shape().0, df_ring);

    // Concate dataframes for saving to CSV
    let mut df_contacts = concat(
        [df_atomic.clone().lazy(), df_ring.clone().lazy()],
        UnionArgs::default(),
    )
    .unwrap()
    .collect()
    .unwrap();

    // Save res to CSV files
    write_df_to_csv(&mut df_contacts, output_path.join("contacts.csv"));
}

pub fn get_contacts(
    pdb: PDB,
    groups: &str,
    vdw_comp: f64,
    dist_cutoff: f64,
) -> (DataFrame, DataFrame, InteractionComplex) {
    let i_complex = InteractionComplex::new(pdb, groups, vdw_comp, dist_cutoff);

    // Find interactions
    let atomic_contacts = i_complex.get_atomic_contacts();
    let df_atomic = results_to_df(&atomic_contacts);

    let mut ring_contacts: Vec<ResultEntry> = Vec::new();
    ring_contacts.extend(i_complex.get_ring_atom_contacts());
    ring_contacts.extend(i_complex.get_ring_ring_contacts());
    let df_ring = results_to_df(&ring_contacts)
        // .drop_many(&["from_atomn", "from_atomi", "to_atomn", "to_atomi"])
        .sort(
            [
                "from_chain",
                "from_resi",
                "from_altloc",
                "to_chain",
                "to_resi",
                "to_altloc",
            ],
            Default::default(),
        )
        .unwrap();

    (df_atomic, df_ring, i_complex)
}

fn results_to_df(res: &[ResultEntry]) -> DataFrame {
    df!(
        "interaction" => res.iter().map(|x| x.interaction.to_string()).collect::<Vec<String>>(),
        "distance" => res.iter().map(|x| x.distance).collect::<Vec<f64>>(),
        "from_chain" => res.iter().map(|x| x.ligand.chain.to_owned()).collect::<Vec<String>>(),
        "from_resn" => res.iter().map(|x| x.ligand.resn.to_owned()).collect::<Vec<String>>(),
        "from_resi" => res.iter().map(|x| x.ligand.resi as i64).collect::<Vec<i64>>(),
        "from_altloc" => res.iter().map(|x| x.ligand.altloc.to_owned()).collect::<Vec<String>>(),
        "from_atomn" => res.iter().map(|x| x.ligand.atomn.to_owned()).collect::<Vec<String>>(),
        "from_atomi" => res.iter().map(|x| x.ligand.atomi as i64).collect::<Vec<i64>>(),
        "to_chain" => res.iter().map(|x| x.receptor.chain.to_owned()).collect::<Vec<String>>(),
        "to_resn" => res.iter().map(|x| x.receptor.resn.to_owned()).collect::<Vec<String>>(),
        "to_resi" => res.iter().map(|x| x.receptor.resi as i64).collect::<Vec<i64>>(),
        "to_altloc" => res.iter().map(|x| x.receptor.altloc.to_owned()).collect::<Vec<String>>(),
        "to_atomn" => res.iter().map(|x| x.receptor.atomn.to_owned()).collect::<Vec<String>>(),
        "to_atomi" => res.iter().map(|x| x.receptor.atomi as i64).collect::<Vec<i64>>(),
    )
    .unwrap()
}
