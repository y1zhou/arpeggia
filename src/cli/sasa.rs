use crate::interactions::InteractingEntity;
use crate::residues::ResidueExt;
use crate::utils::{load_model, write_df_to_file, DataFrameFileType};
use clap::Parser;
use pdbtbx::*;
use polars::prelude::*;
use rust_sasa::calculate_sasa_internal;
use rust_sasa::Atom as SASAAtom;
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

    /// Name of the output file
    #[arg(short, long, default_value_t = String::from("sasa"))]
    name: String,

    /// Output file type
    #[arg(short = 't', long, default_value_t = DataFrameFileType::Csv)]
    output_format: DataFrameFileType,

    /// Probe radius r (smaller r detects more surface details and reports a larger surface)
    #[arg(short = 'r', long = "probe-radius", default_value_t = 1.4)]
    probe_radius: f32,

    /// Distance cutoff when searching for neighboring atoms
    #[arg(short = 'n', long = "num-points", default_value_t = 100)]
    n_points: usize,

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

    let mut df_sasa = get_atom_sasa(&pdb, args.probe_radius, args.n_points);

    // Prepare output directory
    let output_path = Path::new(&args.output).canonicalize().unwrap();
    let _ = std::fs::create_dir_all(output_path.clone());
    let output_file = match output_path.is_dir() {
        true => output_path.join(args.name.clone()),
        false => output_path,
    }
    .with_extension(args.output_format.to_string());

    // Save results and log the identified SASA
    let non_zero_sasa_mask = df_sasa
        .column("sasa")
        .unwrap()
        .as_materialized_series()
        .not_equal(0.0)
        .unwrap();
    let df_sasa_nonzero = df_sasa.filter(&non_zero_sasa_mask).unwrap();
    debug!(
        "Found {} atoms with non-zero SASA\n{}",
        df_sasa_nonzero.shape().0,
        df_sasa_nonzero
    );

    // Save res to CSV files
    write_df_to_file(&mut df_sasa, &output_file, args.output_format);
    let output_file_str = output_file.to_str().unwrap();
    info!("Results saved to {output_file_str}");
}

pub fn get_atom_sasa(pdb: &PDB, probe_radius: f32, n_points: usize) -> DataFrame {
    // Calculate the SASA for each atom
    let atoms = pdb
        .atoms_with_hierarchy()
        .filter(|x| {
            let resn = x.residue().resn().unwrap();
            resn != "O" && resn != "X"
        })
        .map(|x| SASAAtom {
            position: nalgebra::Point3::new(
                x.atom().pos().0 as f32,
                x.atom().pos().1 as f32,
                x.atom().pos().2 as f32,
            ),
            radius: x
                .atom()
                .element()
                .unwrap()
                .atomic_radius()
                .van_der_waals
                .unwrap() as f32,
            id: x.atom().serial_number(),
            parent_id: None,
        })
        .collect::<Vec<_>>();
    let atom_sasa = calculate_sasa_internal(&atoms, Some(probe_radius), Some(n_points));

    // Create a DataFrame with the results
    let atom_annotations = pdb
        .atoms_with_hierarchy()
        .map(|x| InteractingEntity::from_hier(&x))
        .collect::<Vec<InteractingEntity>>();
    let atom_annot_df = df!(
        "chain" => atom_annotations.iter().map(|x| x.chain.to_owned()).collect::<Vec<String>>(),
        "resn" => atom_annotations.iter().map(|x| x.resn.to_owned()).collect::<Vec<String>>(),
        "resi" => atom_annotations.iter().map(|x| x.resi as i32).collect::<Vec<i32>>(),
        "altloc" => atom_annotations.iter().map(|x| x.altloc.to_owned()).collect::<Vec<String>>(),
        "atomn" => atom_annotations.iter().map(|x| x.atomn.to_owned()).collect::<Vec<String>>(),
        "atomi" => atom_annotations.iter().map(|x| x.atomi as i32).collect::<Vec<i32>>(),
    )
    .unwrap();

    df!(
        "atomi" => atoms.iter().map(|x| x.id as i32).collect::<Vec<i32>>(),
        "sasa" => atom_sasa
    )
    .unwrap()
    .join(
        &atom_annot_df,
        ["atomi"],
        ["atomi"],
        JoinArgs::new(JoinType::Inner),
        None,
    )
    .unwrap()
    .sort(["chain", "resi", "altloc", "atomi"], Default::default())
    .unwrap()
}
