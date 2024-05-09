#![warn(missing_docs)]
#![allow(clippy::type_complexity)]
#![doc = include_str!("../README.md")]

mod chains;
mod interactions;
mod residues;
mod utils;

use crate::chains::ChainExt;
use crate::interactions::{Interaction, InteractionComplex, Interactions, ResultEntry};
use crate::utils::load_model;

use clap::Parser;
use polars::prelude::*;
use std::path::{Path, PathBuf};
use tracing::{debug, info, trace, warn};

#[derive(Parser, Debug)]
#[command(version, about)]
struct Args {
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

/// Entry to the CLI tool. Verbosity can be controlled with the `RUST_LOG` environment variable.
fn main() {
    let args = Args::parse();
    tracing_subscriber::fmt::init();
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

    // Prepare output directory
    // let file_id = input_path.file_stem().unwrap().to_str().unwrap();
    let output_path = Path::new(&args.output);
    let _ = std::fs::create_dir_all(output_path);
    let output_dir = output_path.to_str().unwrap();
    debug!("Using input file {input_file}");
    debug!("Results will be saved to {output_dir}");

    // Load file as complex structure
    let (pdb, errors) = load_model(&input_file);
    if !errors.is_empty() {
        warn!("{errors:?}");
    }

    // Information on the sequence of the chains in the model
    info!("Loaded {} chains", pdb.chain_count());
    for chain in pdb.chains() {
        debug!(">{}: {}", chain.id(), chain.pdb_seq().join(""));
    }

    config_polars_output();
    let i_complex = InteractionComplex::new(pdb, &args.groups, args.vdw_comp, args.dist_cutoff);
    debug!(
        "Parsed ligand chains {lig:?}; receptor chains {receptor:?}",
        lig = i_complex.ligand,
        receptor = i_complex.receptor
    );

    // Find interactions
    let atomic_contacts = i_complex.get_atomic_contacts();
    atomic_contacts
        .iter()
        .filter(|x| x.interaction == Interaction::StericClash)
        .for_each(|x| warn!("{x}"));
    let mut df_atomic = results_to_df(&atomic_contacts);

    info!(
        "Found {} atom-atom contacts\n{}",
        atomic_contacts.len(),
        df_atomic
    );

    let mut ring_contacts: Vec<ResultEntry> = Vec::new();
    ring_contacts.extend(i_complex.get_ring_atom_contacts());
    ring_contacts.extend(i_complex.get_ring_ring_contacts());
    let mut df_ring = results_to_df(&ring_contacts)
        .drop_many(&["from_atomn", "from_atomi", "to_atomn", "to_atomi"])
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

    info!("Found {} ring contacts\n{}", ring_contacts.len(), df_ring);
    // ring_contacts.iter().for_each(|h| debug!("{h}"));

    // Save results to CSV files
    let mut file = std::fs::File::create(output_path.join("atomic_contacts.csv")).unwrap();
    CsvWriter::new(&mut file).finish(&mut df_atomic).unwrap();
    let mut file = std::fs::File::create(output_path.join("ring_contacts.csv")).unwrap();
    CsvWriter::new(&mut file).finish(&mut df_ring).unwrap();
}

fn config_polars_output() {
    std::env::set_var("POLARS_FMT_TABLE_HIDE_DATAFRAME_SHAPE_INFORMATION", "1");
    std::env::set_var("POLARS_FMT_TABLE_HIDE_COLUMN_DATA_TYPES", "1");
    std::env::set_var("POLARS_FMT_TABLE_ROUNDED_CORNERS", "1");
    std::env::set_var("POLARS_FMT_MAX_COLS", "14");
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
