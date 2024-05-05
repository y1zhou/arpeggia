#![warn(missing_docs)]
#![allow(clippy::type_complexity)]
#![doc = include_str!("../README.md")]

mod chains;
mod interactions;
mod residues;
mod utils;

use crate::chains::ChainExt;
use crate::interactions::{Interaction, InteractionComplex, Interactions};
use crate::utils::load_model;

use clap::Parser;
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

    let i_complex = InteractionComplex::new(pdb, &args.groups, args.vdw_comp, args.dist_cutoff);
    debug!(
        "Parsed ligand chains {lig:?}; receptor chains {receptor:?}",
        lig = i_complex.ligand,
        receptor = i_complex.receptor
    );

    let atomic_contacts = i_complex.get_atomic_contacts();
    info!("Found {} atom-atom contacts", atomic_contacts.len());
    atomic_contacts.iter().for_each(|h| {
        match h.interaction {
            Interaction::StericClash => warn!("{h}"),
            _ => debug!("{h}"),
        };
    })
}
