mod chains;
mod interactions;
mod residues;
mod utils;

use crate::chains::ChainExt;
use crate::interactions::complex::{InteractionComplex, Interactions};
use crate::utils::load_model;
use clap::Parser;
use std::path::{Path, PathBuf};

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
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

    /// Distance cutoff for grid points to be 'interacting' with the entity
    #[arg(short, long, default_value_t = 5.0)]
    dist_cutoff: f64,

    /// Verbosity of the program:
    /// -v for info, -vv for debug, and -vvv for trace
    #[arg(short, long, action=clap::ArgAction::Count)]
    verbose: u8,
}

fn main() {
    let args = Args::parse();
    if args.verbose > 2 {
        println!("{args:?}");
    }

    // Make sure `input` exists
    let input_path = Path::new(&args.input).canonicalize().unwrap();
    let input_file: String = input_path.to_str().unwrap().parse().unwrap();

    // Prepare output directory
    // let file_id = input_path.file_stem().unwrap().to_str().unwrap();
    let output_path = Path::new(&args.output).canonicalize().unwrap();
    let _ = std::fs::create_dir_all(&output_path);
    let output_dir = output_path.to_str().unwrap();
    if args.verbose > 1 {
        println!("Using input file {input_file}");
        println!("Results will be saved to {output_dir}");
    }

    // Load file as complex structure
    let (pdb, errors) = load_model(&input_file);
    if args.verbose > 0 {
        for err in errors {
            eprintln!("{err}");
        }
    }

    // Information on the sequence of the chains in the model
    if args.verbose > 1 {
        println!("Loaded {} chains", pdb.chain_count());
        for chain in pdb.chains() {
            println!(">{}\n{}", chain.id(), chain.pdb_seq().join(""));
        }
    }

    let i_complex = InteractionComplex::new(pdb, &args.groups, args.vdw_comp, args.dist_cutoff);

    i_complex.get_atomic_contacts();
}
