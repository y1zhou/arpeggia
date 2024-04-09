mod chains;
mod interactions;
mod residues;
mod utils;

use crate::utils::load_model;
use chains::ChainExt;
use clap::Parser;
use interactions::InteractionComplex;
use std::path::{Path, PathBuf};

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Path to the file to be analyzed
    #[arg(short, long)]
    filename: PathBuf,

    /// Output directory
    #[arg(short, long)]
    output: String,

    /// Select the "ligand" for interactions, using selection syntax:
    /// /<chain_id>/<res_num>[<ins_code>]/<atom_name> or
    /// RESNAME:<het_id>. Fields can be omitted.
    #[arg(short, long, num_args=0..)]
    selection: Vec<String>,

    /// Compensation factor for VdW radii dependent interaction types
    #[arg(short = 'c', long = "vdw-comp", default_value_t = 0.1)]
    vdw_comp: f64,

    /// Distance cutoff for grid points to be 'interacting' with the entity
    #[arg(short, long, default_value_t = 5.0)]
    dist_cutoff: f64,

    /// Verbosity of the program
    #[arg(short, long, action=clap::ArgAction::Count)]
    verbose: u8,
}

fn main() {
    let args = Args::parse();

    // Make sure `filename` exists
    let input_path = Path::new(&args.filename).canonicalize().unwrap();
    let input_file: String = input_path.to_str().unwrap().parse().unwrap();

    // Prepare output directory
    // let file_id = input_path.file_stem().unwrap().to_str().unwrap();
    let output_path = Path::new(&args.output).canonicalize().unwrap();
    let _ = std::fs::create_dir_all(&output_path);
    let output_dir = output_path.to_str().unwrap();
    if args.verbose > 0 {
        println!("Using input file {input_file}");
        println!("Saving results to {output_dir}");
    }

    // Load file as complex structure
    let (pdb, errors) = load_model(&input_file);
    if args.verbose > 1 {
        for err in errors {
            eprintln!("{err}");
        }
    }

    if args.verbose > 0 {
        println!("Loaded {} chains", pdb.chain_count());
        for chain in pdb.chains() {
            println!(">{}\n{}", chain.id(), chain.pdb_seq().join(""));
        }
    }

    let _i_complex = InteractionComplex {
        model: pdb,
        vdw_comp_factor: args.vdw_comp,
        interacting_threshold: args.dist_cutoff,
    };
}
