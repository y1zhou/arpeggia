mod utils;

use crate::utils::load_model;
use crate::utils::ResidueExt;
use clap::Parser;
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
    let file_id = input_path.file_stem().unwrap();
    let output_path = Path::new(&args.output)
        .canonicalize()
        .unwrap()
        .join(file_id);
    let _ = std::fs::create_dir_all(&output_path);
    let output_dir: String = output_path.to_str().unwrap().parse().unwrap();
    if args.verbose > 0 {
        println!("Using input file {input_file}");
        println!("Saving results to {output_dir}");
    }

    // Load file as complex structure
    let (mut pdb, errors) = load_model(&input_file);
    if args.verbose > 1 {
        for err in errors {
            eprintln!("{err}");
        }
    }

    // Load the amino acid sequence for each chain
    for chain in pdb.chains() {
        let chain_seq: Vec<_> = chain.residues().map(|res| res.resn().unwrap()).collect();

        println!(">{}\n{}", chain.id(), chain_seq.join(""));
    }

    // let _ = pdbtbx::save(
    //     &pdb,
    //     format!("{output_dir}/{file_id}.pdb"),
    //     pdbtbx::StrictnessLevel::Loose,
    // );
}
