use crate::chains::ChainExt;
use clap::Parser;
use std::path::PathBuf;

use crate::utils::load_model;
use std::path::Path;

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub(crate) struct Args {
    /// Path to the PDB or mmCIF file to be analyzed
    input: Vec<PathBuf>,
}

pub(crate) fn run(args: &Args) {
    for f in &args.input {
        // Make sure `input` exists
        let input_path = Path::new(f).canonicalize().unwrap();
        let input_file: String = input_path.to_str().unwrap().parse().unwrap();

        // Load file and print sequences
        let (pdb, _) = load_model(&input_file);
        println!("File: {}", input_file);
        for chain in pdb.chains() {
            println!("{}: {}", chain.id(), chain.pdb_seq().join(""));
        }
        println!();
    }
}
