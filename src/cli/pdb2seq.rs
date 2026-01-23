use arpeggia::load_model;
use clap::Parser;
use std::path::{Path, PathBuf};

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
        let sequences = arpeggia::get_sequences(&pdb);

        println!("File: {}", input_file);
        for (chain_id, seq) in sequences {
            println!("{}: {}", chain_id, seq);
        }
        println!();
    }
}
