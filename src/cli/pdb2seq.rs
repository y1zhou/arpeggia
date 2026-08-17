use arpeggia::ArpeggiaResult;
use clap::Parser;
use std::path::{Path, PathBuf};

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub(crate) struct Args {
    /// Path to the PDB or mmCIF file to be analyzed
    input: Vec<PathBuf>,
}

pub(crate) fn run(args: &Args) -> ArpeggiaResult<()> {
    for f in &args.input {
        // Make sure `input` exists
        let input_path = Path::new(f).canonicalize()?;
        let input_file = input_path.to_string_lossy();
        let pdb = super::load_input(&input_path)?;
        let sequences = arpeggia::get_sequences(&pdb);

        println!("File: {input_file}");
        for (chain_id, seq) in sequences {
            println!("{chain_id}: {seq}");
        }
        println!();
    }
    Ok(())
}

pub(crate) fn run_declared(args: &Args) -> ArpeggiaResult<()> {
    for input in &args.input {
        let input = input.canonicalize()?;
        let analysis = arpeggia::get_seqres(&input)?;
        for warning in analysis.warnings {
            tracing::warn!("{warning}");
        }
        println!("File: {}", input.display());
        for (chain_id, sequence) in analysis.value {
            println!("{chain_id}: {sequence}");
        }
        println!();
    }
    Ok(())
}
