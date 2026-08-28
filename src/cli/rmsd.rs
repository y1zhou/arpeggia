use arpeggia::{ArpeggiaResult, AtomSubset, get_rmsd};
use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub(crate) struct Args {
    /// First PDB or mmCIF structure
    reference: PathBuf,

    /// Second PDB or mmCIF structure
    mobile: PathBuf,

    /// Model number to select (0 selects the first model)
    #[arg(short = 'm', long = "model", default_value_t = 0)]
    model_num: usize,

    /// Comma-separated chain and author-residue ranges
    #[arg(short = 'r', long, default_value_t = String::new())]
    residues: String,

    /// Atom population used for fitting and RMSD
    #[arg(short = 'a', long, default_value = "ca")]
    atoms: AtomSubset,
}

pub(crate) fn run(args: &Args) -> ArpeggiaResult<()> {
    let reference = super::load_input(&args.reference)?;
    let mobile = super::load_input(&args.mobile)?;
    let analysis = get_rmsd(
        &reference,
        &mobile,
        args.model_num,
        &args.residues,
        args.atoms,
    )?;
    for warning in analysis.warnings {
        tracing::warn!("{warning}");
    }
    println!("{:.6}", analysis.value);
    Ok(())
}
