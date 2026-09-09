use arpeggia::{ArpeggiaResult, AtomSubset, get_rmsd, validate_rmsd_selections};
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

    /// Residues used to determine the rigid-body transform
    #[arg(short = 's', long, default_value_t = String::new())]
    superpose_residues: String,

    /// Residues evaluated after applying the rigid-body transform
    #[arg(short = 'r', long, default_value_t = String::new())]
    rmsd_residues: String,

    /// Atom population used for fitting and RMSD
    #[arg(short = 'a', long, default_value = "ca")]
    atoms: AtomSubset,
}

pub(crate) fn run(args: &Args) -> ArpeggiaResult<()> {
    validate_rmsd_selections(&args.superpose_residues, &args.rmsd_residues)?;
    let reference = super::load_input(&args.reference)?;
    let mobile = super::load_input(&args.mobile)?;
    let analysis = get_rmsd(
        reference,
        mobile,
        &arpeggia::RmsdOptions {
            model_num: args.model_num,
            superpose_residues: args.superpose_residues.clone(),
            rmsd_residues: args.rmsd_residues.clone(),
            atoms: args.atoms,
            ..Default::default()
        },
    )?;
    for warning in analysis.warnings {
        tracing::warn!("{warning}");
    }
    println!("{}", analysis.value.rmsd);
    Ok(())
}
