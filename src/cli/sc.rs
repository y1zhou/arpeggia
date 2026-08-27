use arpeggia::{ArpeggiaResult, run_with_threads};
use clap::Parser;
use std::path::{Path, PathBuf};
use tracing::{debug, info, trace};

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub(crate) struct Args {
    /// Path to the PDB or mmCIF file to be analyzed
    #[arg(short, long)]
    input: PathBuf,

    /// Chain groups for SC calculation:
    /// e.g. A,B/C,D where chains A and B form one surface and C and D form the other.
    /// Both groups must be specified for SC calculation.
    #[arg(short, long)]
    groups: String,

    /// Model number to analyze (default: 0, the first model)
    #[arg(short = 'm', long = "model", default_value_t = 0)]
    model_num: usize,

    /// Number of threads to use for parallel processing
    #[arg(short = 'j', long = "num-threads", default_value_t = 0)]
    num_threads: usize,
}

pub(crate) fn run(args: &Args) -> ArpeggiaResult<()> {
    trace!("{args:?}");

    // Validate groups
    if !args.groups.contains('/') {
        return Err(arpeggia::ArpeggiaError::InvalidArgument(
            "groups must be specified as 'A,B/C,D' with both surfaces defined".into(),
        ));
    }

    // Load file as complex structure
    let pdb = super::load_input(Path::new(&args.input))?;

    // Calculate SC
    let analysis = run_with_threads(args.num_threads as isize, || {
        debug!("Using {} thread(s)", rayon::current_num_threads());
        arpeggia::get_sc_details(&pdb, &args.groups, args.model_num)
    })?;
    for warning in analysis.warnings {
        tracing::warn!("{warning}");
    }
    info!("SC: {:.4}", analysis.value.sc);
    Ok(())
}
