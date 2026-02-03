use arpeggia::load_model;
use clap::Parser;
use std::path::{Path, PathBuf};
use tracing::{error, info, trace, warn};

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

    /// Number of threads to use for parallel calculations (default: 0 for auto)
    #[arg(short = 't', long = "threads", default_value_t = 0)]
    threads: usize,
}

pub(crate) fn run(args: &Args) {
    trace!("{args:?}");

    // Make sure `input` exists
    let input_path = match Path::new(&args.input).canonicalize() {
        Ok(path) => path,
        Err(e) => {
            error!("Failed to retrieve input file: {}", e);
            return;
        }
    };
    let input_file: String = input_path.to_str().unwrap().parse().unwrap();

    // Validate groups
    if !args.groups.contains('/') {
        error!("Groups must be specified as 'A,B/C,D' with both surfaces defined");
        return;
    }

    // Load file as complex structure
    let (pdb, pdb_warnings) = load_model(&input_file);
    if !pdb_warnings.is_empty() {
        pdb_warnings.iter().for_each(|e| match e.level() {
            pdbtbx::ErrorLevel::BreakingError => error!("{e}"),
            pdbtbx::ErrorLevel::InvalidatingError => error!("{e}"),
            _ => warn!("{e}"),
        });
    }

    // Calculate SC
    let sc = arpeggia::get_sc(&pdb, &args.groups, args.model_num, args.threads);

    match sc {
        Ok(score) => {
            info!("SC: {:.4}", score);
        }
        Err(e) => {
            error!("SC calculation failed: {}", e);
        }
    }
}
