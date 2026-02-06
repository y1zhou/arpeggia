use arpeggia::{get_dsasa, load_model, run_with_threads};
use clap::Parser;
use std::path::{Path, PathBuf};
use tracing::{debug, error, info, trace, warn};

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub(crate) struct Args {
    /// Path to the PDB or mmCIF file to be analyzed
    #[arg(short, long)]
    input: PathBuf,

    /// Group chains for interface calculation:
    /// e.g. A,B/C,D
    /// where chains A and B form one side and C and D form the other side.
    /// The buried surface area is calculated as the difference between the
    /// combined SASA and the sum of individual group SASAs.
    #[arg(short, long)]
    groups: String,

    /// Model number to analyze (default: 0, the first model)
    #[arg(short = 'm', long = "model", default_value_t = 0)]
    model_num: usize,

    /// Probe radius r (smaller r detects more surface details and reports a larger surface)
    #[arg(short = 'r', long = "probe-radius", default_value_t = 1.4)]
    probe_radius: f32,

    /// Number of points on the sphere for sampling
    #[arg(short = 'n', long = "num-points", default_value_t = 100)]
    n_points: usize,

    /// Number of threads to use for parallel processing
    #[arg(short = 'j', long = "num-threads", default_value_t = 1)]
    num_threads: usize,
}

pub(crate) fn run(args: &Args) {
    trace!("{args:?}");

    // Make sure `input` exists
    let input_path = Path::new(&args.input).canonicalize().unwrap();
    let input_file: String = input_path.to_str().unwrap().parse().unwrap();

    // Load file as complex structure
    let (pdb, pdb_warnings) = load_model(&input_file);
    if !pdb_warnings.is_empty() {
        for e in &pdb_warnings {
            match e.level() {
                pdbtbx::ErrorLevel::BreakingError => error!("{e}"),
                pdbtbx::ErrorLevel::InvalidatingError => error!("{e}"),
                _ => warn!("{e}"),
            }
        }
    }

    // Convert thread count to isize for rust-sasa
    let dsasa = run_with_threads(args.num_threads as isize, || {
        debug!("Using {} thread(s)", rayon::current_num_threads());
        get_dsasa(
            &pdb,
            &args.groups,
            args.probe_radius,
            args.n_points,
            args.model_num,
        )
    });
    info!(
        "Buried surface area (dSASA) at the interface between chains [{}]: {:.2} Å²",
        args.groups, dsasa
    );
}
