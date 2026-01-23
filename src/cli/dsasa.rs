use arpeggia::{load_model, parse_groups};
use clap::Parser;
use polars::prelude::*;
use std::collections::HashSet;
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

/// Sum the SASA column of a DataFrame using polars aggregation.
fn sum_sasa(df: &DataFrame) -> f32 {
    df.clone()
        .lazy()
        .select([col("sasa").sum()])
        .collect()
        .unwrap()
        .column("sasa")
        .unwrap()
        .f32()
        .unwrap()
        .get(0)
        .unwrap_or(0.0)
}

pub(crate) fn run(args: &Args) {
    trace!("{args:?}");

    // Create Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.num_threads)
        .build_global()
        .unwrap();
    debug!("Using {} thread(s)", rayon::current_num_threads());

    // Make sure `input` exists
    let input_path = Path::new(&args.input).canonicalize().unwrap();
    let input_file: String = input_path.to_str().unwrap().parse().unwrap();

    // Load file as complex structure
    let (pdb, pdb_warnings) = load_model(&input_file);
    if !pdb_warnings.is_empty() {
        pdb_warnings.iter().for_each(|e| match e.level() {
            pdbtbx::ErrorLevel::BreakingError => error!("{e}"),
            pdbtbx::ErrorLevel::InvalidatingError => error!("{e}"),
            _ => warn!("{e}"),
        });
    }

    // Get all chains in the PDB
    let all_chains: HashSet<String> = pdb.chains().map(|c| c.id().to_string()).collect();

    // Parse groups using the utility function
    let (group1_chains, group2_chains) = parse_groups(&all_chains, &args.groups);

    debug!("Group 1 chains: {:?}", group1_chains);
    debug!("Group 2 chains: {:?}", group2_chains);

    // Calculate thread count: convert usize to isize, where 0 means use all cores (-1)
    let num_threads: isize = if args.num_threads == 0 {
        -1
    } else {
        args.num_threads as isize
    };

    // Get combined chains (union of both groups)
    let combined_group_chains: HashSet<String> =
        group1_chains.union(&group2_chains).cloned().collect();

    // Create PDB with only chains from both groups (remove unrelated chains)
    let mut pdb_combined = pdb.clone();
    pdb_combined.remove_chains_by(|chain| !combined_group_chains.contains(chain.id()));

    // Calculate SASA for the combined complex (only chains in groups)
    let combined_sasa = arpeggia::get_chain_sasa(
        &pdb_combined,
        args.probe_radius,
        args.n_points,
        args.model_num,
        num_threads,
    );

    if combined_sasa.is_empty() {
        error!("No SASA data found. Please check the input file and model number.");
        return;
    }

    // Get total SASA for all chains in both groups when together using polars
    let combined_total = sum_sasa(&combined_sasa);

    // Create PDB with only group1 chains and calculate SASA
    let mut pdb_group1 = pdb.clone();
    pdb_group1.remove_chains_by(|chain| !group1_chains.contains(chain.id()));

    let group1_sasa = arpeggia::get_chain_sasa(
        &pdb_group1,
        args.probe_radius,
        args.n_points,
        args.model_num,
        num_threads,
    );

    let group1_total = sum_sasa(&group1_sasa);

    // Create PDB with only group2 chains and calculate SASA
    let mut pdb_group2 = pdb.clone();
    pdb_group2.remove_chains_by(|chain| !group2_chains.contains(chain.id()));

    let group2_sasa = arpeggia::get_chain_sasa(
        &pdb_group2,
        args.probe_radius,
        args.n_points,
        args.model_num,
        num_threads,
    );

    let group2_total = sum_sasa(&group2_sasa);

    // Calculate buried surface area (dSASA)
    // dSASA = SASA_group1 + SASA_group2 - SASA_complex
    // This gives the buried surface area at the interface
    let dsasa = group1_total + group2_total - combined_total;

    debug!(
        "SASA of group 1 ({:?}): {:.2} Å²",
        group1_chains, group1_total
    );
    debug!(
        "SASA of group 2 ({:?}): {:.2} Å²",
        group2_chains, group2_total
    );
    debug!("SASA of combined complex: {:.2} Å²", combined_total);

    let group1_str: Vec<String> = group1_chains.iter().map(|c| c.to_string()).collect();
    let group2_str: Vec<String> = group2_chains.iter().map(|c| c.to_string()).collect();
    info!(
        "Buried surface area (dSASA) at the interface between chains [{}] and [{}]: {:.2} Å²",
        group1_str.join(", "),
        group2_str.join(", "),
        dsasa
    );
}
