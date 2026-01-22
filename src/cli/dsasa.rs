use arpeggia::load_model;
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

    // Parse groups
    let groups_parts: Vec<&str> = args.groups.split('/').collect();
    if groups_parts.len() != 2 {
        error!("Invalid groups format. Expected format: 'A,B/C,D'");
        return;
    }

    let group1_chains: Vec<String> = groups_parts[0]
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.trim().to_string())
        .collect();
    let group2_chains: Vec<String> = groups_parts[1]
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.trim().to_string())
        .collect();

    if group1_chains.is_empty() || group2_chains.is_empty() {
        error!("Both groups must have at least one chain. Got: group1={:?}, group2={:?}", group1_chains, group2_chains);
        return;
    }

    debug!("Group 1 chains: {:?}", group1_chains);
    debug!("Group 2 chains: {:?}", group2_chains);

    // Calculate thread count: convert usize to isize, where 0 means use all cores (-1)
    let num_threads: isize = if args.num_threads == 0 {
        -1
    } else {
        args.num_threads as isize
    };

    // Calculate SASA for the combined complex
    let combined_sasa = arpeggia::get_chain_sasa(
        &pdb,
        args.probe_radius,
        args.n_points,
        args.model_num,
        num_threads,
    );

    if combined_sasa.is_empty() {
        error!("No SASA data found. Please check the input file and model number.");
        return;
    }

    // Get total SASA for all chains in both groups when together
    let combined_group_chains: Vec<String> = group1_chains
        .iter()
        .chain(group2_chains.iter())
        .cloned()
        .collect();

    let combined_total: f32 = combined_sasa
        .column("chain")
        .unwrap()
        .str()
        .unwrap()
        .iter()
        .zip(
            combined_sasa
                .column("sasa")
                .unwrap()
                .f32()
                .unwrap()
                .iter(),
        )
        .filter_map(|(chain_opt, sasa_opt)| {
            if let (Some(chain), Some(sasa)) = (chain_opt, sasa_opt) {
                if combined_group_chains.contains(&chain.to_string()) {
                    Some(sasa)
                } else {
                    None
                }
            } else {
                None
            }
        })
        .sum();

    // Create PDB with only group1 chains and calculate SASA
    let mut pdb_group1 = pdb.clone();
    pdb_group1.remove_chains_by(|chain| !group1_chains.contains(&chain.id().to_string()));

    let group1_sasa = arpeggia::get_chain_sasa(
        &pdb_group1,
        args.probe_radius,
        args.n_points,
        args.model_num,
        num_threads,
    );

    let group1_total: f32 = group1_sasa
        .column("sasa")
        .unwrap()
        .f32()
        .unwrap()
        .iter()
        .filter_map(|v| v)
        .sum();

    // Create PDB with only group2 chains and calculate SASA
    let mut pdb_group2 = pdb.clone();
    pdb_group2.remove_chains_by(|chain| !group2_chains.contains(&chain.id().to_string()));

    let group2_sasa = arpeggia::get_chain_sasa(
        &pdb_group2,
        args.probe_radius,
        args.n_points,
        args.model_num,
        num_threads,
    );

    let group2_total: f32 = group2_sasa
        .column("sasa")
        .unwrap()
        .f32()
        .unwrap()
        .iter()
        .filter_map(|v| v)
        .sum();

    // Calculate buried surface area (dSASA)
    // BSA = (SASA_group1 + SASA_group2 - SASA_complex) / 2
    // This gives the buried surface area at the interface
    let dsasa = (group1_total + group2_total - combined_total) / 2.0;

    debug!("SASA of group 1 ({:?}): {:.2} Å²", group1_chains, group1_total);
    debug!("SASA of group 2 ({:?}): {:.2} Å²", group2_chains, group2_total);
    debug!("SASA of combined complex: {:.2} Å²", combined_total);

    info!("{:.2}", dsasa);
}
