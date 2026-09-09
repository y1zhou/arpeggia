use arpeggia::{
    ArpeggiaError, ArpeggiaResult, AtomSubset, RmsdOptions, get_rmsd, validate_rmsd_selections,
};
use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug, Clone)]
#[command(version, about, after_help = super::RESIDUE_SELECTION_HELP)]
pub(crate) struct Args {
    /// First PDB or mmCIF structure
    reference: PathBuf,

    /// Second PDB or mmCIF structure
    query: PathBuf,

    /// Model serial to select (0 independently selects the first model of each structure)
    #[arg(short = 'm', long = "model", default_value_t = 0)]
    model_num: usize,

    /// Fit residues: comma-separated chains/ranges, e.g. A:1-100,B (empty: all)
    #[arg(short = 's', long, default_value_t = String::new())]
    superpose_residues: String,

    /// Evaluation residues, e.g. A:101-120 (empty: all; independent of fit selection)
    #[arg(short = 'r', long, default_value_t = String::new())]
    rmsd_residues: String,

    /// Atom population for fit/evaluation: ca, backbone (N/CA/C/O/OXT), heavy, or all
    #[arg(short = 'a', long, default_value = "ca")]
    atoms: AtomSubset,
    /// Align observed sequences before applying reference residue selections
    #[arg(long)]
    align_seqs: bool,
    /// Explicit reference=query chain pair; repeat for every selected chain
    #[arg(long,requires="align_seqs",value_parser=parse_chain_pair)]
    chain_map: Vec<(String, String)>,
    #[command(flatten)]
    alignment: super::align_seqs::AlignmentArgs,
    #[command(flatten)]
    display: super::align_seqs::DisplayArgs,
    /// Maximum rejection/refit cycles after the initial fit
    #[arg(long, default_value_t = 0)]
    refine_cycles: usize,
    /// Reject distances exceeding this multiple of current fitting RMSD
    #[arg(long, default_value_t = 2.0)]
    refine_cutoff: f64,
    /// Emit machine-readable results with parameters and mappings
    #[arg(long)]
    json: bool,
}

pub(crate) fn run(args: &Args) -> ArpeggiaResult<()> {
    validate_rmsd_selections(&args.superpose_residues, &args.rmsd_residues)?;
    let reference = super::load_input(&args.reference)?;
    let query = super::load_input(&args.query)?;
    let mut chain_map = std::collections::BTreeMap::new();
    for (a, b) in &args.chain_map {
        if chain_map.insert(a.clone(), b.clone()).is_some() {
            return Err(ArpeggiaError::InvalidArgument(format!(
                "duplicate reference chain mapping {a}"
            )));
        }
    }
    let options = RmsdOptions {
        model_num: args.model_num,
        superpose_residues: args.superpose_residues.clone(),
        rmsd_residues: args.rmsd_residues.clone(),
        atoms: args.atoms,
        align_seqs: args.align_seqs,
        chain_map,
        alignment: args.alignment.options(),
        refine_cycles: args.refine_cycles,
        refine_cutoff: args.refine_cutoff,
    };
    let analysis = get_rmsd(reference, query, &options)?;
    for warning in analysis.warnings {
        tracing::warn!("{warning}");
    }
    let result = analysis.value;
    if args.json {
        return super::align_seqs::print_json(&result);
    }
    println!(
        "RMSD: {} Å ({} evaluation pairs)\nCore RMSD: {} Å ({} retained / {} initial fitting pairs)\nInitial fitting RMSD: {} Å\nRefinement cycles: {}",
        result.rmsd,
        result.evaluation_atoms,
        result.core_rmsd,
        result.retained_fit_atoms,
        result.initial_fit_atoms,
        result.initial_rmsd,
        result.cycles
    );
    for chain in result.chain_alignments {
        println!(
            "Chain {} -> {}: {} corresponding residues",
            chain.reference_chain,
            chain.query_chain,
            chain.residue_pairs.len()
        );
        if let Some(a) = chain.alignment {
            println!("{}", args.display.format(&a)?);
        }
    }
    Ok(())
}

fn parse_chain_pair(value: &str) -> Result<(String, String), String> {
    value
        .split_once('=')
        .filter(|(_, b)| !b.contains('='))
        .map(|(a, b)| (a.to_owned(), b.to_owned()))
        .ok_or_else(|| "chain mapping must be reference=query".into())
}
