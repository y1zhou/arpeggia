use arpeggia::{AlignmentMode, ArpeggiaError, ArpeggiaResult, SeqAlignOptions};
use clap::Parser;

#[derive(clap::Args, Debug, Clone)]
pub(crate) struct AlignmentArgs {
    /// Sequence alignment objective (semi-global consumes the second sequence)
    #[arg(long, alias = "mode", value_enum, default_value = "global")]
    alignment_mode: AlignmentMode,
    /// Positive gap opening cost, at most two decimal places
    #[arg(long, default_value_t = 10.0)]
    gap_open: f64,
    /// Positive extension cost, no larger than opening
    #[arg(long, default_value_t = 0.5)]
    gap_extend: f64,
}
impl AlignmentArgs {
    pub(crate) fn options(&self) -> SeqAlignOptions {
        SeqAlignOptions {
            mode: self.alignment_mode,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
        }
    }
}
#[derive(Parser, Debug, Clone)]
pub(crate) struct Args {
    /// First unaligned protein sequence
    reference: String,
    /// Second unaligned protein sequence
    mobile: String,
    #[command(flatten)]
    alignment: AlignmentArgs,
    /// Emit machine-readable results with parameters and mappings
    #[arg(long)]
    json: bool,
}
pub(crate) fn print_json(value: &impl serde::Serialize) -> ArpeggiaResult<()> {
    println!(
        "{}",
        serde_json::to_string(value).map_err(|e| ArpeggiaError::Calculation(e.to_string()))?
    );
    Ok(())
}
pub(crate) fn run(args: &Args) -> ArpeggiaResult<()> {
    let analysis = arpeggia::align_seqs(&args.reference, &args.mobile, &args.alignment.options())?;
    for warning in analysis.warnings {
        tracing::warn!("{warning}");
    }
    let a = analysis.value;
    if args.json {
        return print_json(&a);
    }
    println!(
        "Mode: {}\nScore: {}\nMatches: {} / {} columns / {} shorter-sequence residues\nPaired residues: {}\nGaps: {} residues in {} runs\nEdit distance: {}",
        a.mode,
        a.score,
        a.matches,
        a.alignment_length,
        a.shorter_length,
        a.paired_residues,
        a.gap_residues,
        a.gap_runs,
        a.edit_distance
    );
    println!(
        "Identity (alignment / shorter): {} / {}\nCoverage (alignment / shorter): {} / {}",
        ratio(a.identity_alignment),
        a.identity_shorter,
        ratio(a.coverage_alignment),
        a.coverage_shorter
    );
    Ok(())
}

pub(crate) fn ratio(value: Option<f64>) -> String {
    value.map_or_else(|| "undefined".into(), |v| format!("{v:.4}"))
}
