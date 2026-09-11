use arpeggia::{AlignmentColor, AlignmentMode, ArpeggiaError, ArpeggiaResult, SeqAlignOptions};
use clap::Parser;

#[derive(clap::Args, Debug, Clone)]
pub(crate) struct AlignmentArgs {
    /// Sequence alignment objective (semi-global consumes the second sequence)
    #[arg(long, alias = "mode", value_enum, default_value = "global")]
    alignment_mode: AlignmentMode,
    /// Positive first-gap-residue cost, at most two decimals; gap cost = open + (length-1)*extend
    #[arg(long, default_value_t = 10.0)]
    gap_open: f64,
    /// Positive cost per additional gap residue, at most two decimals and no larger than opening
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
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct DisplayArgs {
    /// Total alignment line width, including labels (default: terminal width or 80)
    #[arg(long)]
    pub(crate) width: Option<usize>,
    /// Color policy; auto respects terminal capability and NO_COLOR
    #[arg(long, value_enum, default_value = "auto")]
    pub(crate) color: AlignmentColor,
    /// Hide position rulers; retain sequence start/end numbers
    #[arg(long)]
    pub(crate) no_rulers: bool,
}
impl DisplayArgs {
    pub(crate) fn format(&self, alignment: &arpeggia::SeqAlignment) -> ArpeggiaResult<String> {
        alignment.format(self.width, self.color, !self.no_rulers)
    }
}
#[derive(Parser, Debug, Clone)]
#[command(
    after_help = "Example: arpeggia align-seqs GGACDEFGHIKGG ACDEFGHIK --mode semi-global\nInputs are amino-acid strings, without whitespace, gaps, or stop symbols.\nSemi-global consumes the complete second sequence with free first-sequence tails.\nIdentity/coverage are ratios, reported using alignment and shorter-input lengths.\nUse --json for plain structured data; display controls do not affect alignment."
)]
pub(crate) struct Args {
    /// First unaligned protein sequence
    reference: String,
    /// Second unaligned protein sequence
    query: String,
    /// Display name of the reference sequence
    #[arg(long, default_value = "Reference")]
    reference_name: String,
    /// Display name of the query sequence
    #[arg(long, default_value = "Query")]
    query_name: String,
    #[command(flatten)]
    alignment: AlignmentArgs,
    #[command(flatten)]
    display: DisplayArgs,
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
    let analysis = arpeggia::align_seqs(&args.reference, &args.query, &args.alignment.options())?;
    for warning in analysis.warnings {
        tracing::warn!("{warning}");
    }
    let a = arpeggia::SeqAlignment {
        reference_name: args.reference_name.clone(),
        query_name: args.query_name.clone(),
        ..analysis.value
    };
    if args.json {
        return print_json(&a);
    }
    println!("{}", args.display.format(&a)?);
    Ok(())
}
