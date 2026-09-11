use super::align_seqs::{DisplayArgs, print_json};
use arpeggia::{
    ArpeggiaError, ArpeggiaResult, CdrDefinition, GermlineSpecies, NumberedAntibody,
    NumberingOptions, NumberingScheme,
};
use clap::{Args, Parser};

#[derive(Args, Debug, Clone)]
struct NumberingArgs {
    /// Numbering convention (default: imgt; chothia is an alias for martin)
    #[arg(long, value_enum)]
    scheme: Option<NumberingScheme>,
    /// Region convention; non-auto requires an explicit --scheme
    #[arg(long, value_enum, default_value = "auto")]
    cdr_definition: CdrDefinition,
    /// Restrict germline references, comma-separated (default: all bundled species)
    #[arg(long, value_enum, value_delimiter = ',')]
    species: Vec<GermlineSpecies>,
    /// Impute only supported missing beginnings of FR1 and ends of FR4
    #[arg(long)]
    impute: bool,
    /// Exact tied V reference ID to use for imputation (from --json)
    #[arg(long, requires = "impute")]
    v_reference: Option<String>,
    /// Exact tied J reference ID to use for imputation (from --json)
    #[arg(long, requires = "impute")]
    j_reference: Option<String>,
}

impl NumberingArgs {
    fn number(&self, sequence: &str, name: String) -> ArpeggiaResult<NumberedAntibody> {
        let antibody = arpeggia::number_antibody(
            sequence,
            &NumberingOptions {
                name,
                scheme: self.scheme,
                cdr_definition: self.cdr_definition,
                species: self.species.clone(),
            },
        )?;
        if self.impute {
            antibody.impute(self.v_reference.as_deref(), self.j_reference.as_deref())
        } else {
            Ok(antibody)
        }
    }
}

#[derive(Parser, Debug, Clone)]
#[command(
    after_help = "Example: arpeggia number-antibody SEQUENCE --name VHH --scheme imgt --species human,alpaca\nSupply an unaligned amino-acid string, not a FASTA file. One variable domain is required;\ntags and constant tails are allowed. Terminal truncations must retain IMGT 23–118 anchor coverage.\nV/J results are reference similarities, with exact ties retained; they do not infer ancestry.\nThe input appears above stitched V/J germlines. Input rulers preserve original coordinates; the germline\nruler counts V then J continuously, ignoring gaps. Imputed positions are blank; yellow marks imputation.\nConventions and examples: https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md"
)]
pub(crate) struct NumberArgs {
    /// One unaligned antibody amino-acid sequence
    sequence: String,
    /// Display name of the input sequence
    #[arg(long, default_value = "Seq001")]
    name: String,
    #[command(flatten)]
    numbering: NumberingArgs,
    #[command(flatten)]
    display: DisplayArgs,
    /// Emit plain structured results, including tied references and provenance
    #[arg(long)]
    json: bool,
}

#[derive(Parser, Debug, Clone)]
#[command(
    after_help = "Example: arpeggia align-antibodies SEQUENCE1 SEQUENCE2 --names WT,Mutant --reference-index 1\nInputs must share a numbering scheme and be all heavy or all light chains (K/L may mix).\n--names preserves empty entries: --names WT,,Mutant assigns Seq002 to the second input.\nThe selected reference is displayed first and defines CDR bands across all rows.\nRulers use each original input; imputed positions are blank. Germlines are hidden.\nJSON rows stay in input order.\nConventions and examples: https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md"
)]
pub(crate) struct AlignArgs {
    /// One or more unaligned antibody amino-acid sequences
    #[arg(required = true)]
    sequences: Vec<String>,
    /// Comma-separated names; missing/empty names become Seq001, Seq002, ...
    #[arg(long)]
    names: Option<String>,
    /// Zero-based input index to use as the comparison reference
    #[arg(long, default_value_t = 0)]
    reference_index: usize,
    #[command(flatten)]
    numbering: NumberingArgs,
    #[command(flatten)]
    display: DisplayArgs,
    /// Emit plain structured results in input order
    #[arg(long)]
    json: bool,
}

pub(crate) fn run_number(args: &NumberArgs) -> ArpeggiaResult<()> {
    let result = args.numbering.number(&args.sequence, args.name.clone())?;
    if args.json {
        print_json(&result)
    } else {
        println!(
            "{}",
            result.format(
                args.display.width,
                args.display.color,
                !args.display.no_rulers
            )?
        );
        Ok(())
    }
}

pub(crate) fn run_alignment(args: &AlignArgs) -> ArpeggiaResult<()> {
    let names: Vec<_> = args
        .names
        .as_deref()
        .map_or_else(Vec::new, |n| n.split(',').map(str::trim).collect());
    if names.len() > args.sequences.len() {
        return Err(ArpeggiaError::InvalidArgument(
            "--names has more entries than input sequences".into(),
        ));
    }
    if args.reference_index >= args.sequences.len() {
        return Err(ArpeggiaError::InvalidArgument(
            "--reference-index must address an input sequence".into(),
        ));
    }
    let antibodies = args
        .sequences
        .iter()
        .enumerate()
        .map(|(i, sequence)| {
            let name = names
                .get(i)
                .filter(|n| !n.is_empty())
                .map_or_else(|| format!("Seq{:03}", i + 1), |n| (*n).into());
            args.numbering.number(sequence, name)
        })
        .collect::<ArpeggiaResult<Vec<_>>>()?;
    let result = arpeggia::align_antibodies(antibodies, args.reference_index)?;
    if args.json {
        print_json(&result)
    } else {
        for antibody in &result.antibodies {
            for diagnostic in &antibody.diagnostics {
                tracing::warn!("{}: {diagnostic}", antibody.name);
            }
        }
        println!(
            "{}",
            result.format(
                args.display.width,
                args.display.color,
                !args.display.no_rulers,
                None
            )?
        );
        Ok(())
    }
}
