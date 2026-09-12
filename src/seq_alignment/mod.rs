//! Exact pairwise protein sequence correspondence.

use crate::{Analysis, AnalysisWarning, ArpeggiaError, ArpeggiaResult, WarningCode};
use serde::Serialize;
pub(crate) mod display;
mod matrices;
pub use display::AlignmentColor;
#[cfg(feature = "python")]
pub(crate) use display::color_enabled;

/// Pairwise sequence-alignment objective.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum AlignmentMode {
    /// Consume both sequences, penalizing terminal gaps.
    #[default]
    Global,
    /// Select the highest-scoring subsequences.
    Local,
    /// Consume all of the second sequence with free first-sequence overhangs.
    SemiGlobal,
}

impl AlignmentMode {
    fn backend(self) -> hyalite::Mode {
        match self {
            Self::Global => hyalite::Mode::Nw,
            Self::Local => hyalite::Mode::Sw,
            Self::SemiGlobal => hyalite::Mode::Shw,
        }
    }
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Global => "global",
            Self::Local => "local",
            Self::SemiGlobal => "semi-global",
        }
    }
}

/// BLOSUM62 alignment settings. Gap costs have at most two decimal places.
#[derive(Clone, Debug, Serialize)]
pub struct SeqAlignOptions {
    /// Alignment objective.
    pub mode: AlignmentMode,
    /// Positive cost for the first residue of a gap.
    pub gap_open: f64,
    /// Positive cost for each subsequent gap residue, no larger than opening.
    pub gap_extend: f64,
}
impl Default for SeqAlignOptions {
    fn default() -> Self {
        Self {
            mode: AlignmentMode::Global,
            gap_open: 10.0,
            gap_extend: 0.5,
        }
    }
}

/// One exact protein alignment with zero-based, half-open spans.
///
/// Gapped strings contain the scored alignment; terminal clipping is represented
/// by the spans. Identity uses original symbols, including distinct U/C and O/K.
#[derive(Clone, Debug, Serialize)]
#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object, module = "arpeggia")
)]
pub struct SeqAlignment {
    /// Display name of the reference sequence.
    pub reference_name: String,
    /// Display name of the query sequence.
    pub query_name: String,
    /// Normalized first sequence (reference).
    pub reference: String,
    /// Normalized second sequence.
    pub query: String,
    /// Global, local, or semi-global objective.
    pub mode: String,
    /// Substitution matrix name.
    pub matrix: String,
    /// Gap opening cost.
    pub gap_open: f64,
    /// Gap extension cost.
    pub gap_extend: f64,
    /// Optimal unscaled protein alignment score.
    pub score: f64,
    /// Aligned first-sequence span, excluding clipping.
    pub reference_span: (usize, usize),
    /// Aligned second-sequence span, excluding clipping.
    pub query_span: (usize, usize),
    /// Scored reference alignment, with `-` for gaps and without clipped tails.
    pub aligned_reference: String,
    /// Scored query alignment, with `-` for gaps and without clipped tails.
    pub aligned_query: String,
    /// Reference-to-query operations: space (match), + (insertion), - (deletion),
    /// : (positive-score substitution), x (other substitution). All are ASCII.
    pub operations: String,
    /// Number of columns, including gaps.
    pub alignment_length: usize,
    /// Length of the shorter complete input.
    pub shorter_length: usize,
    /// Identical nongap pairs.
    pub matches: usize,
    /// Nonidentical nongap pairs.
    pub mismatches: usize,
    /// All nongap pairs.
    pub paired_residues: usize,
    /// Number of residues opposite gaps.
    pub gap_residues: usize,
    /// Number of contiguous gaps, counting directions separately.
    pub gap_runs: usize,
    /// Matches / alignment length; None for an empty local alignment.
    pub identity_alignment: Option<f64>,
    /// Matches / shorter full input length.
    pub identity_shorter: f64,
    /// Nongap pairs / alignment length; None for an empty local alignment.
    pub coverage_alignment: Option<f64>,
    /// Nongap pairs / shorter full input length.
    pub coverage_shorter: f64,
    /// Full-input unit-cost Levenshtein distance, independent of this traceback.
    pub edit_distance: usize,
}

impl SeqAlignment {
    // Normalized ASCII strings have equal aligned lengths. Recover correspondence
    // only when needed, without storing a second representation of the alignment.
    pub(crate) fn columns(&self) -> impl Iterator<Item = (Option<usize>, Option<usize>)> + '_ {
        let (mut i, mut j) = (self.reference_span.0, self.query_span.0);
        self.aligned_reference
            .bytes()
            .zip(self.aligned_query.bytes())
            .map(move |(a, b)| {
                let left = (a != b'-').then(|| {
                    let index = i;
                    i += 1;
                    index
                });
                let right = (b != b'-').then(|| {
                    let index = j;
                    j += 1;
                    index
                });
                (left, right)
            })
    }
}

pub(crate) fn scaled_cost(value: f64) -> ArpeggiaResult<i32> {
    let scaled = value * 100.0;
    if !value.is_finite()
        || scaled.round() < 1.0
        || scaled > i32::MAX as f64
        || (scaled - scaled.round()).abs() > f64::EPSILON * scaled.abs().max(1.0) * 4.0
    {
        return Err(ArpeggiaError::InvalidArgument(
            "gap costs must be positive, fit integer scoring, and have at most two decimal places"
                .into(),
        ));
    }
    Ok(scaled.round() as i32)
}

pub(crate) fn scoring(options: &SeqAlignOptions) -> ArpeggiaResult<hyalite::Scoring> {
    let open = scaled_cost(options.gap_open)?;
    let extend = scaled_cost(options.gap_extend)?;
    hyalite::Scoring::new(
        23,
        matrices::BLOSUM62.iter().map(|v| v * 100).collect(),
        open,
        extend,
    )
    .map_err(|e| ArpeggiaError::InvalidArgument(e.to_string()))
}

pub(crate) fn encode(sequence: &str) -> ArpeggiaResult<Vec<u8>> {
    if sequence.is_empty() {
        return Err(ArpeggiaError::InvalidArgument(
            "sequences must not be empty".into(),
        ));
    }
    sequence.bytes().map(|b| {
        encode_letter(b)
            .ok_or_else(|| ArpeggiaError::InvalidArgument("sequences must contain only amino-acid letters (including B/Z/X/U/O), without gaps, whitespace, or stops".into()))
    }).collect()
}

fn encode_letter(letter: u8) -> Option<u8> {
    let letter = match letter.to_ascii_uppercase() {
        b'U' => b'C',
        b'O' => b'K',
        c => c,
    };
    matrices::BLOSUM62_ALPHABET
        .iter()
        .position(|a| *a == letter)
        .map(|i| i as u8)
}

// Both callers supply validated amino-acid symbols or alignment gap cells.
pub(crate) fn operation(reference: u8, query: u8) -> u8 {
    if reference == query {
        return b' ';
    }
    if reference == b'-' {
        return b'+';
    }
    if query == b'-' {
        return b'-';
    }
    let row = usize::from(encode_letter(reference).expect("validated reference residue"));
    let column = usize::from(encode_letter(query).expect("validated query residue"));
    if matrices::BLOSUM62[row * 23 + column] > 0 {
        b':'
    } else {
        b'x'
    }
}

pub(crate) fn backend_error(error: hyalite::Error) -> ArpeggiaError {
    ArpeggiaError::Calculation(format!("sequence alignment: {error}"))
}

/// Align two unaligned protein strings with BLOSUM62 and affine gaps.
///
/// Returns one deterministic optimum. No positive local match is a successful
/// empty alignment. U/O score as C/K with a diagnostic but retain their identity.
pub fn align_seqs(
    reference: &str,
    query: &str,
    options: &SeqAlignOptions,
) -> ArpeggiaResult<Analysis<SeqAlignment>> {
    let scoring = scoring(options)?;
    let first = encode(reference)?;
    let second = encode(query)?;
    let reference = reference.to_ascii_uppercase();
    let query = query.to_ascii_uppercase();
    // Limit the full-matrix working set; the backend recomputes checkpoints above
    // this budget without changing the optimum. Retry budget exhaustion with the
    // exact checkpoint requirement so this is not a sequence-length ceiling.
    let alignment = match hyalite::align(
        &first,
        &second,
        &scoring,
        options.mode.backend(),
        64 * 1024 * 1024,
    ) {
        Err(hyalite::Error::TracebackBudgetExceeded { needed_bytes, .. }) => {
            let budget = usize::try_from(needed_bytes).map_err(|_| {
                ArpeggiaError::Calculation(
                    "alignment memory requirement exceeds address space".into(),
                )
            })?;
            hyalite::align(&first, &second, &scoring, options.mode.backend(), budget)
        }
        result => result,
    }
    .map_err(backend_error)?;
    let mut aligned_reference = String::with_capacity(alignment.ops.len());
    let mut aligned_query = String::with_capacity(alignment.ops.len());
    let mut operations = String::with_capacity(alignment.ops.len());
    let (mut i, mut j) = (alignment.query_start, alignment.target_start);
    let (mut matches, mut paired, mut gap_runs) = (0, 0, 0);
    let mut last_gap = None;
    for op in &alignment.ops {
        let gap = match op {
            hyalite::AlignOp::Ins => Some(0),
            hyalite::AlignOp::Del => Some(1),
            _ => None,
        };
        if gap.is_some() && gap != last_gap {
            gap_runs += 1;
        }
        last_gap = gap;
        match op {
            hyalite::AlignOp::Match | hyalite::AlignOp::Mismatch => {
                let same = reference.as_bytes()[i] == query.as_bytes()[j];
                aligned_reference.push(reference.as_bytes()[i] as char);
                aligned_query.push(query.as_bytes()[j] as char);
                operations.push(operation(reference.as_bytes()[i], query.as_bytes()[j]) as char);
                paired += 1;
                matches += usize::from(same);
                i += 1;
                j += 1;
            }
            hyalite::AlignOp::Ins => {
                aligned_reference.push(reference.as_bytes()[i] as char);
                aligned_query.push('-');
                operations.push('-');
                i += 1;
            }
            hyalite::AlignOp::Del => {
                aligned_reference.push('-');
                aligned_query.push(query.as_bytes()[j] as char);
                operations.push('+');
                j += 1;
            }
        }
    }
    let length = operations.len();
    let shorter = reference.len().min(query.len());
    let mut warnings = Vec::new();
    if reference
        .bytes()
        .chain(query.bytes())
        .any(|b| matches!(b, b'U' | b'O'))
    {
        warnings.push(AnalysisWarning::new(
            WarningCode::SequenceScoringAlias,
            "U/O score as C/K in BLOSUM62; identities retain the original symbols",
        ));
    }
    let edit_distance = levenshtein(reference.as_bytes(), query.as_bytes());
    Ok(Analysis::new(
        SeqAlignment {
            reference_name: "Reference".into(),
            query_name: "Query".into(),
            reference,
            query,
            mode: options.mode.name().into(),
            matrix: "BLOSUM62".into(),
            gap_open: options.gap_open,
            gap_extend: options.gap_extend,
            score: f64::from(alignment.score) / 100.0,
            reference_span: (alignment.query_start, alignment.query_end),
            query_span: (alignment.target_start, alignment.target_end),
            aligned_reference,
            aligned_query,
            operations,
            alignment_length: length,
            shorter_length: shorter,
            matches,
            mismatches: paired - matches,
            paired_residues: paired,
            gap_residues: length - paired,
            gap_runs,
            identity_alignment: (length > 0).then(|| matches as f64 / length as f64),
            identity_shorter: matches as f64 / shorter as f64,
            coverage_alignment: (length > 0).then(|| paired as f64 / length as f64),
            coverage_shorter: paired as f64 / shorter as f64,
            edit_distance,
        },
        warnings,
    ))
}

// Unit-cost edit distance is a separate objective from protein alignment. One
// rolling row avoids adding a second alignment dependency or storing traceback.
fn levenshtein(a: &[u8], b: &[u8]) -> usize {
    let (a, b) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    let mut row: Vec<_> = (0..=b.len()).collect();
    for (i, x) in a.iter().enumerate() {
        let mut diagonal = row[0];
        row[0] = i + 1;
        for (j, y) in b.iter().enumerate() {
            let old = row[j + 1];
            row[j + 1] = (diagonal + usize::from(x != y))
                .min(row[j] + 1)
                .min(old + 1);
            diagonal = old;
        }
    }
    row[b.len()]
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn global_statistics_and_semiglobal_orientation() {
        let a = align_seqs("acdefghik", "ACDEYGHIK", &SeqAlignOptions::default())
            .unwrap()
            .value;
        assert_eq!(a.score, 50.0); // Biopython 1.86, BLOSUM62, open=10, extend=.5.
        assert_eq!((a.matches, a.mismatches, a.edit_distance), (8, 1, 1));
        assert_eq!(a.identity_alignment, Some(8.0 / 9.0));
        let opts = SeqAlignOptions {
            mode: AlignmentMode::SemiGlobal,
            ..Default::default()
        };
        let a = align_seqs("GGACDEFGHIKGG", "ACDEFGHIK", &opts)
            .unwrap()
            .value;
        assert_eq!(a.reference_span, (2, 11));
        assert_eq!(a.query_span, (0, 9));
        assert_eq!(a.edit_distance, 4);
        assert_eq!(a.coverage_shorter, 1.0);
        let reverse = align_seqs("ACDEFGHIK", "GGACDEFGHIKGG", &opts)
            .unwrap()
            .value;
        assert_eq!(reverse.query_span, (0, 13));
        assert!(reverse.gap_residues >= 4);
        assert!(reverse.score < a.score);
    }
    #[test]
    fn gaps_clipping_and_original_symbol_identity() {
        let a = align_seqs("ACDEFGHIK", "ACDQQQEFGHIK", &SeqAlignOptions::default())
            .unwrap()
            .value;
        assert_eq!((a.gap_residues, a.gap_runs, a.edit_distance), (3, 1, 3));
        assert_eq!(a.coverage_alignment, Some(9.0 / 12.0));
        assert_eq!(a.coverage_shorter, 1.0);
        assert_eq!(a.aligned_reference, "ACD---EFGHIK");
        assert_eq!(a.aligned_query, "ACDQQQEFGHIK");
        assert_eq!(a.operations, "   +++      ");
        assert_eq!(a.columns().nth(6), Some((Some(3), Some(6))));
        let reversed = align_seqs("ACDQQQEFGHIK", "ACDEFGHIK", &SeqAlignOptions::default())
            .unwrap()
            .value;
        assert_eq!(reversed.operations, "   ---      ");
        let a = align_seqs(
            "AAAA",
            "WWWW",
            &SeqAlignOptions {
                mode: AlignmentMode::Local,
                ..Default::default()
            },
        )
        .unwrap()
        .value;
        assert_eq!(a.score, 0.0);
        assert!(
            a.aligned_reference.is_empty() && a.aligned_query.is_empty() && a.operations.is_empty()
        );
        assert_eq!(a.identity_alignment, None);
        assert_eq!(a.coverage_alignment, None);
        assert_eq!(a.identity_shorter, 0.0);
        assert_eq!(a.edit_distance, 4);
        let a = align_seqs("UO", "CK", &SeqAlignOptions::default()).unwrap();
        assert_eq!(a.value.matches, 0);
        assert_eq!(a.value.operations, "::");
        assert_eq!(a.value.edit_distance, 2);
        assert_eq!(a.warnings.len(), 1);
    }
    #[test]
    fn similarity_uses_positive_scores_without_changing_identity() {
        for (reference, query, operation) in [
            ("F", "Y", ":"),
            ("A", "G", "x"),
            ("A", "X", "x"),
            ("B", "D", ":"),
            ("Z", "E", ":"),
            ("U", "C", ":"),
            ("O", "K", ":"),
            ("X", "X", " "),
        ] {
            let a = align_seqs(reference, query, &SeqAlignOptions::default())
                .unwrap()
                .value;
            assert_eq!(a.operations, operation, "{reference}/{query}");
            let identical = usize::from(reference == query);
            assert_eq!((a.matches, a.mismatches), (identical, 1 - identical));
            assert_eq!(a.edit_distance, 1 - identical);
            assert_eq!(a.identity_alignment, Some(identical as f64));
            assert_eq!(
                (a.alignment_length, a.paired_residues, a.gap_residues),
                (1, 1, 0)
            );
        }
    }
    #[test]
    fn validate_alphabet_and_precise_costs() {
        for input in ["", "A-C", "A*C", "A C", "é", "AJ"] {
            assert!(align_seqs(input, "AC", &SeqAlignOptions::default()).is_err());
        }
        for value in [
            0.0,
            -1.0,
            0.001,
            1e-20,
            f64::MIN_POSITIVE,
            f64::from_bits(1),
            f64::NAN,
            f64::INFINITY,
            1e15,
        ] {
            assert!(scaled_cost(value).is_err());
        }
        assert_eq!(scaled_cost(0.01).unwrap(), 1);
        assert_eq!(scaled_cost(0.29).unwrap(), 29);
        assert!(
            scoring(&SeqAlignOptions {
                gap_open: 0.5,
                gap_extend: 1.0,
                ..Default::default()
            })
            .is_err()
        );
        assert_eq!(levenshtein(b"ABC", b""), 3);
    }
    #[test]
    fn tied_paths_and_checkpoint_traceback_are_stable() {
        let scoring = scoring(&SeqAlignOptions::default()).unwrap();
        let a = encode(&"A".repeat(80)).unwrap();
        let b = encode(&"A".repeat(40)).unwrap();
        let full = hyalite::align(&a, &b, &scoring, hyalite::Mode::Shw, usize::MAX).unwrap();
        let checkpoint = hyalite::align(&a, &b, &scoring, hyalite::Mode::Shw, 20000).unwrap();
        assert_eq!(full, checkpoint);
        assert_eq!(full.query_start, 0);
    }
}
