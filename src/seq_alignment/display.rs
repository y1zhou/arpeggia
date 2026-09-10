//! Shared alignment renderer for Rust, CLI, and Python.
use super::SeqAlignment;
use crate::{ArpeggiaError, ArpeggiaResult};
use clap::builder::styling::{AnsiColor, Style};
use std::fmt::Write;
use std::io::IsTerminal;
use unicode_width::UnicodeWidthStr;

/// Color policy for human-readable sequence alignments.
pub use clap::ColorChoice as AlignmentColor;

pub(crate) fn color_enabled(color: AlignmentColor, terminal: bool) -> bool {
    match color {
        AlignmentColor::Always => true,
        AlignmentColor::Never => false,
        AlignmentColor::Auto => {
            terminal && !anstyle_query::no_color() && anstyle_query::term_supports_color()
        }
    }
}

impl SeqAlignment {
    /// Render statistics and the complete alignment within a terminal width (80
    /// if unavailable) or explicit total line width. Rulers label every tenth
    /// input residue. Widths too small for labels and one residue fail.
    pub fn format(
        &self,
        width: Option<usize>,
        color: AlignmentColor,
        rulers: bool,
    ) -> ArpeggiaResult<String> {
        let width = width.unwrap_or_else(|| {
            terminal_size::terminal_size().map_or(80, |(w, _)| usize::from(w.0))
        });
        self.render(
            width,
            color_enabled(color, std::io::stdout().is_terminal()),
            rulers,
        )
    }

    pub(crate) fn render(&self, width: usize, color: bool, rulers: bool) -> ArpeggiaResult<String> {
        let digits = self.reference.len().max(self.query.len()).to_string().len();
        let names = [
            display_name(&self.reference_name),
            display_name(&self.query_name),
        ];
        let label_width = names.iter().map(|n| n.width()).max().unwrap_or(0);
        let prefix = label_width + 2 + digits; // label, start coordinate, spaces
        let overhead = prefix + 1 + digits;
        let block_width = width
            .checked_sub(overhead)
            .filter(|v| *v > 0)
            .ok_or_else(|| {
                ArpeggiaError::InvalidArgument(format!(
                    "alignment width must be at least {} columns",
                    overhead + 1
                ))
            })?;
        let ratio =
            |value: Option<f64>| value.map_or_else(|| "undefined".into(), |v| format!("{v:.4}"));
        let summary = format!(
            "{} alignment; {} (gap open {}, extend {})\nScore: {}; edit distance: {}\nMatches: {}; mismatches: {}; paired residues: {}\nGaps: {} residues in {} runs\nIdentity (alignment / shorter): {} / {:.4}\nCoverage (alignment / shorter): {} / {:.4}\nLengths (alignment / shorter): {} / {}",
            self.mode,
            self.matrix,
            self.gap_open,
            self.gap_extend,
            self.score,
            self.edit_distance,
            self.matches,
            self.mismatches,
            self.paired_residues,
            self.gap_residues,
            self.gap_runs,
            ratio(self.identity_alignment),
            self.identity_shorter,
            ratio(self.coverage_alignment),
            self.coverage_shorter,
            self.alignment_length,
            self.shorter_length
        );
        let mut output = String::new();
        // ASCII summaries are unstyled; wrapping also bounds long numeric values.
        for line in summary
            .lines()
            .chain((self.alignment_length == 0).then_some("No positive-scoring alignment"))
        {
            for chunk in line.as_bytes().chunks(width) {
                output.push_str(std::str::from_utf8(chunk).expect("ASCII summary"));
                output.push('\n');
            }
        }
        let leading = self.reference_span.0.max(self.query_span.0);
        let trailing = (self.reference.len() - self.reference_span.1)
            .max(self.query.len() - self.query_span.1);
        let padded = |sequence: &str, span: (usize, usize), aligned: &str| {
            format!(
                "{}{}{}{}{}",
                " ".repeat(leading - span.0),
                &sequence[..span.0],
                aligned,
                &sequence[span.1..],
                " ".repeat(trailing - (sequence.len() - span.1))
            )
        };
        let first = padded(
            &self.reference,
            self.reference_span,
            &self.aligned_reference,
        );
        let second = padded(&self.query, self.query_span, &self.aligned_query);
        // Dot marks display-only clipping, never a public alignment operation.
        let operations = format!(
            "{}{}{}",
            ".".repeat(leading),
            self.operations,
            ".".repeat(trailing)
        );
        let mut positions = [0usize; 2];
        for start in (0..operations.len()).step_by(block_width) {
            let end = start.saturating_add(block_width).min(operations.len());
            let ops = &operations.as_bytes()[start..end];
            output.push('\n');
            for (row, (label, sequence)) in
                [(names[0].as_str(), &first), (names[1].as_str(), &second)]
                    .into_iter()
                    .enumerate()
            {
                let cells = &sequence.as_bytes()[start..end];
                let mut ruler = vec![b' '; prefix + cells.len()];
                let mut begin = None;
                for (column, &cell) in cells.iter().enumerate() {
                    if cell != b'-' && cell != b' ' {
                        positions[row] += 1;
                        begin.get_or_insert(positions[row]);
                        if rulers && positions[row].is_multiple_of(10) {
                            let number = positions[row].to_string();
                            let stop = prefix + column + 1;
                            // Left gutter accommodates ticks at wrap boundaries.
                            ruler[stop - number.len()..stop].copy_from_slice(number.as_bytes());
                        }
                    }
                }
                if rulers {
                    output.push_str(std::str::from_utf8(&ruler).expect("ASCII ruler").trim_end());
                    output.push('\n');
                }
                let first_number = begin.map_or_else(String::new, |n| n.to_string());
                let last_number = begin.map_or_else(String::new, |_| positions[row].to_string());
                output.push_str(label);
                output.push_str(&" ".repeat(label_width - label.width()));
                output.push_str(&format!(" {first_number:>digits$} "));
                for (&cell, &op) in cells.iter().zip(ops) {
                    paint(&mut output, cell, op, color);
                }
                output.push_str(&format!(" {last_number:>digits$}\n"));
            }
            output.push_str(&" ".repeat(prefix));
            for &op in ops {
                paint(&mut output, if op == b'.' { b' ' } else { op }, op, color);
            }
            output.push('\n');
        }
        Ok(output.trim_end_matches('\n').into())
    }
}

// Labels are metadata; escape controls only in the human-readable rendering.
fn display_name(name: &str) -> String {
    let mut label = String::new();
    for c in name.chars() {
        if c.is_control() {
            label.extend(c.escape_default());
        } else {
            label.push(c);
        }
    }
    label
}

fn paint(output: &mut String, cell: u8, operation: u8, color: bool) {
    let style = match operation {
        b'+' => AnsiColor::Green.on_default(),
        b'-' => AnsiColor::Red.on_default(),
        b':' => AnsiColor::Blue.on_default(),
        b'x' => AnsiColor::Yellow.on_default(),
        b'.' => AnsiColor::BrightBlack.on_default(),
        _ => Style::new(),
    };
    if color && cell != b' ' && !style.is_plain() {
        write!(
            output,
            "{style}{}{reset}",
            cell as char,
            reset = style.render_reset()
        )
        .expect("writing to a String cannot fail");
    } else {
        output.push(cell as char);
    }
}

#[cfg(test)]
mod tests {
    use crate::{AlignmentMode, SeqAlignOptions, align_seqs};
    use unicode_width::UnicodeWidthStr;
    #[test]
    fn custom_names_preserve_unicode_width_and_escape_controls() {
        let mut a = align_seqs("ACDEFGHIKLMN", "ACDFGHIKLMN", &SeqAlignOptions::default())
            .unwrap()
            .value;
        a.reference_name = "野生型".into();
        a.query_name = "Mutant\nβ".into();
        let text = a.render(40, false, true).unwrap();
        assert!(text.contains("野生型"));
        assert!(text.contains("Mutant\\nβ"));
        assert!(text.lines().all(|line| line.width() <= 40));
        let lines: Vec<_> = text.lines().collect();
        let row = lines.iter().position(|l| l.starts_with("野生型")).unwrap();
        let tenth = lines[row][..lines[row].find('L').unwrap()].width();
        assert_eq!(&lines[row - 1][tenth - 1..=tenth], "10");
        assert_eq!(a.query_name, "Mutant\nβ");
        assert!(a.render(10, false, true).is_err());
    }
    #[test]
    fn wraps_and_numbers_residues_across_gaps() {
        let a = align_seqs(
            "ACDEFGHIKLMNPQRSTVWY",
            "ACDFGHIKLMNPQRSTVWY",
            &SeqAlignOptions::default(),
        )
        .unwrap()
        .value;
        let text = a.render(80, false, true).unwrap();
        let lines: Vec<_> = text.lines().collect();
        for label in ["Reference", "Query"] {
            let row = lines
                .iter()
                .position(|line| line.starts_with(label))
                .unwrap();
            let tenth = lines[row]
                .find(if label == "Reference" { 'L' } else { 'M' })
                .unwrap();
            assert_eq!(&lines[row - 1][tenth - 1..=tenth], "10");
        }
        for width in [18, 23, 40, 80] {
            let text = a.render(width, false, true).unwrap();
            assert!(text.lines().all(|line| line.len() <= width));
            let recovered: String = text
                .lines()
                .filter(|l| l.starts_with("Reference"))
                .map(|l| l.split_whitespace().nth(2).unwrap())
                .collect();
            assert_eq!(recovered, a.aligned_reference);
        }
        assert!(a.render(1, false, true).is_err());
    }
    #[test]
    fn colors_operations_and_clipping_without_changing_data() {
        let a = align_seqs(
            "GGACDEFGHIKGG",
            "ACDEYGHIK",
            &SeqAlignOptions {
                mode: AlignmentMode::SemiGlobal,
                ..Default::default()
            },
        )
        .unwrap()
        .value;
        let colored = a.render(80, true, false).unwrap();
        assert!(colored.contains("\x1b[90mG\x1b[0m"));
        for cell in ['F', 'Y', ':'] {
            assert!(colored.contains(&format!("\x1b[34m{cell}\x1b[0m")));
        }
        assert!(!a.render(80, false, false).unwrap().contains('\x1b'));
        for (reference, query, marker, escape) in [
            ("ACDEFGHIK", "ACDQQQEFGHIK", '+', "\x1b[32m"),
            ("ACDQQQEFGHIK", "ACDEFGHIK", '-', "\x1b[31m"),
            ("A", "G", 'x', "\x1b[33m"),
        ] {
            let a = align_seqs(reference, query, &SeqAlignOptions::default())
                .unwrap()
                .value;
            assert!(
                a.render(80, true, false)
                    .unwrap()
                    .contains(&format!("{escape}{marker}\x1b[0m"))
            );
        }
        let empty = align_seqs(
            "AAAA",
            "WWWW",
            &SeqAlignOptions {
                mode: AlignmentMode::Local,
                ..Default::default()
            },
        )
        .unwrap()
        .value;
        let text = empty.render(80, true, false).unwrap();
        assert!(text.contains("No positive-scoring alignment"));
        assert!(text.contains("\x1b[90mW\x1b[0m"));
    }
}
