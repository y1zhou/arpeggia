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
        let rows: Vec<_> = [(&self.reference_name, first), (&self.query_name, second)]
            .into_iter()
            .enumerate()
            .map(|(index, (name, cells))| {
                let mut position = 0usize;
                let labels = cells
                    .bytes()
                    .map(|cell| {
                        if cell == b'-' || cell == b' ' {
                            None
                        } else {
                            position += 1;
                            Some(position)
                        }
                    })
                    .collect();
                DisplayRow {
                    name: name.clone(),
                    right_name: String::new(),
                    cells: cells.into_bytes(),
                    positions: labels,
                    operations: operations.as_bytes().to_vec(),
                    imputed: Vec::new(),
                    show_operations: index == 1,
                }
            })
            .collect();
        let body = render_rows(&rows, width, color, rulers, &[])?;
        let mut output = wrap_summary(&summary, width);
        if self.alignment_length == 0 {
            output.push_str(&wrap_summary("No positive-scoring alignment", width));
        }
        output.push_str(&body);
        Ok(output.trim_end_matches('\n').into())
    }
}

/// Internal display data. Cells, positions and operations have equal lengths;
/// imputed is either empty or one flag per cell. Source coordinates are one-based.
pub(crate) struct DisplayRow {
    pub(crate) name: String,
    pub(crate) right_name: String,
    pub(crate) cells: Vec<u8>,
    pub(crate) positions: Vec<Option<usize>>,
    pub(crate) operations: Vec<u8>,
    pub(crate) imputed: Vec<bool>,
    pub(crate) show_operations: bool,
}

pub(crate) fn render_rows(
    rows: &[DisplayRow],
    width: usize,
    color: bool,
    rulers: bool,
    regions: &[u8],
) -> ArpeggiaResult<String> {
    let names: Vec<_> = rows.iter().map(|r| display_name(&r.name)).collect();
    let right_names: Vec<_> = rows.iter().map(|r| display_name(&r.right_name)).collect();
    let label_width = names.iter().map(|n| n.width()).max().unwrap_or(0);
    let right_width = right_names.iter().map(|n| n.width()).max().unwrap_or(0);
    let digits = rows
        .iter()
        .flat_map(|r| r.positions.iter().flatten())
        .map(|p| p.ilog10() as usize + 1)
        .max()
        .unwrap_or(1);
    let prefix = label_width + 2 + digits;
    let overhead = prefix + 1 + digits + if right_width == 0 { 0 } else { right_width + 1 };
    let block_width = width
        .checked_sub(overhead)
        .filter(|v| *v > 0)
        .ok_or_else(|| {
            ArpeggiaError::InvalidArgument(format!(
                "alignment width must be at least {} columns",
                overhead + 1
            ))
        })?;
    let length = rows.first().map_or(0, |r| r.cells.len());
    debug_assert!(rows.iter().all(|r| r.cells.len() == length
        && r.operations.len() == length
        && r.positions.len() == length
        && (r.imputed.is_empty() || r.imputed.len() == length)));
    debug_assert!(regions.is_empty() || regions.len() == length);
    let mut output = String::new();
    for start in (0..length).step_by(block_width) {
        let end = start.saturating_add(block_width).min(length);
        output.push('\n');
        // Antibody callers supply one reference-derived band map for the whole
        // block. Ordinary pairwise alignments pass no region annotations.
        if !regions.is_empty() {
            output.push_str(&" ".repeat(prefix));
            let mut labels = vec![b' '; end - start];
            let mut offset = 0;
            for run in regions[start..end].chunk_by(|a, b| a == b) {
                if run[0] > 0 {
                    labels[offset..offset + run.len()].fill(b'-');
                    let name = format!("CDR{}", run[0]);
                    if run.len() >= name.len() {
                        let at = offset + (run.len() - name.len()) / 2;
                        labels[at..at + name.len()].copy_from_slice(name.as_bytes());
                    } else {
                        labels[offset..offset + run.len()].fill(b'0' + run[0]);
                    }
                }
                offset += run.len();
            }
            for (column, cell) in labels.into_iter().enumerate() {
                paint(
                    &mut output,
                    cell,
                    b' ',
                    regions[start + column],
                    false,
                    color,
                );
            }
            output.push('\n');
        }
        for (index, row) in rows.iter().enumerate() {
            let cells = &row.cells[start..end];
            let coordinates = &row.positions[start..end];
            if rulers {
                let mut ruler = vec![b' '; prefix + cells.len()];
                for (column, label) in coordinates.iter().enumerate() {
                    if let Some(position) = label.filter(|p| p.is_multiple_of(10)) {
                        let text = position.to_string();
                        let stop = prefix + column + 1;
                        ruler[stop - text.len()..stop].copy_from_slice(text.as_bytes());
                    }
                }
                let visible = if regions.is_empty() {
                    ruler.iter().rposition(|b| *b != b' ').map_or(0, |i| i + 1)
                } else {
                    ruler.len()
                };
                for (column, &cell) in ruler[..visible].iter().enumerate() {
                    let region = column
                        .checked_sub(prefix)
                        .and_then(|i| regions.get(start + i))
                        .copied()
                        .unwrap_or(0);
                    paint(&mut output, cell, b' ', region, false, color);
                }
                output.push('\n');
            }
            let begin = coordinates
                .iter()
                .flatten()
                .next()
                .map_or_else(String::new, |p| p.to_string());
            let last = coordinates
                .iter()
                .flatten()
                .next_back()
                .map_or_else(String::new, |p| p.to_string());
            output.push_str(&names[index]);
            output.push_str(&" ".repeat(label_width - names[index].width()));
            write!(output, " {begin:>digits$} ").expect("String write");
            for (column, &cell) in cells.iter().enumerate() {
                paint(
                    &mut output,
                    cell,
                    row.operations[start + column],
                    regions.get(start + column).copied().unwrap_or(0),
                    row.imputed.get(start + column).copied().unwrap_or(false),
                    color,
                );
            }
            write!(output, " {last:>digits$}").expect("String write");
            if right_width > 0 {
                write!(output, " {}", right_names[index]).expect("String write");
            }
            output.push('\n');
            if row.show_operations {
                output.push_str(&" ".repeat(prefix));
                for column in start..end {
                    let op = row.operations[column];
                    paint(
                        &mut output,
                        if op == b'.' { b' ' } else { op },
                        op,
                        0,
                        false,
                        color,
                    );
                }
                output.push('\n');
            }
        }
    }
    Ok(output)
}

pub(crate) fn wrap_summary(summary: &str, width: usize) -> String {
    use unicode_width::UnicodeWidthChar;
    let mut output = String::new();
    for line in summary.lines() {
        let mut used = 0;
        for c in line.chars() {
            let size = c.width().unwrap_or(0);
            if used + size > width {
                output.push('\n');
                used = 0;
            }
            output.push(c);
            used += size;
        }
        output.push('\n');
    }
    output
}

// Labels are metadata; escape controls only in the human-readable rendering.
pub(crate) fn display_name(name: &str) -> String {
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

pub(crate) fn region_style(region: u8) -> Style {
    let background = match region {
        1 => AnsiColor::White,
        2 => AnsiColor::BrightMagenta,
        3 => AnsiColor::BrightCyan,
        _ => return Style::new(),
    };
    AnsiColor::Black.on(background)
}

fn paint(output: &mut String, cell: u8, operation: u8, region: u8, imputed: bool, color: bool) {
    let mut style = region_style(region);
    let foreground = match operation {
        b'+' => Some(AnsiColor::Green),
        b'-' => Some(AnsiColor::Red),
        b':' => Some(AnsiColor::Blue),
        b'x' => Some(AnsiColor::Yellow),
        b'.' => Some(AnsiColor::BrightBlack),
        _ => None,
    };
    if let Some(foreground) = foreground {
        style = style.fg_color(Some(foreground.into()));
    }
    if imputed {
        style = style.invert();
    }
    if color && (cell != b' ' || region > 0) && !style.is_plain() {
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
