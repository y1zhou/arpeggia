use super::*;
use crate::AlignmentColor;
use crate::seq_alignment::display::{
    DisplayRow, color_enabled, display_name, imputed_style, region_style, render_rows,
    wrap_styled_summary,
};
use clap::builder::styling::Style;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::IsTerminal;

fn dimensions(width: Option<usize>, color: AlignmentColor) -> (usize, bool) {
    (
        width.unwrap_or_else(|| {
            terminal_size::terminal_size().map_or(80, |(w, _)| usize::from(w.0))
        }),
        color_enabled(color, std::io::stdout().is_terminal()),
    )
}

fn region(name: &str) -> u8 {
    match name {
        "CDR1" => 1,
        "CDR2" => 2,
        "CDR3" => 3,
        _ => 0,
    }
}

// Reference insertions/deletions inside a CDR must not break its vertical band.
fn cdr_bands(mut regions: Vec<u8>) -> Vec<u8> {
    for cdr in 1..=3 {
        if let Some(start) = regions.iter().position(|r| *r == cdr) {
            let end = regions.iter().rposition(|r| *r == cdr).unwrap();
            regions[start..=end].fill(cdr);
        }
    }
    regions
}

fn summary(intro: &str, imputed: usize, details: &str, width: usize, color: bool) -> String {
    let count = format!("imputed residues: {imputed}");
    wrap_styled_summary(
        [
            (intro, Style::new()),
            (count.as_str(), imputed_style()),
            ("\nCDR regions: ", Style::new()),
            ("1 = CDR1", region_style(1)),
            (", ", Style::new()),
            ("2 = CDR2", region_style(2)),
            (", ", Style::new()),
            ("3 = CDR3", region_style(3)),
            ("\n", Style::new()),
            (details, Style::new()),
        ],
        width,
        color,
    )
}

fn gene_name(reference: &GermlineReference) -> String {
    format!("{}*{}", reference.gene, reference.allele)
}

impl AntibodyAlignment {
    /// Format antibody comparisons with the selected reference first and germlines hidden.
    ///
    /// Width includes names and original input-position labels. A rendering-only
    /// reference override does not change the stored reference index or row order.
    pub fn format(
        &self,
        width: Option<usize>,
        color: AlignmentColor,
        rulers: bool,
        reference_index: Option<usize>,
    ) -> ArpeggiaResult<String> {
        let (width, color) = dimensions(width, color);
        self.render(width, color, rulers, reference_index)
    }

    pub(crate) fn render(
        &self,
        width: usize,
        color: bool,
        rulers: bool,
        reference_index: Option<usize>,
    ) -> ArpeggiaResult<String> {
        let reference = reference_index.unwrap_or(self.reference_index);
        let reference_row = self
            .aligned_sequences
            .get(reference)
            .ok_or_else(|| {
                ArpeggiaError::InvalidArgument(
                    "reference_index must address an input antibody".into(),
                )
            })?
            .as_bytes();
        // Alignment construction guarantees each row has numbered residues.
        let span = |row: &[u8]| {
            row.iter().position(|b| *b != b'-').unwrap()
                ..=row.iter().rposition(|b| *b != b'-').unwrap()
        };
        let reference_span = span(reference_row);
        let mut regions = Vec::new();
        let rows: Vec<_> = std::iter::once(reference)
            .chain((0..self.antibodies.len()).filter(|i| *i != reference))
            .map(|index| {
                let antibody = &self.antibodies[index];
                let by_position: HashMap<_, _> =
                    antibody.residues.iter().map(|r| (r.position, r)).collect();
                let cells = self.aligned_sequences[index].as_bytes().to_vec();
                let query_span = span(&cells);
                let operations = cells
                    .iter()
                    .enumerate()
                    .map(|(i, cell)| {
                        if index == reference {
                            b' '
                        } else if !reference_span.contains(&i) || !query_span.contains(&i) {
                            b'.'
                        } else {
                            crate::seq_alignment::operation(reference_row[i], *cell)
                        }
                    })
                    .collect();
                let positions = self
                    .positions
                    .iter()
                    .map(|p| {
                        by_position
                            .get(p)
                            .and_then(|r| r.input_index)
                            .map(|i| i + 1)
                    })
                    .collect();
                if index == reference {
                    regions = self
                        .positions
                        .iter()
                        .map(|p| by_position.get(p).map_or(0, |r| region(&r.region)))
                        .collect();
                }
                let imputed = self
                    .positions
                    .iter()
                    .map(|p| by_position.get(p).is_some_and(|r| r.input_index.is_none()))
                    .collect();
                DisplayRow {
                    name: antibody.name.clone(),
                    right_name: String::new(),
                    cells,
                    positions,
                    operations,
                    imputed,
                    show_operations: index != reference,
                }
            })
            .collect();
        let body = render_rows(&rows, width, color, rulers, &cdr_bands(regions))?;
        let selected = &self.antibodies[reference];
        let intro = format!(
            "{} antibodies; {} positions\nReference: {}; {} numbering; {} CDR definition\nTotal ",
            self.antibodies.len(),
            self.positions.len(),
            display_name(&selected.name),
            selected.scheme,
            selected.cdr_definition,
        );
        let imputed = self
            .antibodies
            .iter()
            .flat_map(|a| &a.residues)
            .filter(|r| r.input_index.is_none())
            .count();
        let mut output = summary(&intro, imputed, "", width, color);
        output.push_str(&body);
        Ok(output.trim_end_matches('\n').into())
    }
}

#[derive(Clone)]
struct Column {
    query: u8,
    reference: u8,
    operation: u8,
    input_index: Option<usize>,
    query_position: Option<NumberedPosition>,
    reference_position: Option<NumberedPosition>,
    imputed: bool,
    region: u8,
}

impl NumberedAntibody {
    /// Format the input above combined V/J germlines, with distinct CDR backgrounds.
    ///
    /// Gray junction gaps mark unavailable reference sequence. V/J scores remain
    /// separate similarities measured on the supplied input, including after imputation.
    /// Rulers use original input positions and a continuous V-then-J residue count.
    /// Imputed input positions are blank and their residues have yellow backgrounds.
    pub fn format(
        &self,
        width: Option<usize>,
        color: AlignmentColor,
        rulers: bool,
    ) -> ArpeggiaResult<String> {
        let (width, color) = dimensions(width, color);
        self.render(width, color, rulers)
    }

    pub(crate) fn render(&self, width: usize, color: bool, rulers: bool) -> ArpeggiaResult<String> {
        let by_input: HashMap<_, _> = self
            .residues
            .iter()
            .filter_map(|r| r.input_index.map(|i| (i, r)))
            .collect();
        let mut columns: Vec<_> = self
            .input_sequence
            .bytes()
            .enumerate()
            .map(|(i, query)| Column {
                query,
                reference: b'-',
                operation: b'.',
                input_index: Some(i),
                query_position: by_input.get(&i).map(|r| r.position),
                reference_position: None,
                imputed: false,
                region: by_input.get(&i).map_or(0, |r| region(&r.region)),
            })
            .collect();
        let mut before = vec![Vec::new(); columns.len() + 1];
        let scheme: NumberingScheme = clap::ValueEnum::from_str(&self.scheme, true)
            .map_err(ArpeggiaError::InvalidArgument)?;
        let chain: Chain = self
            .chain
            .parse()
            .map_err(|_| ArpeggiaError::InvalidArgument("unsupported antibody chain".into()))?;
        for matching in [&self.v_match, &self.j_match].into_iter().flatten() {
            let hit = &matching.hits[0];
            // Position labels are optional: an unrepresentable reference position
            // must not prevent rendering its valid pairwise sequence alignment.
            let reference_positions = germline::reference_positions(hit, scheme.backend(), chain)
                .unwrap_or_else(|_| vec![None; hit.alignment.reference.len()]);
            let alignment = &hit.alignment;
            let mut cursor = hit.query_input_start + alignment.query_span.0;
            for ((r, q), op) in alignment.columns().zip(alignment.operations.bytes()) {
                if let Some(q) = q {
                    let index = hit.query_input_start + q;
                    columns[index].reference =
                        r.map_or(b'-', |i| alignment.reference.as_bytes()[i]);
                    columns[index].reference_position = r.and_then(|i| reference_positions[i]);
                    columns[index].operation = op;
                    cursor = index + 1;
                } else if let Some(r) = r {
                    before[cursor].push(Column {
                        query: b'-',
                        reference: alignment.reference.as_bytes()[r],
                        operation: op,
                        input_index: None,
                        query_position: None,
                        reference_position: reference_positions[r],
                        imputed: false,
                        region: 0,
                    });
                }
            }
            // Unaligned reference ends use the same right-/left-aligned clipping
            // policy as SeqAlignment; their blank markers do not assert edits.
            for (r, q) in (0..alignment.reference_span.0)
                .map(|r| {
                    (
                        r,
                        r as isize + alignment.query_span.0 as isize
                            - alignment.reference_span.0 as isize,
                    )
                })
                .chain(
                    (alignment.reference_span.1..alignment.reference.len()).map(|r| {
                        (
                            r,
                            (alignment.query_span.1 + r - alignment.reference_span.1) as isize,
                        )
                    }),
                )
            {
                if q >= 0 && (q as usize) < alignment.query.len() {
                    let column = &mut columns[hit.query_input_start + q as usize];
                    column.reference = alignment.reference.as_bytes()[r];
                    column.reference_position = reference_positions[r];
                } else {
                    let at = hit.query_input_start + if q < 0 { 0 } else { alignment.query.len() };
                    before[at].push(Column {
                        query: b'-',
                        reference: alignment.reference.as_bytes()[r],
                        operation: b'.',
                        input_index: None,
                        query_position: None,
                        reference_position: reference_positions[r],
                        imputed: false,
                        region: 0,
                    });
                }
            }
        }
        let mut combined = Vec::new();
        for (i, column) in columns.into_iter().enumerate() {
            combined.append(&mut before[i]);
            combined.push(column);
        }
        combined.append(before.last_mut().unwrap());
        for residue in self.residues.iter().filter(|r| r.input_index.is_none()) {
            if let Some(column) = combined
                .iter_mut()
                .find(|c| c.query == b'-' && c.reference_position == Some(residue.position))
            {
                column.query = residue.amino_acid as u8;
                column.query_position = Some(residue.position);
                column.imputed = true;
                column.operation = crate::seq_alignment::operation(column.reference, column.query);
            } else {
                let at = combined
                    .iter()
                    .position(|c| {
                        c.query_position.is_some_and(|p| {
                            p.order(scheme.backend(), chain)
                                > residue.position.order(scheme.backend(), chain)
                        }) || c.input_index.is_some_and(|i| i >= self.domain_span.1)
                    })
                    .unwrap_or(combined.len());
                combined.insert(
                    at,
                    Column {
                        query: residue.amino_acid as u8,
                        reference: b'-',
                        operation: b'.',
                        input_index: None,
                        query_position: Some(residue.position),
                        reference_position: None,
                        imputed: true,
                        region: 0,
                    },
                );
            }
        }
        let first_reference = combined.iter().position(|c| c.reference != b'-');
        let last_reference = combined.iter().rposition(|c| c.reference != b'-');
        for (index, column) in combined.iter_mut().enumerate() {
            if first_reference.is_none_or(|first| index < first)
                || last_reference.is_none_or(|last| index > last)
            {
                column.reference = b' ';
                column.operation = b'.';
            }
        }
        let regions = cdr_bands(combined.iter().map(|c| c.region).collect());
        let name = |m: &Option<GermlineMatch>, absent: &str| {
            m.as_ref().map_or_else(
                || absent.into(),
                |m| {
                    let r = &m.hits[0].references[0];
                    gene_name(r)
                },
            )
        };
        let rows: Vec<_> = [false, true]
            .into_iter()
            .map(|germline| DisplayRow {
                name: if germline {
                    name(&self.v_match, "V unavailable")
                } else {
                    self.name.clone()
                },
                right_name: if germline {
                    name(&self.j_match, "J unavailable")
                } else {
                    String::new()
                },
                cells: combined
                    .iter()
                    .map(|c| if germline { c.reference } else { c.query })
                    .collect(),
                positions: combined
                    .iter()
                    .scan(0, |position, c| {
                        Some(if germline {
                            c.reference.is_ascii_alphabetic().then(|| {
                                *position += 1;
                                *position
                            })
                        } else {
                            c.input_index.map(|i| i + 1)
                        })
                    })
                    .collect(),
                operations: combined
                    .iter()
                    .map(|c| {
                        if germline {
                            match c.operation {
                                b'+' => b'-',
                                b'-' => b'+',
                                op => op,
                            }
                        } else if c.query_position.is_none() {
                            b'.'
                        } else {
                            b' '
                        }
                    })
                    .collect(),
                imputed: if germline {
                    Vec::new()
                } else {
                    combined.iter().map(|c| c.imputed).collect()
                },
                show_operations: germline,
            })
            .collect();
        let body = render_rows(&rows, width, color, rulers, &regions)?;
        let intro = format!(
            "Reference: {}; {} chain; {} numbering; {} CDR definition\nConfidence: {:.3}; matched profile positions: {}\nNumbered domain in supplied input: {}–{} (1-based)\n",
            display_name(&self.name),
            self.chain,
            self.scheme,
            self.cdr_definition,
            self.confidence,
            self.matched_profile_positions,
            self.domain_span.0 + 1,
            self.domain_span.1,
        );
        let imputed = self
            .residues
            .iter()
            .filter(|r| r.input_index.is_none())
            .count();
        let mut details = String::new();
        for (segment, matching) in [("V", &self.v_match), ("J", &self.j_match)] {
            if let Some(matching) = matching {
                let hit = &matching.hits[0];
                let shown = &hit.references[0];
                // Names may identify several accession records. Group display
                // names by species while retaining every source in the result.
                let mut names: BTreeMap<&str, BTreeSet<String>> = BTreeMap::new();
                let mut records = 0;
                for reference in matching.hits.iter().flat_map(|h| &h.references) {
                    names
                        .entry(&reference.species)
                        .or_default()
                        .insert(gene_name(reference));
                    records += 1;
                }
                let names = names
                    .into_iter()
                    .map(|(species, genes)| {
                        format!(
                            "{}: {}",
                            display_name(species),
                            genes
                                .into_iter()
                                .map(|gene| display_name(&gene))
                                .collect::<Vec<_>>()
                                .join(", ")
                        )
                    })
                    .collect::<Vec<_>>()
                    .join("; ");
                let record_label = if records == 1 { "record" } else { "records" };
                details.push_str(&format!(
                    "{segment} input similarity (shown: {} {}): score {}; known identity {:.3}; germline/input coverage {:.3}/{:.3}\nTied {segment} references ({records} {record_label}): {names}\n",
                    display_name(&shown.species), display_name(&gene_name(shown)),
                    matching.score, hit.known_identity, hit.reference_coverage, hit.query_coverage,
                ));
            }
        }
        for diagnostic in &self.diagnostics {
            details.push_str(&display_name(diagnostic));
            details.push('\n');
        }
        let mut output = summary(&intro, imputed, &details, width, color);
        output.push_str(&body);
        Ok(output.trim_end_matches('\n').into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use unicode_width::UnicodeWidthStr;
    const SEQUENCE: &str = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

    #[test]
    fn antibody_display_has_two_gene_labels_and_bounded_region_annotations() {
        let ab = number_antibody(
            SEQUENCE,
            &NumberingOptions {
                name: "野生型\nα".into(),
                ..Default::default()
            },
        )
        .unwrap();
        for width in [40, 80] {
            let text = ab.render(width, false, true).unwrap();
            assert!(text.lines().all(|l| l.width() <= width));
            assert!(text.contains("野生型\\nα"));
            assert!(text.contains("IGHV") && text.contains("IGHJ"));
            assert!(!text.contains('\x1b'));
        }
        let colored = ab.render(80, true, false).unwrap();
        assert!(colored.contains("\x1b["));
        assert!(!serde_json::to_string(&ab).unwrap().contains("\\u001b"));
        assert!(ab.render(1, false, true).is_err());
    }

    #[test]
    fn comparison_reference_is_first_without_germline_rows_or_stored_reordering() {
        let first = number_antibody(
            SEQUENCE,
            &NumberingOptions {
                name: "first".into(),
                ..Default::default()
            },
        )
        .unwrap();
        let mut second = first.clone();
        second.name = "second".into();
        let alignment = align_antibodies(vec![first, second], 0).unwrap();
        let text = alignment.render(80, false, false, Some(1)).unwrap();
        assert!(text.find("second").unwrap() < text.find("first").unwrap());
        assert!(!text.contains("IGHV") && !text.contains("IGHJ"));
        assert_eq!(alignment.antibodies[0].name, "first");
        assert_eq!(alignment.reference_index, 0);
        assert!(alignment.render(80, false, true, Some(2)).is_err());
    }
    // Decode the terminal cells so assertions cover visible columns and active
    // backgrounds, rather than depending on a particular ANSI escape grouping.
    fn styled_cells(line: &str) -> Vec<(char, Option<u8>, bool)> {
        let mut result = Vec::new();
        let (mut background, mut inverse) = (None, false);
        let mut text = line;
        while !text.is_empty() {
            if let Some(escape) = text.strip_prefix("\x1b[") {
                let end = escape.find('m').unwrap();
                for code in escape[..end].split(';').map(|c| c.parse::<u8>().unwrap()) {
                    match code {
                        0 => {
                            background = None;
                            inverse = false;
                        }
                        7 => inverse = true,
                        27 => inverse = false,
                        40..=47 | 100..=107 => background = Some(code),
                        49 => background = None,
                        _ => {}
                    }
                }
                text = &escape[end + 1..];
            } else {
                let c = text.chars().next().unwrap();
                result.push((c, background, inverse));
                text = &text[c.len_utf8()..];
            }
        }
        result
    }

    #[test]
    fn summary_legends_and_all_tied_names_survive_styled_wrapping() {
        let ab = number_antibody(
            &SEQUENCE.replace("REGTTGKPIGAFAH", "KDRGGYFDY"),
            &NumberingOptions {
                name: "野生型".into(),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            ab.j_match
                .as_ref()
                .unwrap()
                .hits
                .iter()
                .map(|h| h.references.len())
                .sum::<usize>(),
            3
        );
        let before = serde_json::to_string(&ab).unwrap();
        for width in [40, 80] {
            let colored = ab.render(width, true, false).unwrap();
            let plain = ab.render(width, false, false).unwrap();
            let decoded: String = styled_cells(&colored).iter().map(|c| c.0).collect();
            assert_eq!(decoded, plain);
            assert!(plain.lines().all(|l| l.width() <= width));
            let summary = colored.split("\n\n").next().unwrap();
            // Each wrapped line must restore terminal styles before its newline.
            let cells: Vec<_> = summary.lines().flat_map(styled_cells).collect();
            let text: String = cells.iter().map(|c| c.0).collect();
            // Locate by characters because names can occupy several UTF-8 bytes.
            let byte = text.find("imputed residues: 0").unwrap();
            let start = text[..byte].chars().count();
            assert!(
                cells[start..start + "imputed residues: 0".len()]
                    .iter()
                    .all(|c| c.1 == Some(43) && !c.2)
            );
            assert!(cells[..start].iter().all(|c| c.1 != Some(43) && !c.2));
            assert!(
                cells[start + "imputed residues: 0".len()..]
                    .iter()
                    .all(|c| c.1 != Some(43) && !c.2)
            );
            for (label, background) in [("1 = CDR1", 47), ("2 = CDR2", 105), ("3 = CDR3", 106)] {
                let byte = text.find(label).unwrap();
                let start = text[..byte].chars().count();
                assert!(
                    cells[start..start + label.len()]
                        .iter()
                        .all(|c| c.1 == Some(background))
                );
            }
            for (segment, matching) in [("V", &ab.v_match), ("J", &ab.j_match)] {
                let references: Vec<_> = matching
                    .as_ref()
                    .unwrap()
                    .hits
                    .iter()
                    .flat_map(|h| &h.references)
                    .collect();
                assert!(text.contains(&format!(
                    "Tied {segment} references ({} record",
                    references.len()
                )));
                for r in references {
                    assert!(text.contains(&r.species));
                    assert!(text.contains(&format!("{}*{}", r.gene, r.allele)));
                }
            }
            assert!(text.contains(&format!(
                "Numbered domain in supplied input: {}–{} (1-based)",
                ab.domain_span.0 + 1,
                ab.domain_span.1,
            )));
        }
        assert_eq!(serde_json::to_string(&ab).unwrap(), before);
    }

    #[test]
    fn single_display_places_source_rulers_above_input_and_stitched_germlines() {
        let sequence = format!("AAAAAA{}AAAAAA", SEQUENCE.replace("GGSFSTY", "GGGSGGSFSTY"));
        let ab = number_antibody(
            &sequence,
            &NumberingOptions {
                name: "input".into(),
                ..Default::default()
            },
        )
        .unwrap();
        let before = serde_json::to_string(&ab).unwrap();
        for width in [80, 250] {
            let text = ab.render(width, false, true).unwrap();
            let lines: Vec<_> = text.lines().collect();
            let mut input_position = 0usize;
            let mut reference_sequence = String::new();
            let mut operations = String::new();
            for (line_index, line) in lines
                .iter()
                .enumerate()
                .filter(|(_, l)| l.starts_with("input "))
            {
                let cells = line.split_whitespace().nth(2).unwrap();
                let prefix = line.find(cells).unwrap();
                let germline = lines[line_index + 2];
                assert!(
                    germline.starts_with(
                        &ab.v_match.as_ref().unwrap().hits[0]
                            .alignment
                            .reference_name
                    )
                );
                assert!(!lines[line_index + 1].contains("CDR"));
                let (row, j_name) = germline.rsplit_once("  ").unwrap();
                let (row, last) = row.rsplit_once("  ").unwrap();
                assert_eq!(
                    j_name,
                    ab.j_match.as_ref().unwrap().hits[0]
                        .alignment
                        .reference_name
                );
                assert!(!row.ends_with(' '));
                let visible = &row[prefix..];
                let reference_cells = format!("{visible:width$}", width = cells.len());
                let ops = &lines[line_index + 3][prefix..prefix + cells.len()];
                reference_sequence.push_str(&reference_cells);
                assert_eq!(
                    last.parse::<usize>().unwrap(),
                    reference_sequence
                        .bytes()
                        .filter(u8::is_ascii_alphabetic)
                        .count()
                );
                operations.push_str(ops);
                for (column, residue) in cells.bytes().enumerate() {
                    if residue != b'-' {
                        input_position += 1;
                        if input_position.is_multiple_of(10) {
                            let label = input_position.to_string();
                            let stop = prefix + column + 1;
                            assert_eq!(&lines[line_index - 1][stop - label.len()..stop], label);
                        }
                    }
                    match ops.as_bytes()[column] {
                        b'+' => assert_eq!(residue, b'-'),
                        b'-' => assert_eq!(reference_cells.as_bytes()[column], b'-'),
                        _ => {}
                    }
                }
            }
            assert_eq!(input_position, sequence.len());
            assert!(reference_sequence.starts_with("      "));
            assert!(reference_sequence.ends_with("      "));
            assert!(
                operations.contains('-'),
                "insertions in the top input are germline deletions"
            );
            assert!(operations.starts_with("      ") && operations.ends_with("      "));
            assert!(text.lines().all(|line| line.width() <= width));
            if width == 250 {
                let index = lines.iter().position(|l| l.starts_with("input ")).unwrap();
                let input_cells = lines[index].split_whitespace().nth(2).unwrap();
                let prefix = lines[index].find(input_cells).unwrap();
                let mut column = 0;
                let mut position = 0usize;
                for matching in [&ab.v_match, &ab.j_match] {
                    let hit = &matching.as_ref().unwrap().hits[0];
                    for residue in hit.alignment.reference.bytes() {
                        while reference_sequence.as_bytes()[column] == b' '
                            || reference_sequence.as_bytes()[column] == b'-'
                        {
                            column += 1;
                        }
                        assert_eq!(reference_sequence.as_bytes()[column], residue);
                        position += 1;
                        if position.is_multiple_of(10) {
                            let label = position.to_string();
                            let stop = prefix + column + 1;
                            assert_eq!(&lines[index + 1][stop - label.len()..stop], label);
                        }
                        column += 1;
                    }
                }
            }
        }
        assert_eq!(serde_json::to_string(&ab).unwrap(), before);
    }

    #[test]
    fn shared_bands_follow_reference_but_rulers_and_imputation_follow_each_input() {
        let first = number_antibody(
            SEQUENCE,
            &NumberingOptions {
                name: "original".into(),
                ..Default::default()
            },
        )
        .unwrap();
        let partial = number_antibody(
            &SEQUENCE[5..],
            &NumberingOptions {
                name: "shortened".into(),
                scheme: Some(NumberingScheme::Imgt),
                cdr_definition: CdrDefinition::Chothia,
                ..Default::default()
            },
        )
        .unwrap();
        let v = &partial.v_match.as_ref().unwrap().hits[0].references[0].id;
        let second = partial.impute(Some(v), None).unwrap();
        assert!(second.residues.iter().any(|r| r.input_index.is_none()));
        let alignment = align_antibodies(vec![first, second], 0).unwrap();
        for reference in [0, 1] {
            let text = alignment.render(250, true, true, Some(reference)).unwrap();
            let lines: Vec<_> = text.lines().collect();
            let selected = &alignment.antibodies[reference];
            let imputed = alignment
                .antibodies
                .iter()
                .flat_map(|a| &a.residues)
                .filter(|r| r.input_index.is_none())
                .count();
            let header: String = styled_cells(lines[2]).iter().map(|c| c.0).collect();
            assert_eq!(header, format!("Total imputed residues: {imputed}"));
            assert!(text.contains(&format!(
                "{} numbering; {} CDR definition",
                selected.scheme, selected.cdr_definition
            )));
            let top = lines
                .iter()
                .position(|l| l.starts_with(&format!("{} ", selected.name)))
                .unwrap();
            let decoded = styled_cells(lines[top]);
            let plain: String = decoded.iter().map(|c| c.0).collect();
            let cells = plain.split_whitespace().nth(2).unwrap();
            let prefix = plain.find(cells).unwrap();
            let marker = styled_cells(lines[top - 2]);
            for offset in [top - 1, top, top + 1, top + 2] {
                let row = styled_cells(lines[offset]);
                for column in prefix..prefix + alignment.positions.len() {
                    if row[column].1 != Some(43) {
                        assert_eq!(
                            row[column].1, marker[column].1,
                            "CDR band at column {column}"
                        );
                    }
                }
            }
            assert!(marker.iter().any(|c| c.1 == Some(47)));
            assert!(marker.iter().any(|c| c.1 == Some(105)));
            assert!(marker.iter().any(|c| c.1 == Some(106)));
            assert!(
                styled_cells(lines[top + 3])
                    .iter()
                    .all(|c| c.1.is_none() && !c.2)
            );
            for antibody in &alignment.antibodies {
                let index = lines
                    .iter()
                    .position(|l| l.starts_with(&format!("{} ", antibody.name)))
                    .unwrap();
                let row = styled_cells(lines[index]);
                let ruler: String = styled_cells(lines[index - 1]).iter().map(|c| c.0).collect();
                for residue in &antibody.residues {
                    let column = alignment
                        .positions
                        .iter()
                        .position(|p| *p == residue.position)
                        .unwrap();
                    assert_eq!(
                        row[prefix + column].1 == Some(43),
                        residue.input_index.is_none()
                    );
                    assert!(!row[prefix + column].2);
                    if let Some(position) = residue
                        .input_index
                        .map(|i| i + 1)
                        .filter(|i| i.is_multiple_of(10))
                    {
                        let label = position.to_string();
                        let stop = prefix + column + 1;
                        assert_eq!(&ruler[stop - label.len()..stop], label);
                    }
                }
            }
        }
        assert_eq!(alignment.reference_index, 0);
        assert_eq!(alignment.antibodies[0].cdr_definition, "imgt");
        assert_eq!(alignment.antibodies[1].cdr_definition, "chothia");
    }
}
