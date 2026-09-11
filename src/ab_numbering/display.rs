use super::*;
use crate::AlignmentColor;
use crate::seq_alignment::display::{
    DisplayRow, PositionLabel, color_enabled, render_rows, wrap_summary,
};
use std::collections::HashMap;
use std::io::IsTerminal;

fn dimensions(width: Option<usize>, color: AlignmentColor) -> (usize, bool) {
    (
        width.unwrap_or_else(|| {
            terminal_size::terminal_size().map_or(80, |(w, _)| usize::from(w.0))
        }),
        color_enabled(color, std::io::stdout().is_terminal()),
    )
}

fn label(position: Option<NumberedPosition>) -> Option<PositionLabel> {
    position.map(|p| PositionLabel {
        text: p.to_string(),
        tick: p.insertion.is_none() && p.number.is_multiple_of(10),
    })
}

fn region(name: &str) -> u8 {
    match name {
        "CDR1" => 1,
        "CDR2" => 2,
        "CDR3" => 3,
        _ => 0,
    }
}

impl AntibodyAlignment {
    /// Format antibody comparisons with the selected reference first and germlines hidden.
    ///
    /// Width includes names and numbered-position labels. A rendering-only
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
                    .map(|p| label(by_position.contains_key(p).then_some(*p)))
                    .collect();
                let regions = self
                    .positions
                    .iter()
                    .map(|p| by_position.get(p).map_or(0, |r| region(&r.region)))
                    .collect();
                DisplayRow {
                    name: antibody.name.clone(),
                    right_name: String::new(),
                    cells,
                    positions,
                    operations,
                    regions,
                    show_operations: index != reference,
                }
            })
            .collect();
        let body = render_rows(&rows, width, color, rulers)?;
        let mut output = wrap_summary(
            &format!(
                "{} antibodies; {} positions; {} numbering\nCDR regions: 1 = CDR1, 2 = CDR2, 3 = CDR3",
                self.antibodies.len(),
                self.positions.len(),
                self.antibodies[reference].scheme
            ),
            width,
        );
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
    region: u8,
}

impl NumberedAntibody {
    /// Format the input against a combined V/J reference, with distinct CDR backgrounds.
    ///
    /// Gray junction gaps mark unavailable reference sequence. V/J scores remain
    /// separate similarities measured on the supplied input, including after imputation.
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
                        region: 0,
                    },
                );
            }
        }
        let name = |m: &Option<GermlineMatch>, absent: &str| {
            m.as_ref().map_or_else(
                || absent.into(),
                |m| {
                    let r = &m.hits[0].references[0];
                    format!("{}*{}", r.gene, r.allele)
                },
            )
        };
        let rows: Vec<_> = [true, false]
            .into_iter()
            .map(|reference| DisplayRow {
                name: if reference {
                    name(&self.v_match, "V unavailable")
                } else {
                    self.name.clone()
                },
                right_name: if reference {
                    name(&self.j_match, "J unavailable")
                } else {
                    String::new()
                },
                cells: combined
                    .iter()
                    .map(|c| if reference { c.reference } else { c.query })
                    .collect(),
                positions: combined
                    .iter()
                    .map(|c| {
                        label(if reference {
                            c.reference_position
                        } else {
                            c.query_position
                        })
                    })
                    .collect(),
                operations: combined.iter().map(|c| c.operation).collect(),
                regions: combined.iter().map(|c| c.region).collect(),
                show_operations: !reference,
            })
            .collect();
        let body = render_rows(&rows, width, color, rulers)?;
        let mut summary = format!(
            "{} chain; {} numbering; {} CDR definition\nConfidence: {:.3}; matched profile positions: {}\nDomain input span: [{}, {}); imputed residues: {}\nCDR regions: 1 = CDR1, 2 = CDR2, 3 = CDR3",
            self.chain,
            self.scheme,
            self.cdr_definition,
            self.confidence,
            self.matched_profile_positions,
            self.domain_span.0,
            self.domain_span.1,
            self.residues
                .iter()
                .filter(|r| r.input_index.is_none())
                .count()
        );
        for (segment, matching) in [("V", &self.v_match), ("J", &self.j_match)] {
            if let Some(matching) = matching {
                let hit = &matching.hits[0];
                summary.push_str(&format!("\n{segment} input similarity: score {}; known identity {:.3}; reference/query coverage {:.3}/{:.3}; {} tied references",matching.score,hit.known_identity,hit.reference_coverage,hit.query_coverage,matching.hits.iter().map(|h| h.references.len()).sum::<usize>()));
            }
        }
        for diagnostic in &self.diagnostics {
            summary.push('\n');
            summary.push_str(diagnostic);
        }
        let mut output = wrap_summary(&summary, width);
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
}
