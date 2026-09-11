use super::*;
use immunum::alignment::{AlignBuffer, AlignedPosition, Alignment};
use immunum::{Insertion, NumberingRule, ScoringMatrix};
use std::sync::OnceLock;

fn failure(message: impl Into<String>) -> ArpeggiaError {
    ArpeggiaError::Calculation(format!("antibody numbering: {}", message.into()))
}

fn best_alignment(sequence: &str, buffer: &mut AlignBuffer) -> ArpeggiaResult<(Chain, Alignment)> {
    let scoring_sequence = if sequence.contains(['U', 'O']) {
        std::borrow::Cow::Owned(sequence.replace('U', "C").replace('O', "K"))
    } else {
        std::borrow::Cow::Borrowed(sequence)
    };
    static PROFILES: OnceLock<Result<Vec<(Chain, ScoringMatrix)>, String>> = OnceLock::new();
    let profiles = PROFILES
        .get_or_init(|| {
            [Chain::IGH, Chain::IGK, Chain::IGL]
                .into_iter()
                .map(|chain| {
                    ScoringMatrix::load(chain)
                        .map(|matrix| (chain, matrix))
                        .map_err(|e| e.to_string())
                })
                .collect()
        })
        .as_ref()
        .map_err(failure)?;
    let mut best: Option<(Chain, Alignment)> = None;
    for (chain, matrix) in profiles {
        let candidate = immunum::align(&scoring_sequence, &matrix.positions, Some(buffer));
        if best
            .as_ref()
            .is_none_or(|(_, previous)| candidate.score > previous.score)
        {
            best = Some((*chain, candidate));
        }
    }
    Ok(best.expect("three bundled antibody profiles"))
}

fn evidence(alignment: &Alignment) -> (f32, usize) {
    let confidence = if alignment.max_confidence_score > 0.0 {
        (alignment.confidence_score / alignment.max_confidence_score).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let mut seen = [false; 256];
    for position in &alignment.positions[alignment.query_start..=alignment.query_end] {
        if let AlignedPosition::Aligned(n) = position {
            seen[usize::from(*n)] = true;
        }
    }
    (confidence, seen.into_iter().filter(|p| *p).count())
}

pub(super) fn rules(scheme: Scheme, chain: Chain) -> &'static [NumberingRule] {
    use immunum::{aho::*, chothia::*, imgt::*, kabat::*, martin::*};
    match (scheme, chain) {
        (Scheme::IMGT, _) => IMGT_RULES,
        (Scheme::Martin, Chain::IGH) => MARTIN_HEAVY_RULES,
        (Scheme::Martin, _) => MARTIN_LIGHT_RULES,
        (Scheme::Kabat, Chain::IGH) => KABAT_HEAVY_RULES,
        (Scheme::Kabat, _) => KABAT_LIGHT_RULES,
        (Scheme::Chothia, Chain::IGH) => CHOTHIA_HEAVY_RULES,
        (Scheme::Chothia, _) => CHOTHIA_LIGHT_RULES,
        (Scheme::Aho, Chain::IGH) => AHO_HEAVY_RULES,
        (Scheme::Aho, Chain::IGK) => AHO_KAPPA_RULES,
        (Scheme::Aho, _) => AHO_LAMBDA_RULES,
    }
}

pub(super) fn convert(
    alignment: &Alignment,
    scheme: Scheme,
    chain: Chain,
    input_length: usize,
) -> ArpeggiaResult<Vec<NumberedPosition>> {
    // Count-based rules cannot distinguish a cut loop from an internal deletion.
    // Require their whole source window for both numbering and CDR conversion.
    // AHo's rule starting at 1 explicitly handles missing N-terminal FR1 residues.
    // Germline V/J segments use convert_states directly: they are not full domains.
    for rule in rules(scheme, chain) {
        if rule.align_start > 1
            && !matches!(rule.insertion, Insertion::None)
            && (alignment.cons_start > rule.align_start || alignment.cons_end < rule.align_end)
        {
            return Err(failure(format!(
                "partial domain lacks {scheme} conversion context at IMGT positions {}–{}; supply the surrounding framework sequence",
                rule.align_start, rule.align_end
            )));
        }
    }
    let aligned = &alignment.positions[alignment.query_start..=alignment.query_end];
    let mut converted = convert_states(aligned, scheme, chain)?;
    // Preserve Immunum's AHo light-chain tail rule; this residue has no raw
    // profile state when the alignment ends at IMGT 127.
    if scheme == Scheme::Aho
        && chain != Chain::IGH
        && converted.last()
            == Some(&NumberedPosition {
                number: 148,
                insertion: None,
            })
        && alignment.query_start + converted.len() < input_length
    {
        converted.push(NumberedPosition {
            number: 149,
            insertion: None,
        });
    }
    Ok(converted)
}

pub(super) fn convert_states(
    aligned: &[AlignedPosition],
    scheme: Scheme,
    chain: Chain,
) -> ArpeggiaResult<Vec<NumberedPosition>> {
    // Match upstream's insertion inheritance to validate rule capacities before
    // calling its unchecked converter. States must be ordered H/K/L profile or
    // germline positions, with insertions following an aligned position.
    let mut previous = 0;
    let inherited: Vec<_> = aligned
        .iter()
        .map(|state| {
            if let AlignedPosition::Aligned(n) = state {
                previous = *n;
            }
            previous
        })
        .collect();
    let rules = rules(scheme, chain);
    let end = inherited.partition_point(|n| *n <= rules.last().unwrap().align_end);
    let mut index = 0;
    for rule in rules {
        let start = index;
        while index < end && rule.contains(inherited[index]) {
            index += 1;
        }
        let run = &inherited[start..index];
        if run.is_empty() {
            continue;
        }
        let base = usize::from(rule.num_end - rule.num_start + 1);
        let extra = run.len().saturating_sub(base);
        let supported = match rule.insertion {
            Insertion::Sequential(_) => extra <= 26,
            Insertion::Symmetric { .. } => extra <= 52,
            Insertion::None => run.chunk_by(|a, b| a == b).all(|r| r.len() <= 27),
        };
        if !supported {
            return Err(failure(format!(
                "unsupported insertion length in {scheme} source positions {}–{} (single-letter insertion codes only)",
                rule.align_start, rule.align_end
            )));
        }
        if !matches!(rule.insertion, Insertion::None)
            && base.saturating_sub(run.len()) > rule.deletion_order.len()
        {
            return Err(failure(format!("unsupported deletion in {scheme}")));
        }
    }
    if index != end {
        return Err(failure(format!(
            "{scheme} cannot represent an internal profile position"
        )));
    }
    let converted: Vec<NumberedPosition> =
        immunum::numbering::apply_numbering(&aligned[..end], scheme, chain)
            .into_iter()
            .map(Into::into)
            .collect();
    if converted.len() != end
        || converted
            .iter()
            .any(|p| p.number == 0 || p.insertion.is_some_and(|c| !c.is_ascii_uppercase()))
        || converted
            .windows(2)
            .any(|p| p[0].order(scheme, chain) >= p[1].order(scheme, chain))
    {
        return Err(failure(format!("inconsistent {scheme} conversion")));
    }
    Ok(converted)
}

pub(super) fn region(position: u8, scheme: Scheme, chain: Chain) -> String {
    // Published AHo structural loops; upstream explicitly marks its table
    // unverified. Other definitions reuse upstream's cited region tables.
    let definition = if scheme == Scheme::Aho {
        immunum::types::RegionDefinition {
            fr1_end: 24,
            cdr1_end: 40,
            fr2_end: 57,
            cdr2_end: 77,
            fr3_end: 108,
            cdr3_end: 137,
            fr4_end: 149,
        }
    } else {
        immunum::numbering::regions_for(scheme, chain)
    };
    definition
        .region(position)
        .expect("validated numbered position within region table")
        .to_string()
}

pub(super) fn number(
    sequence: &str,
    options: &NumberingOptions,
) -> ArpeggiaResult<NumberedAntibody> {
    crate::seq_alignment::encode(sequence)?;
    if !(30..=10000).contains(&sequence.len()) {
        return Err(ArpeggiaError::InvalidArgument(
            "antibody inputs must contain 30–10000 amino acids".into(),
        ));
    }
    let sequence = sequence.to_ascii_uppercase();
    let mut buffer = AlignBuffer::new();
    let (chain, alignment) = best_alignment(&sequence, &mut buffer)?;
    let (confidence, matched) = evidence(&alignment);
    if confidence < 0.5 || matched < 30 {
        return Err(failure(format!(
            "no supported antibody domain (confidence {confidence:.3}, {matched} matched profile positions; requires 0.5 and 30)"
        )));
    }
    for flank in [
        &sequence[..alignment.query_start],
        &sequence[alignment.query_end + 1..],
    ] {
        if flank.len() >= 30 {
            let (_, other) = best_alignment(flank, &mut buffer)?;
            let (confidence, matched) = evidence(&other);
            if confidence >= 0.5 && matched >= 30 {
                return Err(failure(
                    "multiple variable domains detected; supply one domain per input",
                ));
            }
        }
    }
    if alignment.cons_start > 26 || alignment.cons_end < 118 {
        return Err(failure(
            "partial domain lacks the complete FR1-through-FR4 core; only terminal FR1/FR4 truncations are supported",
        ));
    }
    let scheme = options.scheme.unwrap_or_default();
    let definition = options.cdr_definition.backend(scheme);
    let positions = convert(&alignment, scheme.backend(), chain, sequence.len())?;
    let definition_positions = if definition == scheme.backend() {
        positions.clone()
    } else {
        convert(&alignment, definition, chain, sequence.len())?
    };
    let start = alignment.query_start;
    let end = start + positions.len();
    let residues = positions
        .into_iter()
        .enumerate()
        .map(|(offset, position)| {
            // Schemes can differ only in the terminal light-chain residue. Any
            // residue beyond the definition's supported suffix is still FR4.
            let region = definition_positions
                .get(offset)
                .map_or_else(|| "FR4".into(), |p| region(p.number, definition, chain));
            NumberedResidue {
                position,
                amino_acid: sequence.as_bytes()[start + offset] as char,
                input_index: Some(start + offset),
                region,
                imputed_from: Vec::new(),
            }
        })
        .collect();
    let mut diagnostics = Vec::new();
    if alignment.cons_start > 1 || alignment.cons_end < 127 {
        diagnostics
            .push("PARTIAL_DOMAIN: only observed terminal framework residues were numbered".into());
    }
    if sequence.contains(['U', 'O']) {
        diagnostics
            .push("SEQUENCE_SCORING_ALIAS: U/O score as C/K; original symbols are retained".into());
    }
    let mut antibody = NumberedAntibody {
        name: if options.name.is_empty() {
            "Seq001".into()
        } else {
            options.name.clone()
        },
        input_sequence: sequence,
        domain_span: (start, end),
        chain: chain.to_string(),
        scheme: scheme.name().into(),
        cdr_definition: definition.to_string().to_ascii_lowercase(),
        residues,
        confidence,
        matched_profile_positions: matched,
        diagnostics,
        v_match: None,
        j_match: None,
    };
    super::germline::matches(&mut antibody, &alignment, chain, &options.species)?;
    Ok(antibody)
}

#[cfg(test)]
mod tests {
    use super::*;
    const HEAVY: &str = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";
    const KAPPA: &str = "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIK";
    // AntPack's pinned numbering fixture 788; see the engine comparison report.
    const LAMBDA: &str = "QSALTQPASVSGSPGQSITISCTGTTSDVGTYNFVSWYQQHPGKAPKAIIFDVTNRPSGISNRFSGSKFGNTASLTISGLQAEDEADYYCAAYTVASTLLFGGGTKVTVL";

    #[test]
    fn domain_spans_and_regions_follow_the_requested_convention() {
        for scheme in [
            NumberingScheme::Imgt,
            NumberingScheme::Martin,
            NumberingScheme::Aho,
            NumberingScheme::Kabat,
        ] {
            for sequence in [HEAVY, KAPPA, LAMBDA] {
                let options = NumberingOptions {
                    scheme: Some(scheme),
                    ..Default::default()
                };
                let result = number_antibody(&format!("AAAAAA{sequence}AAAAAA"), &options).unwrap();
                assert_eq!(result.domain_span.0, 6);
                assert_eq!(result.residues.len(), result.domain_span.1 - 6);
                assert!(result.matched_profile_positions >= 80);
                assert_eq!(
                    result.region_sequence("FR1").chars().next(),
                    sequence.chars().next()
                );
                assert!(result.residues.iter().all(|r| r.amino_acid as u8
                    == result.input_sequence.as_bytes()[r.input_index.unwrap()]));
            }
        }
        let martin = number_antibody(
            HEAVY,
            &NumberingOptions {
                scheme: Some(NumberingScheme::Martin),
                ..Default::default()
            },
        )
        .unwrap();
        let chothia = number_antibody(
            HEAVY,
            &NumberingOptions {
                scheme: Some(NumberingScheme::Martin),
                cdr_definition: CdrDefinition::Chothia,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            martin
                .residues
                .iter()
                .map(|r| r.position)
                .collect::<Vec<_>>(),
            chothia
                .residues
                .iter()
                .map(|r| r.position)
                .collect::<Vec<_>>()
        );
        assert!(martin.region_sequence("CDR1").len() > chothia.region_sequence("CDR1").len());
        assert!(
            number_antibody(
                HEAVY,
                &NumberingOptions {
                    cdr_definition: CdrDefinition::Chothia,
                    ..Default::default()
                }
            )
            .is_err()
        );
    }

    #[test]
    fn recognition_rejects_tandem_domains_and_severe_fragments() {
        for input in [
            "A".repeat(100),
            HEAVY[..60].into(),
            format!("{HEAVY}GGGGSGGGGS{KAPPA}"),
        ] {
            assert!(number_antibody(&input, &NumberingOptions::default()).is_err());
        }
        let full = number_antibody(HEAVY, &NumberingOptions::default()).unwrap();
        let partial =
            number_antibody(&HEAVY[5..HEAVY.len() - 5], &NumberingOptions::default()).unwrap();
        assert_eq!(
            partial
                .residues
                .iter()
                .map(|r| r.position)
                .collect::<Vec<_>>(),
            full.residues[5..full.residues.len() - 5]
                .iter()
                .map(|r| r.position)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn terminal_truncations_preserve_surviving_numbering_and_regions() {
        for sequence in [HEAVY, KAPPA, LAMBDA] {
            for scheme in [
                NumberingScheme::Imgt,
                NumberingScheme::Martin,
                NumberingScheme::Aho,
                NumberingScheme::Kabat,
            ] {
                let options = NumberingOptions {
                    scheme: Some(scheme),
                    ..Default::default()
                };
                let full = number_antibody(sequence, &options).unwrap();
                for removed in [1, 5, 9, 20, 22, 23, 24, 25, 26] {
                    let partial = number_antibody(&sequence[removed..sequence.len() - 3], &options);
                    if removed <= 5 {
                        assert!(partial.is_ok(), "{scheme:?}: {partial:?}");
                    }
                    if let Ok(partial) = partial {
                        for residue in &partial.residues {
                            let original_index = residue.input_index.unwrap() + removed;
                            let original = full
                                .residues
                                .iter()
                                .find(|r| r.input_index == Some(original_index))
                                .unwrap();
                            assert_eq!(
                                (residue.position, &residue.region),
                                (original.position, &original.region),
                                "{scheme:?} {} truncated by {removed}, input index {original_index}",
                                full.chain
                            );
                        }
                    }
                }
            }
        }
        let kabat = NumberingOptions {
            scheme: Some(NumberingScheme::Kabat),
            ..Default::default()
        };
        assert_eq!(number_antibody(HEAVY, &kabat).unwrap().cdr1(), "TYALS");
        assert!(
            number_antibody(&HEAVY[23..], &kabat)
                .unwrap_err()
                .to_string()
                .contains("IMGT positions 24–40")
        );
        // A safe numbering scheme must not bypass an unsafe CDR-definition conversion.
        for scheme in [
            NumberingScheme::Imgt,
            NumberingScheme::Martin,
            NumberingScheme::Aho,
            NumberingScheme::Kabat,
        ] {
            let options = NumberingOptions {
                scheme: Some(scheme),
                cdr_definition: CdrDefinition::Kabat,
                ..Default::default()
            };
            assert!(number_antibody(&HEAVY[23..], &options).is_err());
        }
    }

    #[test]
    fn conversion_checks_insertion_limits_before_upstream_arithmetic() {
        let (_, mut alignment) = best_alignment(HEAVY, &mut AlignBuffer::new()).unwrap();
        for extra in [52, 53, 400] {
            alignment.positions = (1..=128).map(AlignedPosition::Aligned).collect();
            alignment
                .positions
                .splice(111..111, vec![AlignedPosition::Insertion(); extra]);
            alignment.query_start = 0;
            alignment.query_end = alignment.positions.len() - 1;
            let converted = convert(
                &alignment,
                Scheme::IMGT,
                Chain::IGH,
                alignment.positions.len(),
            );
            assert_eq!(converted.is_ok(), extra == 52);
        }
    }

    #[test]
    fn light_chain_conversion_retains_the_correct_terminal_span() {
        // AntPack COVID fixture 17116: upstream reports 107 residues but emits
        // only 106 Martin/Kabat positions, leaving the terminal R unsupported.
        let sequence = "DIQMTQSPSSLSASVGDRVTITCQASQDISNYLNWYQQKPGKAPKLLIYDASNLETGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQYDNLPRFGPGTKVDIKR";
        let (_, mut alignment) = best_alignment(sequence, &mut AlignBuffer::new()).unwrap();
        assert_eq!(alignment.cons_end, 128);
        for scheme in [Scheme::Martin, Scheme::Kabat] {
            assert_eq!(
                convert(&alignment, scheme, Chain::IGK, sequence.len())
                    .unwrap()
                    .len(),
                sequence.len() - 1
            );
        }
        let aho = convert(&alignment, Scheme::Aho, Chain::IGK, sequence.len()).unwrap();
        assert_eq!(aho.len(), sequence.len());
        assert_eq!(aho.last().unwrap().number, 149);
        alignment.positions.pop();
        alignment.query_end -= 1;
        assert_eq!(
            convert(&alignment, Scheme::Aho, Chain::IGK, sequence.len()).unwrap(),
            aho
        );
    }

    #[test]
    fn scoring_aliases_preserve_input_symbols() {
        let aliased = HEAVY.replace('C', "U").replace('K', "O");
        let full = number_antibody(HEAVY, &NumberingOptions::default()).unwrap();
        let result =
            number_antibody(&aliased.to_ascii_lowercase(), &NumberingOptions::default()).unwrap();
        assert_eq!(result.input_sequence, aliased);
        assert_eq!(result.confidence, full.confidence);
        assert_eq!(result.domain_span, full.domain_span);
        assert!(result.sequence().contains('U'));
    }
}
