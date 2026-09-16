use super::*;
use std::collections::{HashMap, HashSet};

impl NumberedAntibody {
    /// Fill supported missing beginnings of FR1 and ends of FR4 in a new object.
    ///
    /// All tied references must agree on presence and a standard amino acid.
    /// Optional selectors are exact `GermlineReference.id` values among the
    /// corresponding tied hits. Unknown input residues and internal gaps remain
    /// unchanged. Inferred residues have no original input index and retain their
    /// source IDs; the original input, domain span and object are preserved.
    /// Germline matching must have run when constructing the antibody.
    pub fn impute(
        &self,
        v_reference: Option<&str>,
        j_reference: Option<&str>,
    ) -> ArpeggiaResult<Self> {
        if !self.germlines_searched {
            return Err(ArpeggiaError::InvalidArgument(
                "imputation requires germline matching; number the input again with match_germlines enabled".into(),
            ));
        }
        let scheme = clap::ValueEnum::from_str(&self.scheme, true)
            .map_err(ArpeggiaError::InvalidArgument)?;
        let scheme: NumberingScheme = scheme;
        let definition: CdrDefinition = clap::ValueEnum::from_str(&self.cdr_definition, true)
            .map_err(ArpeggiaError::InvalidArgument)?;
        let chain: Chain = self
            .chain
            .parse()
            .map_err(|_| ArpeggiaError::InvalidArgument("unsupported antibody chain".into()))?;
        let first = self
            .residues
            .first()
            .ok_or_else(|| {
                ArpeggiaError::InvalidArgument("cannot impute an empty numbered antibody".into())
            })?
            .position
            .order(scheme.backend(), chain);
        let last = self
            .residues
            .last()
            .unwrap()
            .position
            .order(scheme.backend(), chain);
        let mut result = self.clone();
        result.imputation_attempted = true;
        for (matching, selector, region) in [
            (&self.v_match, v_reference, "FR1"),
            (&self.j_match, j_reference, "FR4"),
        ] {
            let hits: Vec<_> = matching
                .iter()
                .flat_map(|m| &m.hits)
                .filter(|hit| selector.is_none_or(|id| hit.references.iter().any(|r| r.id == id)))
                .collect();
            if hits.is_empty() {
                if selector.is_some() {
                    return Err(ArpeggiaError::InvalidArgument(format!(
                        "{region} reference selector is not among the tied germline references"
                    )));
                }
                result.diagnostics.push(format!(
                    "IMPUTATION_UNAVAILABLE: {region} has no supported germline match"
                ));
                continue;
            }
            let mut maps = Vec::new();
            let mut ids = Vec::new();
            let mut unavailable = false;
            for hit in hits {
                match project(
                    hit,
                    scheme.backend(),
                    definition.backend(scheme),
                    chain,
                    region,
                ) {
                    Ok(map) => maps.push(map),
                    Err(error) => {
                        result.diagnostics.push(format!("IMPUTATION_UNAVAILABLE: {region} reference cannot be converted ({error})"));
                        unavailable = true;
                    }
                }
                ids.extend(
                    hit.references
                        .iter()
                        .filter(|r| selector.is_none_or(|id| r.id == id))
                        .map(|r| r.id.clone()),
                );
            }
            if unavailable {
                continue;
            }
            let endpoint = NumberedPosition {
                number: if region == "FR1" {
                    1
                } else {
                    core::rules(scheme.backend(), chain).last().unwrap().num_end
                },
                insertion: None,
            };
            let outside = if region == "FR1" {
                endpoint.order(scheme.backend(), chain) < first
            } else {
                endpoint.order(scheme.backend(), chain) > last
            };
            if outside && maps.iter().all(|map| !map.contains_key(&endpoint)) {
                result.diagnostics.push(format!(
                    "IMPUTATION_REFERENCE_COVERAGE: no tied reference covers {region} endpoint {endpoint}; no extrapolation beyond available sequence"
                ));
            }
            ids.sort();
            ids.dedup();
            let candidates: HashSet<_> = maps.iter().flat_map(|map| map.keys().copied()).collect();
            let mut unresolved = 0;
            for position in candidates {
                let order = position.order(scheme.backend(), chain);
                if (region == "FR1" && order >= first) || (region == "FR4" && order <= last) {
                    continue;
                }
                let amino_acid = maps[0].get(&position).copied();
                if let Some(amino_acid) = amino_acid.filter(|aa| {
                    germline::known(*aa as u8)
                        && maps.iter().all(|map| map.get(&position) == Some(aa))
                }) {
                    result.residues.push(NumberedResidue {
                        position,
                        amino_acid,
                        input_index: None,
                        region: region.into(),
                        imputed_from: ids.clone(),
                    });
                } else {
                    unresolved += 1;
                }
            }
            if unresolved > 0 {
                result.diagnostics.push(format!("IMPUTATION_UNRESOLVED: {unresolved} {region} positions lack unanimous known reference residues"));
            }
        }
        result
            .residues
            .sort_by_key(|r| r.position.order(scheme.backend(), chain));
        Ok(result)
    }
}

fn project(
    hit: &GermlineHit,
    scheme: Scheme,
    definition: Scheme,
    chain: Chain,
    region: &str,
) -> ArpeggiaResult<HashMap<NumberedPosition, char>> {
    let positions = germline::reference_positions(hit, scheme, chain)?;
    let regions = if definition == scheme {
        positions.clone()
    } else {
        germline::reference_positions(hit, definition, chain)?
    };
    Ok(positions
        .into_iter()
        .enumerate()
        .filter_map(|(i, position)| {
            let position = position?;
            let assigned = regions[i].map_or_else(
                || "FR4".into(),
                |p| core::region(p.number, definition, chain),
            );
            (assigned == region).then(|| (position, hit.alignment.reference.as_bytes()[i] as char))
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    const SEQUENCE: &str = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQRLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARGGGYFDYWGQGTLVTVSS";

    #[test]
    fn terminal_imputation_preserves_supplied_sequence_and_provenance() {
        for scheme in [
            NumberingScheme::Imgt,
            NumberingScheme::Chothia,
            NumberingScheme::Martin,
            NumberingScheme::Aho,
            NumberingScheme::Kabat,
        ] {
            let partial = number_antibody(
                &SEQUENCE[5..SEQUENCE.len() - 3],
                &NumberingOptions {
                    scheme: Some(scheme),
                    species: vec![GermlineSpecies::Human],
                    ..Default::default()
                },
            )
            .unwrap();
            let v = &partial.v_match.as_ref().unwrap().hits[0].references[0].id;
            let j = &partial.j_match.as_ref().unwrap().hits[0].references[0].id;
            let imputed = partial.impute(Some(v), Some(j)).unwrap();
            assert!(
                imputed.sequence().starts_with("QVQLV"),
                "{scheme:?}: {:?}",
                imputed.diagnostics
            );
            assert!(imputed.sequence().ends_with("VSS"));
            assert_eq!(imputed.input_sequence, partial.input_sequence);
            assert_eq!(imputed.domain_span, partial.domain_span);
            assert_eq!(
                imputed
                    .residues
                    .iter()
                    .filter(|r| r.input_index.is_some())
                    .map(|r| r.amino_acid)
                    .collect::<String>(),
                partial.sequence()
            );
            assert!(
                imputed
                    .residues
                    .iter()
                    .filter(|r| r.input_index.is_none())
                    .all(|r| !r.imputed_from.is_empty()
                        && ["FR1", "FR4"].contains(&r.region.as_str()))
            );
            assert!(partial.residues.iter().all(|r| r.input_index.is_some()));
        }
    }

    #[test]
    fn imputation_does_not_replace_unknown_input_or_guess_between_ties() {
        let mut partial = number_antibody(
            &SEQUENCE[5..].replacen('Y', "X", 1),
            &NumberingOptions {
                species: vec![GermlineSpecies::Human],
                ..Default::default()
            },
        )
        .unwrap();
        let matching = partial.v_match.as_mut().unwrap();
        let mut alternative = matching.hits[0].clone();
        alternative.references.truncate(1);
        alternative.references[0].id = "different-terminal-reference".into();
        alternative.alignment.reference.replace_range(..1, "X");
        matching.hits.push(alternative);
        let result = partial.impute(None, None).unwrap();
        assert!(result.sequence().contains('X'));
        assert!(result.residues.iter().all(|r| r.position.number != 1));
        assert!(
            result
                .diagnostics
                .iter()
                .any(|d| d.starts_with("IMPUTATION_UNRESOLVED"))
        );
        assert!(partial.impute(Some("unknown-reference"), None).is_err());
    }

    #[test]
    fn uncovered_reference_ends_are_reported_without_extrapolation() {
        let mut partial = number_antibody(
            &SEQUENCE[5..SEQUENCE.len() - 3],
            &NumberingOptions::default(),
        )
        .unwrap();
        for matching in [&mut partial.v_match, &mut partial.j_match]
            .into_iter()
            .flatten()
        {
            for hit in &mut matching.hits {
                // Model references whose terminal coverage is unavailable, as
                // in the bundled partial V and short alpaca J records.
                for position in &mut hit.imgt_positions {
                    if position.is_some_and(|p| p.number == 1 || p.number == 128) {
                        *position = None;
                    }
                }
            }
        }
        let result = partial.impute(None, None).unwrap();
        assert!(
            result
                .residues
                .iter()
                .all(|r| ![1, 128].contains(&r.position.number))
        );
        assert_eq!(
            result
                .diagnostics
                .iter()
                .filter(|d| d.starts_with("IMPUTATION_REFERENCE_COVERAGE"))
                .count(),
            2
        );
        assert_eq!(result.input_sequence, partial.input_sequence);
    }
}
