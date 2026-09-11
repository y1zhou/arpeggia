use super::*;
use std::collections::{HashMap, HashSet};

/// Antibodies aligned by numbered positions, with rows stored in input order.
#[derive(Clone, Debug, Serialize)]
#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object, module = "arpeggia")
)]
pub struct AntibodyAlignment {
    /// Original numbered antibodies, retaining their CDR definitions and germlines.
    pub antibodies: Vec<NumberedAntibody>,
    /// Ordered union of numbered positions.
    pub positions: Vec<NumberedPosition>,
    /// One plain gapped sequence per input, in original row order.
    pub aligned_sequences: Vec<String>,
    /// Zero-based input row used as the display comparison reference.
    pub reference_index: usize,
}

/// Align one or more antibodies through their shared numbering scheme.
///
/// Inputs must all be heavy or all light chains; kappa/lambda mixtures are
/// supported. Each row retains its CDR convention. This position correspondence
/// has no optimized multiple-sequence-alignment score. The comparison reference
/// changes display direction, without changing stored row/column order.
pub fn align_antibodies(
    antibodies: Vec<NumberedAntibody>,
    reference_index: usize,
) -> ArpeggiaResult<AntibodyAlignment> {
    let reference = antibodies.get(reference_index).ok_or_else(|| {
        ArpeggiaError::InvalidArgument(
            "reference_index must address an input antibody; at least one antibody is required"
                .into(),
        )
    })?;
    let scheme: NumberingScheme = clap::ValueEnum::from_str(&reference.scheme, true)
        .map_err(ArpeggiaError::InvalidArgument)?;
    let chain: Chain = reference
        .chain
        .parse()
        .map_err(|_| ArpeggiaError::InvalidArgument("unsupported antibody chain".into()))?;
    let heavy = reference.chain == "H";
    let mut positions = HashSet::new();
    for antibody in &antibodies {
        if antibody.scheme != reference.scheme
            || !["H", "K", "L"].contains(&antibody.chain.as_str())
            || (antibody.chain == "H") != heavy
        {
            return Err(ArpeggiaError::InvalidArgument("antibody alignment requires one numbering scheme and either all heavy or all light chains".into()));
        }
        if antibody.residues.is_empty()
            || antibody.residues.iter().any(|r| {
                r.position.number == 0
                    || r.position
                        .insertion
                        .is_some_and(|c| !c.is_ascii_uppercase())
                    || !r.amino_acid.is_ascii_uppercase()
            })
        {
            return Err(ArpeggiaError::InvalidArgument(
                "antibodies must contain valid numbered residues".into(),
            ));
        }
        if antibody.residues.windows(2).any(|r| {
            r[0].position.order(scheme.backend(), chain)
                >= r[1].position.order(scheme.backend(), chain)
        }) {
            return Err(ArpeggiaError::InvalidArgument(
                "antibody positions must be unique and in scheme order".into(),
            ));
        }
        positions.extend(antibody.residues.iter().map(|r| r.position));
    }
    let mut positions: Vec<_> = positions.into_iter().collect();
    positions.sort_by_key(|p| p.order(scheme.backend(), chain));
    let indices: HashMap<_, _> = positions.iter().enumerate().map(|(i, p)| (*p, i)).collect();
    let aligned_sequences = antibodies
        .iter()
        .map(|antibody| {
            let mut row = vec![b'-'; positions.len()];
            for residue in &antibody.residues {
                row[indices[&residue.position]] = residue.amino_acid as u8;
            }
            String::from_utf8(row).expect("validated ASCII residues")
        })
        .collect();
    Ok(AntibodyAlignment {
        antibodies,
        positions,
        aligned_sequences,
        reference_index,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    const HEAVY: &str = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

    #[test]
    fn numbered_columns_preserve_input_order_and_reference_independence() {
        let first = number_antibody(
            HEAVY,
            &NumberingOptions {
                name: "first".into(),
                ..Default::default()
            },
        )
        .unwrap();
        let second = number_antibody(
            &HEAVY.replace("REGTT", "REGGGGTT"),
            &NumberingOptions {
                name: "second".into(),
                cdr_definition: CdrDefinition::Chothia,
                scheme: Some(NumberingScheme::Imgt),
                ..Default::default()
            },
        )
        .unwrap();
        let alignment = align_antibodies(vec![first.clone(), second.clone()], 1).unwrap();
        let other = align_antibodies(vec![first, second], 0).unwrap();
        assert_eq!(alignment.antibodies[0].name, "first");
        assert_eq!(alignment.reference_index, 1);
        assert_eq!(alignment.positions, other.positions);
        assert_eq!(alignment.aligned_sequences, other.aligned_sequences);
        assert!(
            alignment
                .positions
                .iter()
                .any(|p| p.number == 112 && p.insertion.is_some())
        );
        assert!(
            alignment.positions.windows(2).all(
                |p| p[0].order(Scheme::IMGT, Chain::IGH) < p[1].order(Scheme::IMGT, Chain::IGH)
            )
        );
        for (row, antibody) in alignment
            .aligned_sequences
            .iter()
            .zip(&alignment.antibodies)
        {
            assert_eq!(row.replace('-', ""), antibody.sequence());
        }
    }

    #[test]
    fn numbered_alignment_accepts_light_mixtures_and_rejects_incompatible_inputs() {
        let kappa = number_antibody("DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIK", &NumberingOptions::default()).unwrap();
        let lambda = number_antibody("QSALTQPASVSGSPGQSITISCTGTTSDVGTYNFVSWYQQHPGKAPKAIIFDVTNRPSGISNRFSGSKFGNTASLTISGLQAEDEADYYCAAYTVASTLLFGGGTKVTVL", &NumberingOptions::default()).unwrap();
        assert!(align_antibodies(vec![kappa.clone(), lambda], 0).is_ok());
        assert!(align_antibodies(vec![kappa.clone()], 0).is_ok());
        assert!(align_antibodies(vec![kappa.clone()], 1).is_err());
        assert!(align_antibodies(Vec::new(), 0).is_err());
        let heavy = number_antibody(HEAVY, &NumberingOptions::default()).unwrap();
        assert!(align_antibodies(vec![kappa, heavy.clone()], 0).is_err());
        let mut other = heavy.clone();
        other.scheme = "martin".into();
        assert!(align_antibodies(vec![heavy, other], 0).is_err());
    }
}
