//! Rigid protein-structure superposition and atom correspondence.

mod correspondence;
mod kabsch;
pub(crate) mod pairwise;
mod selection;

use crate::structure::{select_conformers, selected_model};
use crate::{Analysis, ArpeggiaError, ArpeggiaResult};
pub use correspondence::{ChainAlignment, ResiduePair};
pub use kabsch::kabsch_rmsd;
use kabsch::{
    fit_prepared_transform, kabsch_prepared_rmsd, kabsch_prepared_selected_rmsd,
    prepare_coordinate_union,
};
use nalgebra as na;
use pdbtbx::PDB;
pub use selection::{AtomSubset, validate_rmsd_selections};
use selection::{ResidueSelector, select_coordinate_union, validate_selection_correspondence};
use serde::Serialize;
use std::collections::BTreeMap;

/// Options for correspondence, rigid fitting, and independent RMSD evaluation.
#[derive(Clone, Debug, Serialize)]
pub struct RmsdOptions {
    /// Model serial; zero independently selects the first model in each input.
    pub model_num: usize,
    /// Reference residue selection for fitting (both inputs when alignment is off).
    pub superpose_residues: String,
    /// Reference residue selection for evaluation (both inputs when alignment is off).
    pub rmsd_residues: String,
    /// Shared atom subset.
    pub atoms: AtomSubset,
    /// Establish residue correspondence through observed sequences.
    pub align_seqs: bool,
    /// Complete explicit reference-to-query chain map; empty enables inference.
    pub chain_map: BTreeMap<String, String>,
    /// Final residue alignment settings; inference always uses semi-global scores.
    pub alignment: crate::SeqAlignOptions,
    /// Maximum rejection/refit cycles after the initial fit.
    pub refine_cycles: usize,
    /// Reject distances greater than this multiple of current fitting RMSD.
    pub refine_cutoff: f64,
}
impl Default for RmsdOptions {
    fn default() -> Self {
        Self {
            model_num: 0,
            superpose_residues: String::new(),
            rmsd_residues: String::new(),
            atoms: AtomSubset::Ca,
            align_seqs: false,
            chain_map: BTreeMap::new(),
            alignment: Default::default(),
            refine_cycles: 0,
            refine_cutoff: 2.0,
        }
    }
}

/// Compact result separating retained fitting pairs from the full evaluation set.
#[derive(Clone, Debug, Serialize)]
#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object, module = "arpeggia")
)]
pub struct RmsdResult {
    /// Final RMSD over all mapped evaluation pairs, in Angstroms.
    pub rmsd: f64,
    /// Final RMSD over retained fitting pairs, in Angstroms.
    pub core_rmsd: f64,
    /// Fitting RMSD before any rejection, in Angstroms.
    pub initial_rmsd: f64,
    /// Initial number of fitting atom pairs.
    pub initial_fit_atoms: usize,
    /// Retained number of fitting atom pairs.
    pub retained_fit_atoms: usize,
    /// Number of evaluation atom pairs, unaffected by rejection.
    pub evaluation_atoms: usize,
    /// Rejection inspections performed, including a final unchanged inspection.
    pub cycles: usize,
    /// Requested maximum number of rejection cycles.
    pub refine_cycles: usize,
    /// Relative rejection threshold.
    pub refine_cutoff: f64,
    /// Whether sequence-derived correspondence was requested.
    pub align_seqs: bool,
    /// Atom subset name.
    pub atoms: String,
    /// Fitting residue selector.
    pub superpose_residues: String,
    /// Evaluation residue selector.
    pub rmsd_residues: String,
    /// Selected reference model serial.
    pub reference_model: usize,
    /// Selected query model serial.
    pub query_model: usize,
    /// Chain and residue correspondence; sequence results are absent in exact mode.
    pub chain_alignments: Vec<ChainAlignment>,
}

/// Calculate detailed RMSD with optional sequence correspondence and refinement.
///
/// Consumes structures to select conformers. Evaluation pairs never participate
/// in rejection unless they also belong to the fitting selection; even then they
/// remain in evaluation after rejection from fitting.
///
/// Selections use reference chain IDs and author residue numbers. With sequence
/// alignment enabled, corresponding query residues come from that alignment;
/// otherwise the same selectors apply to the query with exact atom identities.
/// Both selections independently default to all eligible residues.
/// See <https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md>.
pub fn get_rmsd(
    mut reference: PDB,
    mut query: PDB,
    options: &RmsdOptions,
) -> ArpeggiaResult<Analysis<RmsdResult>> {
    if !options.refine_cutoff.is_finite() || options.refine_cutoff <= 0.0 {
        return Err(ArpeggiaError::InvalidArgument(
            "refine_cutoff must be finite and positive".into(),
        ));
    }
    if !options.align_seqs && !options.chain_map.is_empty() {
        return Err(ArpeggiaError::InvalidArgument(
            "chain_map requires align_seqs".into(),
        ));
    }
    crate::seq_alignment::scoring(&options.alignment)?;
    let mut warnings = select_conformers(&mut reference);
    warnings.extend(select_conformers(&mut query));
    let superpose_selector = ResidueSelector::parse(&options.superpose_residues)?;
    let rmsd_selector = ResidueSelector::parse(&options.rmsd_residues)?;
    let reference_model = selected_model(&reference, options.model_num)?.serial_number();
    let query_model = selected_model(&query, options.model_num)?.serial_number();
    let (reference, query, chains) = if options.align_seqs {
        correspondence::aligned_coordinates(
            &reference,
            &query,
            options,
            &superpose_selector,
            &rmsd_selector,
            &mut warnings,
        )?
    } else {
        let first = select_coordinate_union(
            &reference,
            options.model_num,
            &superpose_selector,
            &rmsd_selector,
            options.atoms,
        )?;
        let second = select_coordinate_union(
            &query,
            options.model_num,
            &superpose_selector,
            &rmsd_selector,
            options.atoms,
        )?;
        validate_selection_correspondence(&first, &second)?;
        let chains = correspondence::exact_chains(&first.keys);
        (first, second, chains)
    };
    warnings.extend(reference.warnings.iter().cloned());
    warnings.extend(query.warnings.iter().cloned());
    let initial_count = reference.superpose_end;
    let eval_count = reference.coordinates.len() - reference.rmsd_start;
    let eval_reference = &reference.coordinates[reference.rmsd_start..];
    let eval_query = &query.coordinates[query.rmsd_start..];
    let mut retained: Vec<_> = (0..initial_count).collect();
    let mut first =
        prepare_coordinate_union(&reference.coordinates, initial_count, reference.rmsd_start)?;
    let mut second = prepare_coordinate_union(&query.coordinates, initial_count, query.rmsd_start)?;
    let initial_rmsd = kabsch_prepared_rmsd(&first, &second)?;
    let mut core_rmsd = initial_rmsd;
    let mut cycles = 0;
    for _ in 0..options.refine_cycles {
        if core_rmsd <= f64::EPSILON * first.scale.max(second.scale) * 64.0 {
            break;
        }
        let transform = fit_prepared_transform(&first, &second)?;
        let mut survivors = Vec::with_capacity(retained.len());
        for (position, index) in retained.iter().enumerate() {
            let residual = na::Vector3::from(first.points[position]) * transform.reference_factor
                - transform.rotation
                    * na::Vector3::from(second.points[position])
                    * transform.query_factor
                - transform.residual_centroid;
            // Compare normalized residuals to avoid overflow in cutoff * RMSD.
            let distance = residual.iter().fold(0.0_f64, |n, v| n.hypot(*v));
            let rms = transform.fit_residual_norm / (retained.len() as f64).sqrt();
            if distance / rms <= options.refine_cutoff {
                survivors.push(*index);
            }
        }
        cycles += 1;
        if survivors.len() == retained.len() {
            break;
        }
        let count = survivors.len();
        let pack = |coords: &[[f64; 3]], eval: &[[f64; 3]]| {
            let mut packed: Vec<_> = survivors.iter().map(|i| coords[*i]).collect();
            packed.extend_from_slice(eval);
            prepare_coordinate_union(&packed, count, count).map_err(|error| {
                ArpeggiaError::Calculation(format!(
                    "refinement retained {count} fitting pairs: {error}"
                ))
            })
        };
        first = pack(&reference.coordinates, eval_reference)?;
        second = pack(&query.coordinates, eval_query)?;
        retained = survivors;
        core_rmsd = kabsch_prepared_rmsd(&first, &second)?;
    }
    // Unchanged inspections leave fit and evaluation populations identical.
    // Preserve the core solver's reflection handling and rounding in that case.
    let rmsd = if retained.len() == initial_count && superpose_selector == rmsd_selector {
        core_rmsd
    } else {
        kabsch_prepared_selected_rmsd(&first, &second)?
    };
    let atoms = match options.atoms {
        AtomSubset::Ca => "ca",
        AtomSubset::Backbone => "backbone",
        AtomSubset::Heavy => "heavy",
        AtomSubset::All => "all",
    };
    Ok(Analysis::new(
        RmsdResult {
            rmsd,
            core_rmsd,
            initial_rmsd,
            initial_fit_atoms: initial_count,
            retained_fit_atoms: retained.len(),
            evaluation_atoms: eval_count,
            cycles,
            refine_cycles: options.refine_cycles,
            refine_cutoff: options.refine_cutoff,
            align_seqs: options.align_seqs,
            atoms: atoms.into(),
            superpose_residues: options.superpose_residues.clone(),
            rmsd_residues: options.rmsd_residues.clone(),
            reference_model,
            query_model,
            chain_alignments: chains,
        },
        warnings,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::WarningCode;
    #[test]
    fn identical_structure_has_zero_rmsd() {
        let input = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
        let pdb = crate::load_model(&input).unwrap().value;
        let analysis = get_rmsd(
            pdb.clone(),
            pdb,
            &RmsdOptions {
                superpose_residues: "A:1-20".into(),
                rmsd_residues: "A:1-20".into(),
                ..Default::default()
            },
        )
        .unwrap();
        assert!(analysis.value.rmsd < 1e-12);
    }

    #[test]
    fn public_rmsd_superposes_one_chain_and_scores_others() {
        let directory = std::env::temp_dir().join(format!(
            "arpeggia-relative-domain-rmsd-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&directory).unwrap();
        let reference_path = directory.join("reference.pdb");
        let query_path = directory.join("query.pdb");
        std::fs::write(
            &reference_path,
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
             ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00 20.00           C  \n\
             ATOM      3  CA  ALA A   3       0.000   1.000   0.000  1.00 20.00           C  \n\
             ATOM      4  CA  ALA B   1       0.000   0.000   1.000  1.00 20.00           C  \n\
             ATOM      5  CA  ALA C   1       0.000   0.000   2.000  1.00 20.00           C  \n\
             END\n",
        )
        .unwrap();
        std::fs::write(
            &query_path,
            "ATOM      1  CA  ALA A   1      10.000  -4.000   2.000  1.00 20.00           C  \n\
             ATOM      2  CA  ALA A   2      10.000  -3.000   2.000  1.00 20.00           C  \n\
             ATOM      3  CA  ALA A   3       9.000  -4.000   2.000  1.00 20.00           C  \n\
             ATOM      4  CA  ALA B   1      10.000  -4.000   5.000  1.00 20.00           C  \n\
             ATOM      5  CA  ALA C   1      10.000  -4.000   6.000  1.00 20.00           C  \n\
             END\n",
        )
        .unwrap();
        let reference = crate::load_model(reference_path.to_str().unwrap())
            .unwrap()
            .value;
        let query = crate::load_model(query_path.to_str().unwrap())
            .unwrap()
            .value;
        let selected = get_rmsd(
            reference.clone(),
            query.clone(),
            &RmsdOptions {
                superpose_residues: "A".into(),
                rmsd_residues: "B,C".into(),
                ..Default::default()
            },
        )
        .unwrap()
        .value;
        assert!((selected.rmsd - 2.0).abs() < 1e-12);

        let default_all = get_rmsd(
            reference,
            query,
            &RmsdOptions {
                superpose_residues: "A".into(),
                ..Default::default()
            },
        )
        .unwrap()
        .value;
        assert!((default_all.rmsd - (8.0_f64 / 5.0).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn public_rmsd_selects_conformers() {
        let input = std::env::temp_dir().join(format!(
            "arpeggia-rmsd-conformers-{}.pdb",
            std::process::id()
        ));
        std::fs::write(
            &input,
            "ATOM      1  CA AALA A   1       0.000   0.000   0.000  0.60 20.00           C  \n\
             ATOM      2  CA BALA A   1      10.000   0.000   0.000  0.40 20.00           C  \n\
             ATOM      3  CA  ALA A   2       1.000   0.000   0.000  1.00 20.00           C  \n\
             ATOM      4  CA  ALA A   3       0.000   1.000   0.000  1.00 20.00           C  \n\
             END\n",
        )
        .unwrap();
        let reference = pdbtbx::ReadOptions::default()
            .read(input.to_str().unwrap())
            .unwrap()
            .0;
        let query = reference.clone();
        let analysis = get_rmsd(reference, query, &RmsdOptions::default()).unwrap();
        assert_eq!(analysis.value.rmsd, 0.0);
        assert_eq!(
            analysis
                .warnings
                .iter()
                .filter(|warning| warning.code == WarningCode::ConformerSelected)
                .count(),
            2
        );
    }
}
