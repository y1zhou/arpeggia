//! Spatial Aggregation Propensity (SAP) score calculations.
//!
//! This module provides functions for calculating SAP scores, which predict
//! aggregation-prone regions in proteins based on hydrophobic surface patches.
//!
//! The SAP score was developed by Chennamsetty et al. and is defined in:
//! "Developability Index: A Rapid In Silico Tool for the Screening of
//! Antibody Aggregation Propensity" (J Pharm Sci, 2012).
//!
//! The formula is:
//! SAP(i) = Σ{j ∈ neighbors(i, R)} [ Hydrophobicity(j) × (SASA(j) / SASA_max(j)) ]
//!
//! Where:
//! - Neighbors are atoms/residues within radius R of atom/residue i
//! - Hydrophobicity uses the Black & Mould (1991) scale, normalized so glycine = 0
//! - SASA is the side-chain solvent accessible surface area
//! - SASA_max is the maximum SASA for that residue type

use crate::sasa::{get_atom_sasa, get_max_asa};
use pdbtbx::*;
use polars::prelude::*;
use std::collections::HashMap;

/// Black & Mould (1991) hydrophobicity scale values for amino acids.
/// These values are normalized so that the scale ranges from 0 to 1,
/// with Arginine being least hydrophobic (0.0) and Phenylalanine most hydrophobic (1.0).
///
/// For SAP calculations, we normalize this scale so that Glycine = 0.
/// This means hydrophobic residues (more hydrophobic than Gly) have positive values,
/// and hydrophilic residues (less hydrophobic than Gly) have negative values.
///
/// Original Black & Mould values (0-1 normalized):
/// - PHE: 1.000, ILE: 0.943, LEU: 0.943, TYR: 0.880, TRP: 0.878, VAL: 0.825
/// - MET: 0.738, PRO: 0.711, CYS: 0.680, ALA: 0.616
/// - GLY: 0.501, THR: 0.450, SER: 0.359, LYS: 0.283
/// - GLN: 0.251, ASN: 0.236, HIS: 0.165, GLU: 0.043
/// - ASP: 0.028, ARG: 0.000
///
/// We subtract 0.501 (Glycine's value) to get the normalized scale.
fn get_hydrophobicity(resn: &str) -> Option<f32> {
    // Black & Mould values minus Glycine (0.501)
    match resn.to_uppercase().as_str() {
        "ALA" => Some(0.616 - 0.501),  // 0.115
        "ARG" => Some(0.000 - 0.501),  // -0.501
        "ASN" => Some(0.236 - 0.501),  // -0.265
        "ASP" => Some(0.028 - 0.501),  // -0.473
        "CYS" => Some(0.680 - 0.501),  // 0.179
        "GLU" => Some(0.043 - 0.501),  // -0.458
        "GLN" => Some(0.251 - 0.501),  // -0.250
        "GLY" => Some(0.0),            // 0.0 (reference)
        "HIS" => Some(0.165 - 0.501),  // -0.336
        "ILE" => Some(0.943 - 0.501),  // 0.442
        "LEU" => Some(0.943 - 0.501),  // 0.442
        "LYS" => Some(0.283 - 0.501),  // -0.218
        "MET" => Some(0.738 - 0.501),  // 0.237
        "PHE" => Some(1.000 - 0.501),  // 0.499
        "PRO" => Some(0.711 - 0.501),  // 0.210
        "SER" => Some(0.359 - 0.501),  // -0.142
        "THR" => Some(0.450 - 0.501),  // -0.051
        "TRP" => Some(0.878 - 0.501),  // 0.377
        "TYR" => Some(0.880 - 0.501),  // 0.379
        "VAL" => Some(0.825 - 0.501),  // 0.324
        _ => None,
    }
}

/// Calculate the SAP score for each atom in a PDB structure.
///
/// The SAP score quantifies the aggregation propensity by combining the
/// solvent-accessible hydrophobic surface area of neighboring residues.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms for SASA calculation (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `sap_radius` - Radius in Ångströms for neighbor search (typically 5.0)
/// * `num_threads` - Number of threads for parallel processing
///
/// # Returns
///
/// A Polars DataFrame with columns:
/// - chain, resn, resi, insertion, atomn, atomi, sasa, sap_score
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_per_atom_sap_score};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
/// let sap_df = get_per_atom_sap_score(&pdb, 1.4, 100, 0, 5.0, 1);
/// println!("Calculated SAP score for {} atoms", sap_df.height());
/// ```
pub fn get_per_atom_sap_score(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    sap_radius: f32,
    num_threads: isize,
) -> DataFrame {
    // Use get_atom_sasa from sasa.rs for SASA calculation
    // Note: This already handles hydrogen stripping and solvent/ion removal
    let atom_sasa_df = get_atom_sasa(pdb, probe_radius, n_points, model_num, num_threads);

    // Create a lookup map from atom serial number to SASA value
    let atomi_col = atom_sasa_df.column("atomi").unwrap();
    let sasa_col = atom_sasa_df.column("sasa").unwrap();
    let atomi_values: Vec<i32> = atomi_col.i32().unwrap().into_iter().flatten().collect();
    let sasa_values: Vec<f32> = sasa_col.f32().unwrap().into_iter().flatten().collect();

    let sasa_map: HashMap<usize, f32> = atomi_values
        .iter()
        .zip(sasa_values.iter())
        .map(|(&atomi, &sasa)| (atomi as usize, sasa))
        .collect();

    // Also create a map of atom serial to residue name for hydrophobicity lookup
    let resn_col = atom_sasa_df.column("resn").unwrap();
    let resn_values: Vec<String> = resn_col
        .str()
        .unwrap()
        .into_iter()
        .flatten()
        .map(|s| s.to_string())
        .collect();

    let resn_map: HashMap<usize, String> = atomi_values
        .iter()
        .zip(resn_values.iter())
        .map(|(&atomi, resn)| (atomi as usize, resn.clone()))
        .collect();

    // Use pdbtbx's R-tree for spatial indexing (similar to InteractionComplex::get_atomic_contacts)
    let tree = pdb.create_hierarchy_rtree();
    let sap_radius_sq = (sap_radius * sap_radius) as f64;

    // Calculate SAP score for each atom in the DataFrame
    let sap_scores: Vec<f32> = atomi_values
        .iter()
        .map(|&atomi| {
            // Find the atom in the hierarchy to get its position
            let atom_hier = pdb
                .atoms_with_hierarchy()
                .find(|h| h.atom().serial_number() == atomi as usize);

            if let Some(hier) = atom_hier {
                // Find neighbors using pdbtbx's R-tree
                let neighbors = tree.locate_within_distance(hier.atom().pos(), sap_radius_sq);

                // Calculate SAP contribution from each neighbor
                let mut sap = 0.0f32;
                for neighbor in neighbors {
                    // Skip self
                    if neighbor.atom().serial_number() == atomi as usize {
                        continue;
                    }

                    let neighbor_serial = neighbor.atom().serial_number();

                    // Get SASA and residue name from our maps
                    if let (Some(&neighbor_sasa), Some(neighbor_resn)) =
                        (sasa_map.get(&neighbor_serial), resn_map.get(&neighbor_serial))
                    {
                        // Get hydrophobicity and max SASA for the neighbor's residue
                        // Using get_max_asa from sasa.rs as requested
                        if let (Some(hydrop), Some(max_sasa)) = (
                            get_hydrophobicity(neighbor_resn),
                            get_max_asa(neighbor_resn),
                        ) {
                            // Only consider side-chain atoms for SAP
                            // Use pdbtbx's is_backbone method
                            if !neighbor.atom().is_backbone() {
                                // SAP contribution = hydrophobicity * (SASA / max_SASA)
                                // Clamp to 1.0 because observed SASA can exceed theoretical max
                                // due to structural context and calculation parameters
                                let sasa_fraction = (neighbor_sasa / max_sasa).min(1.0);
                                sap += hydrop * sasa_fraction;
                            }
                        }
                    }
                }
                sap
            } else {
                0.0
            }
        })
        .collect();

    // Add SAP scores to the DataFrame
    atom_sasa_df
        .clone()
        .lazy()
        .with_column(Series::new("sap_score".into(), sap_scores).lit())
        .collect()
        .unwrap()
        .select(["chain", "resn", "resi", "insertion", "atomn", "atomi", "sasa", "sap_score"])
        .unwrap()
        .sort(["atomi"], Default::default())
        .unwrap()
}

/// Calculate the SAP score aggregated by residue.
///
/// The per-residue SAP score is calculated by summing the per-atom SAP scores
/// for all atoms in each residue, weighted by their SASA contribution.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms for SASA calculation (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `sap_radius` - Radius in Ångströms for neighbor search (typically 5.0)
/// * `num_threads` - Number of threads for parallel processing
///
/// # Returns
///
/// A Polars DataFrame with columns:
/// - chain, resn, resi, insertion, sasa, sap_score
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_per_residue_sap_score};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
/// let sap_df = get_per_residue_sap_score(&pdb, 1.4, 100, 0, 5.0, 1);
/// println!("Calculated SAP score for {} residues", sap_df.height());
/// ```
pub fn get_per_residue_sap_score(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    sap_radius: f32,
    num_threads: isize,
) -> DataFrame {
    // Get per-atom SAP scores
    let atom_sap = get_per_atom_sap_score(pdb, probe_radius, n_points, model_num, sap_radius, num_threads);

    // Aggregate by residue
    atom_sap
        .lazy()
        .group_by([col("chain"), col("resn"), col("resi"), col("insertion")])
        .agg([
            col("sasa").sum(),
            col("sap_score").sum(),
        ])
        .sort(
            ["chain", "resi", "insertion"],
            Default::default(),
        )
        .collect()
        .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::load_model;

    fn load_ubiquitin() -> PDB {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/1ubq.pdb");
        let (pdb, _) = load_model(&path);
        pdb
    }

    fn load_multi_chain() -> PDB {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/6bft.pdb");
        let (pdb, _) = load_model(&path);
        pdb
    }

    #[test]
    fn test_hydrophobicity_values() {
        // Glycine should be 0 (reference)
        assert!((get_hydrophobicity("GLY").unwrap() - 0.0).abs() < 1e-6);

        // Phenylalanine should be most hydrophobic (positive)
        let phe = get_hydrophobicity("PHE").unwrap();
        assert!(phe > 0.4, "PHE hydrophobicity should be high: {}", phe);

        // Arginine should be most hydrophilic (negative)
        let arg = get_hydrophobicity("ARG").unwrap();
        assert!(arg < -0.4, "ARG hydrophobicity should be low: {}", arg);

        // Unknown residues should return None
        assert!(get_hydrophobicity("XXX").is_none());
    }

    #[test]
    fn test_max_asa_values() {
        // All standard amino acids should have values (using get_max_asa from sasa.rs)
        let amino_acids = [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        ];

        for aa in amino_acids.iter() {
            let max_sasa = get_max_asa(aa);
            assert!(max_sasa.is_some(), "Should have max SASA for {}", aa);
            assert!(max_sasa.unwrap() > 0.0, "Max SASA for {} should be positive", aa);
        }

        // Larger residues should have larger max SASA
        let gly_sasa = get_max_asa("GLY").unwrap();
        let trp_sasa = get_max_asa("TRP").unwrap();
        assert!(trp_sasa > gly_sasa, "TRP should have larger max SASA than GLY");
    }

    #[test]
    fn test_per_atom_sap_returns_data() {
        let pdb = load_ubiquitin();
        let df = get_per_atom_sap_score(&pdb, 1.4, 100, 0, 5.0, 1);

        // Check that we get results
        assert!(!df.is_empty(), "SAP DataFrame should not be empty");

        // Check columns
        let columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(columns.contains(&"chain".to_string()));
        assert!(columns.contains(&"resn".to_string()));
        assert!(columns.contains(&"resi".to_string()));
        assert!(columns.contains(&"atomn".to_string()));
        assert!(columns.contains(&"atomi".to_string()));
        assert!(columns.contains(&"sasa".to_string()));
        assert!(columns.contains(&"sap_score".to_string()));
    }

    #[test]
    fn test_per_atom_sap_values_reasonable() {
        let pdb = load_ubiquitin();
        let df = get_per_atom_sap_score(&pdb, 1.4, 100, 0, 5.0, 1);

        // Get SAP scores
        let sap_values: Vec<f32> = df
            .column("sap_score")
            .unwrap()
            .f32()
            .unwrap()
            .into_iter()
            .flatten()
            .collect();

        // SAP values should have both positive and negative values
        // (hydrophobic and hydrophilic patches)
        let has_positive = sap_values.iter().any(|&v| v > 0.0);
        let has_negative = sap_values.iter().any(|&v| v < 0.0);

        assert!(has_positive, "Should have some positive SAP scores");
        assert!(has_negative, "Should have some negative SAP scores");
    }

    #[test]
    fn test_per_residue_sap_returns_data() {
        let pdb = load_ubiquitin();
        let df = get_per_residue_sap_score(&pdb, 1.4, 100, 0, 5.0, 1);

        // Check that we get results
        assert!(!df.is_empty(), "Residue SAP DataFrame should not be empty");

        // Check columns
        let columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(columns.contains(&"chain".to_string()));
        assert!(columns.contains(&"resn".to_string()));
        assert!(columns.contains(&"resi".to_string()));
        assert!(columns.contains(&"sasa".to_string()));
        assert!(columns.contains(&"sap_score".to_string()));
    }

    #[test]
    fn test_per_residue_sap_aggregation() {
        let pdb = load_ubiquitin();
        let atom_sap = get_per_atom_sap_score(&pdb, 1.4, 100, 0, 5.0, 1);
        let residue_sap = get_per_residue_sap_score(&pdb, 1.4, 100, 0, 5.0, 1);

        // Should have fewer rows in residue level
        assert!(
            residue_sap.height() < atom_sap.height(),
            "Residue SAP should have fewer rows: {} vs {}",
            residue_sap.height(),
            atom_sap.height()
        );

        // Ubiquitin has 76 residues
        assert!(
            residue_sap.height() >= 70 && residue_sap.height() <= 80,
            "Ubiquitin should have ~76 residues: {}",
            residue_sap.height()
        );
    }

    #[test]
    fn test_sap_radius_effect() {
        let pdb = load_ubiquitin();

        // Larger radius should include more neighbors, potentially higher absolute SAP
        let small_radius = get_per_residue_sap_score(&pdb, 1.4, 100, 0, 3.0, 1);
        let large_radius = get_per_residue_sap_score(&pdb, 1.4, 100, 0, 10.0, 1);

        // Sum of absolute SAP scores should be larger with larger radius
        let small_sum: f32 = small_radius
            .column("sap_score")
            .unwrap()
            .f32()
            .unwrap()
            .into_iter()
            .flatten()
            .map(|x| x.abs())
            .sum();

        let large_sum: f32 = large_radius
            .column("sap_score")
            .unwrap()
            .f32()
            .unwrap()
            .into_iter()
            .flatten()
            .map(|x| x.abs())
            .sum();

        assert!(
            large_sum > small_sum,
            "Larger radius should give larger absolute SAP sum: {} vs {}",
            large_sum,
            small_sum
        );
    }

    #[test]
    fn test_sap_multi_chain() {
        let pdb = load_multi_chain();
        let df = get_per_residue_sap_score(&pdb, 1.4, 100, 0, 5.0, 1);

        // Should have results
        assert!(!df.is_empty(), "Multi-chain SAP should not be empty");

        // Should have multiple chains
        let chain_count = df
            .column("chain")
            .unwrap()
            .unique()
            .unwrap()
            .len();
        assert!(chain_count > 1, "Should have multiple chains: {}", chain_count);
    }
}
