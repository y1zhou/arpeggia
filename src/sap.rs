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

use crate::sasa::{filter_pdb_by_model, prepare_pdb_for_sasa};
use pdbtbx::*;
use polars::prelude::*;
use rstar::{PointDistance, RTree, RTreeObject, AABB};

/// Black & Mould (1991) hydrophobicity scale values for amino acids.
/// These values are normalized so that the scale ranges from 0 to 1,
/// with Arginine being least hydrophobic (0.0) and Phenylalanine most hydrophobic (1.0).
///
/// For SAP calculations, we normalize this scale so that Glycine = 0.
/// This means hydrophobic residues (more hydrophobic than Gly) have positive values,
/// and hydrophilic residues (less hydrophobic than Gly) have negative values.
///
/// Original Black & Mould values (0-1 normalized):
/// - PHE: 1.000, TYR: 0.880, TRP: 0.878, VAL: 0.825
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

/// Maximum side-chain SASA values for amino acids in a Gly-X-Gly tripeptide context.
/// These values are used to normalize the side-chain SASA for SAP calculations.
///
/// Values are derived from Tien et al. (2013) "Maximum Allowed Solvent Accessibilities
/// of Residues in Proteins" PLOS ONE, specifically using side-chain values.
/// For Glycine and Alanine, we use small positive values since they have minimal
/// side chains but still contribute to the surface.
fn get_max_sidechain_sasa(resn: &str) -> Option<f32> {
    // Side-chain SASA values (approximate, based on literature)
    // These represent the maximum exposed side-chain surface area
    match resn.to_uppercase().as_str() {
        "ALA" => Some(67.0),   // Small side chain (CB only)
        "ARG" => Some(196.0),
        "ASN" => Some(113.0),
        "ASP" => Some(106.0),
        "CYS" => Some(104.0),
        "GLU" => Some(138.0),
        "GLN" => Some(144.0),
        "GLY" => Some(25.0),   // No side chain, but CA contributes
        "HIS" => Some(151.0),
        "ILE" => Some(140.0),
        "LEU" => Some(137.0),
        "LYS" => Some(167.0),
        "MET" => Some(160.0),
        "PHE" => Some(175.0),
        "PRO" => Some(105.0),
        "SER" => Some(80.0),
        "THR" => Some(102.0),
        "TRP" => Some(217.0),
        "TYR" => Some(187.0),
        "VAL" => Some(117.0),
        _ => None,
    }
}

/// Helper struct for R-tree spatial indexing of atoms
#[derive(Debug, Clone)]
struct SpatialAtom {
    /// Atom position
    position: [f32; 3],
    /// Atom serial number
    atom_serial: usize,
    /// Residue name (3-letter code)
    resn: String,
    /// Is this a side-chain atom (not backbone N, CA, C, O)?
    is_sidechain: bool,
    /// SASA value for this atom
    sasa: f32,
}

impl RTreeObject for SpatialAtom {
    type Envelope = AABB<[f32; 3]>;

    fn envelope(&self) -> Self::Envelope {
        AABB::from_point(self.position)
    }
}

impl PointDistance for SpatialAtom {
    fn distance_2(&self, point: &[f32; 3]) -> f32 {
        let dx = self.position[0] - point[0];
        let dy = self.position[1] - point[1];
        let dz = self.position[2] - point[2];
        dx * dx + dy * dy + dz * dz
    }
}

/// Check if an atom is a backbone atom
fn is_backbone_atom(atom_name: &str) -> bool {
    matches!(atom_name, "N" | "CA" | "C" | "O" | "OXT" | "H" | "HA" | "HA2" | "HA3")
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
    use rust_sasa::{Atom as SASAAtom, calculate_sasa_internal};

    // Prepare PDB: remove solvent, ions, and hydrogens, then filter by model
    let pdb_prepared = prepare_pdb_for_sasa(pdb);
    let pdb_filtered = filter_pdb_by_model(&pdb_prepared, model_num);

    // Get the actual model number we're using
    let actual_model_num = if model_num == 0 {
        pdb_filtered
            .models()
            .collect::<Vec<_>>()
            .first()
            .map_or(0, |m| m.serial_number())
    } else {
        model_num
    };

    // Collect atoms with hierarchy info for SASA calculation
    let atoms_with_hier: Vec<_> = pdb_filtered
        .atoms_with_hierarchy()
        .filter(|x| x.model().serial_number() == actual_model_num)
        .collect();

    // Calculate SASA for each atom
    let sasa_atoms: Vec<SASAAtom> = atoms_with_hier
        .iter()
        .map(|x| SASAAtom {
            position: [
                x.atom().pos().0 as f32,
                x.atom().pos().1 as f32,
                x.atom().pos().2 as f32,
            ],
            radius: x
                .atom()
                .element()
                .unwrap()
                .atomic_radius()
                .van_der_waals
                .unwrap() as f32,
            id: x.atom().serial_number(),
            parent_id: None,
        })
        .collect();

    let atom_sasa_values = calculate_sasa_internal(&sasa_atoms, probe_radius, n_points, num_threads);

    // Build SpatialAtom list for R-tree
    let spatial_atoms: Vec<SpatialAtom> = atoms_with_hier
        .iter()
        .zip(atom_sasa_values.iter())
        .map(|(hier, &sasa)| {
            let atom_name = hier.atom().name().to_string();
            SpatialAtom {
                position: [
                    hier.atom().pos().0 as f32,
                    hier.atom().pos().1 as f32,
                    hier.atom().pos().2 as f32,
                ],
                atom_serial: hier.atom().serial_number(),
                resn: hier.residue().name().unwrap_or("").to_string(),
                is_sidechain: !is_backbone_atom(&atom_name),
                sasa,
            }
        })
        .collect();

    // Build R-tree for spatial indexing
    let rtree = RTree::bulk_load(spatial_atoms.clone());

    // Calculate SAP score for each atom
    let sap_radius_sq = sap_radius * sap_radius;
    let sap_scores: Vec<f32> = spatial_atoms
        .iter()
        .map(|atom| {
            // Find neighbors within sap_radius
            let neighbors: Vec<_> = rtree
                .locate_within_distance(atom.position, sap_radius_sq)
                .collect();

            // Calculate SAP contribution from each neighbor
            let mut sap = 0.0f32;
            for neighbor in neighbors {
                // Skip self
                if neighbor.atom_serial == atom.atom_serial {
                    continue;
                }

                // Get hydrophobicity and max SASA for the neighbor's residue
                if let (Some(hydrop), Some(max_sasa)) = (
                    get_hydrophobicity(&neighbor.resn),
                    get_max_sidechain_sasa(&neighbor.resn),
                ) {
                    // Only consider side-chain atoms for SAP
                    if neighbor.is_sidechain && max_sasa > 0.0 {
                        // SAP contribution = hydrophobicity * (SASA / max_SASA)
                        let sasa_fraction = (neighbor.sasa / max_sasa).min(1.0);
                        sap += hydrop * sasa_fraction;
                    }
                }
            }
            sap
        })
        .collect();

    // Build result DataFrame
    let chains: Vec<String> = atoms_with_hier
        .iter()
        .map(|x| x.chain().id().to_string())
        .collect();
    let resns: Vec<String> = atoms_with_hier
        .iter()
        .map(|x| x.residue().name().unwrap_or("").to_string())
        .collect();
    let resis: Vec<i32> = atoms_with_hier
        .iter()
        .map(|x| x.residue().id().0 as i32)
        .collect();
    let insertions: Vec<String> = atoms_with_hier
        .iter()
        .map(|x| x.residue().id().1.unwrap_or("").to_string())
        .collect();
    let atomns: Vec<String> = atoms_with_hier
        .iter()
        .map(|x| x.atom().name().to_string())
        .collect();
    let atomis: Vec<i32> = atoms_with_hier
        .iter()
        .map(|x| x.atom().serial_number() as i32)
        .collect();

    df!(
        "chain" => chains,
        "resn" => resns,
        "resi" => resis,
        "insertion" => insertions,
        "atomn" => atomns,
        "atomi" => atomis,
        "sasa" => atom_sasa_values,
        "sap_score" => sap_scores,
    )
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
    fn test_max_sidechain_sasa() {
        // All standard amino acids should have values
        let amino_acids = [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        ];

        for aa in amino_acids.iter() {
            let max_sasa = get_max_sidechain_sasa(aa);
            assert!(max_sasa.is_some(), "Should have max SASA for {}", aa);
            assert!(max_sasa.unwrap() > 0.0, "Max SASA for {} should be positive", aa);
        }

        // Larger residues should have larger max SASA
        let gly_sasa = get_max_sidechain_sasa("GLY").unwrap();
        let trp_sasa = get_max_sidechain_sasa("TRP").unwrap();
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
