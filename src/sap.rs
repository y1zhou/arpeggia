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
//! - Neighbors are atoms within radius R of atom i
//! - Hydrophobicity uses the Black & Mould (1991) scale, normalized so glycine = 0
//! - SASA is the side-chain solvent accessible surface area
//! - SASA_max is the maximum SASA for that residue type

use crate::sasa::{calculate_sap_atom_sasa_records, validate_sasa_input};
use crate::structure::prepare_structure;
use pdbtbx::*;
use polars::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

#[derive(Clone, Debug)]
struct AtomSapRecord {
    chain: String,
    resn: String,
    resi: i64,
    insertion: String,
    atomn: String,
    atomi: u64,
    sasa: f32,
    sap_score: f32,
}

#[derive(Clone, Debug)]
struct ResidueSapRecord {
    chain: String,
    resn: String,
    resi: i64,
    insertion: String,
    sc_sasa: f32,
    sap_score: f32,
    max_sc_asa: f32,
    relative_sc_sasa: f32,
}

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
        "ALA" => Some(0.616 - 0.501), // 0.115
        "ARG" => Some(0.000 - 0.501), // -0.501
        "ASN" => Some(0.236 - 0.501), // -0.265
        "ASP" => Some(0.028 - 0.501), // -0.473
        "CYS" => Some(0.680 - 0.501), // 0.179
        "GLU" => Some(0.043 - 0.501), // -0.458
        "GLN" => Some(0.251 - 0.501), // -0.250
        "GLY" => Some(0.000),         // 0.0 (reference)
        "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP" => Some(0.165 - 0.501), // -0.336
        "ILE" | "LEU" => Some(0.943 - 0.501), // 0.442
        "LYS" => Some(0.283 - 0.501), // -0.218
        "MET" | "MSE" => Some(0.738 - 0.501), // 0.237
        "PHE" => Some(1.000 - 0.501), // 0.499
        "PRO" => Some(0.711 - 0.501), // 0.210
        "SER" => Some(0.359 - 0.501), // -0.142
        "THR" => Some(0.450 - 0.501), // -0.051
        "TRP" => Some(0.878 - 0.501), // 0.377
        "TYR" => Some(0.880 - 0.501), // 0.379
        "VAL" => Some(0.825 - 0.501), // 0.324
        _ => None,
    }
}

/// Get the maximum side-chain solvent accessible surface area (SASA).
///
/// These values are adapted from the Rosetta database file
/// `database/scoring/score_functions/sap_sasa_calib.dat` and represent
/// the maximum SASA contribution from side-chain atoms only.
///
/// Note: Glycine uses the H atom contribution as it has no side-chain.
///
/// Returns the maximum side-chain SASA in Å² for a given 3-letter amino acid code,
/// or None if the residue is not a standard amino acid.
fn get_sc_max_asa(resn: &str) -> Option<f32> {
    match resn.to_uppercase().as_str() {
        "ALA" => Some(15.395),
        "ARG" => Some(124.338),
        "ASN" => Some(90.303),
        "ASP" => Some(87.601),
        "CYS" => Some(46.456),
        "GLN" => Some(99.186),
        "GLY" => Some(3.229), // Using H
        "GLU" => Some(95.534),
        "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP" => Some(96.532),
        "ILE" => Some(31.448),
        "LEU" => Some(30.271),
        "LYS" => Some(61.962),
        "MET" | "MSE" => Some(65.233),
        "PHE" => Some(67.945),
        "PRO" => Some(17.812),
        "SER" => Some(39.355),
        "THR" => Some(42.648),
        "TRP" => Some(101.491),
        "TYR" => Some(94.478),
        "VAL" => Some(26.702),
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
/// * `chains` - Comma-separated chain IDs to include (e.g., "A,B,C"). Empty string includes all chains.
///
/// # Returns
///
/// A Polars `DataFrame` with columns:
/// - `chain`, `resn`, `resi`, `insertion`, `atomn`, `atomi`, `sasa`, `sap_score`
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_per_atom_sap_score};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let pdb = load_model(&input_file).unwrap().value;
///
/// // Calculate SAP for all chains
/// let sap_df = get_per_atom_sap_score(&pdb, 1.4, 100, 0, 5.0, "").unwrap().value;
/// println!("Calculated SAP score for {} atoms", sap_df.height());
///
/// // Calculate SAP for only chains H and L (antibody heavy and light chain)
/// let sap_hl = get_per_atom_sap_score(&pdb, 1.4, 100, 0, 5.0, "H,L").unwrap().value;
/// ```
pub fn get_per_atom_sap_score(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    sap_radius: f32,
    chains: &str,
) -> crate::ArpeggiaResult<crate::Analysis<DataFrame>> {
    if !sap_radius.is_finite() || sap_radius <= 0.0 {
        return Err(crate::ArpeggiaError::InvalidArgument(
            "sap_radius must be finite and positive".into(),
        ));
    }
    let analysis = validate_sasa_input(pdb, probe_radius, n_points, model_num, chains)?;
    let (records, unsupported) = calculate_per_atom_sap_records(
        &analysis.value,
        probe_radius,
        n_points,
        model_num,
        sap_radius,
        chains,
    )?;
    let mut warnings = analysis.warnings;
    append_unsupported_monomer_warning(&unsupported, &mut warnings);
    Ok(crate::Analysis::new(
        atom_sap_records_to_dataframe(&records),
        warnings,
    ))
}

fn calculate_per_atom_sap_records(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    sap_radius: f32,
    chains: &str,
) -> crate::ArpeggiaResult<(Vec<AtomSapRecord>, BTreeSet<String>)> {
    // [ACCEPTED DEVIATION] SAP retains elemental van der Waals radii to preserve
    // its Rosetta-calibrated exposure scale; standard SASA uses ProtOr radii.
    let atom_sasa_records =
        calculate_sap_atom_sasa_records(pdb, probe_radius, n_points, model_num, chains)?;
    let unsupported = atom_sasa_records
        .iter()
        .filter(|record| {
            get_hydrophobicity(&record.resn).is_none() || get_sc_max_asa(&record.resn).is_none()
        })
        .map(|record| record.resn.clone())
        .collect::<BTreeSet<_>>();
    // [WARNING] Rosetta publishes no SAP calibration for these monomers. They
    // are omitted and diagnosed instead of borrowing an unjustified constant.
    let atom_sasa_records = atom_sasa_records
        .into_iter()
        .filter(|record| !unsupported.contains(&record.resn))
        .collect::<Vec<_>>();
    let atom_sasa_map: HashMap<u64, f32> = atom_sasa_records
        .iter()
        .map(|record| (record.atomi, record.sasa))
        .collect();

    let pdb_no_hydrogens = prepare_structure(pdb, model_num, true, chains);
    let tree = pdb_no_hydrogens.create_hierarchy_rtree();
    let sap_radius_sq = f64::from(sap_radius * sap_radius);

    let sap_scores_map: HashMap<u64, f32> = pdb_no_hydrogens
        .atoms_with_hierarchy()
        .filter(|h| h.is_sidechain())
        .map(|x| {
            let x_atomi = x.atom().serial_number() as u64;
            let atom_sap_score = tree
                .locate_within_distance(x.atom().pos(), sap_radius_sq)
                .filter(|y| y.is_sidechain())
                .filter_map(|y| {
                    let neighbor_resn = y.residue().name()?;
                    let neighbor_atom_sasa =
                        atom_sasa_map.get(&(y.atom().serial_number() as u64))?;
                    let neighbor_res_hydrophobicity = get_hydrophobicity(neighbor_resn)?;
                    let max_res_asa = get_sc_max_asa(neighbor_resn)?;
                    Some(
                        neighbor_res_hydrophobicity
                            * (neighbor_atom_sasa / max_res_asa).clamp(0.0, 1.0),
                    )
                })
                .sum::<f32>();
            (x_atomi, atom_sap_score)
        })
        .collect();

    let sidechain_atomi: HashSet<u64> = pdb_no_hydrogens
        .atoms_with_hierarchy()
        .filter(|h| h.is_sidechain())
        .map(|h| h.atom().serial_number() as u64)
        .collect();

    let mut records = atom_sasa_records
        .into_iter()
        .filter(|record| sidechain_atomi.contains(&record.atomi))
        .map(|record| AtomSapRecord {
            sap_score: sap_scores_map.get(&record.atomi).copied().unwrap_or(0.0),
            chain: record.chain,
            resn: record.resn,
            resi: record.resi,
            insertion: record.insertion,
            atomn: record.atomn,
            atomi: record.atomi,
            sasa: record.sasa,
        })
        .collect::<Vec<_>>();
    records.sort_by_key(|record| record.atomi);
    Ok((records, unsupported))
}

fn append_unsupported_monomer_warning(
    unsupported: &BTreeSet<String>,
    warnings: &mut Vec<crate::AnalysisWarning>,
) {
    if !unsupported.is_empty() {
        warnings.push(crate::AnalysisWarning::new(
            crate::WarningCode::UnsupportedMonomer,
            format!(
                "SAP omitted monomers without Rosetta calibration: {}",
                unsupported.iter().cloned().collect::<Vec<_>>().join(", ")
            ),
        ));
    }
}

fn atom_sap_records_to_dataframe(records: &[AtomSapRecord]) -> DataFrame {
    df!(
        "chain" => records.iter().map(|r| r.chain.clone()).collect::<Vec<String>>(),
        "resn" => records.iter().map(|r| r.resn.clone()).collect::<Vec<String>>(),
        "resi" => records.iter().map(|r| r.resi).collect::<Vec<i64>>(),
        "insertion" => records.iter().map(|r| r.insertion.clone()).collect::<Vec<String>>(),
        "atomn" => records.iter().map(|r| r.atomn.clone()).collect::<Vec<String>>(),
        "atomi" => records.iter().map(|r| r.atomi).collect::<Vec<u64>>(),
        "sasa" => records.iter().map(|r| r.sasa).collect::<Vec<f32>>(),
        "sap_score" => records.iter().map(|r| r.sap_score).collect::<Vec<f32>>(),
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
/// * `chains` - Comma-separated chain IDs to include (e.g., "A,B,C"). Empty string includes all chains.
///
/// # Returns
///
/// A Polars `DataFrame` with columns:
/// - `chain`, `resn`, `resi`, `insertion`, `sc_sasa`, `sap_score`
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_per_residue_sap_score};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let pdb = load_model(&input_file).unwrap().value;
///
/// // Calculate SAP for all chains
/// let sap_df = get_per_residue_sap_score(&pdb, 1.4, 100, 0, 5.0, "").unwrap().value;
/// println!("Calculated SAP score for {} residues", sap_df.height());
///
/// // Calculate SAP for only chain A
/// let sap_a = get_per_residue_sap_score(&pdb, 1.4, 100, 0, 5.0, "A").unwrap().value;
/// ```
pub fn get_per_residue_sap_score(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    sap_radius: f32,
    chains: &str,
) -> crate::ArpeggiaResult<crate::Analysis<DataFrame>> {
    if !sap_radius.is_finite() || sap_radius <= 0.0 {
        return Err(crate::ArpeggiaError::InvalidArgument(
            "sap_radius must be finite and positive".into(),
        ));
    }
    let analysis = validate_sasa_input(pdb, probe_radius, n_points, model_num, chains)?;
    let (records, unsupported) = calculate_per_residue_sap_records(
        &analysis.value,
        probe_radius,
        n_points,
        model_num,
        sap_radius,
        chains,
    )?;
    let mut warnings = analysis.warnings;
    append_unsupported_monomer_warning(&unsupported, &mut warnings);
    Ok(crate::Analysis::new(
        residue_sap_records_to_dataframe(&records),
        warnings,
    ))
}

fn calculate_per_residue_sap_records(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    sap_radius: f32,
    chains: &str,
) -> crate::ArpeggiaResult<(Vec<ResidueSapRecord>, BTreeSet<String>)> {
    let mut grouped: BTreeMap<(String, String, i64, String), (f32, f32)> = BTreeMap::new();

    let (atom_records, unsupported) =
        calculate_per_atom_sap_records(pdb, probe_radius, n_points, model_num, sap_radius, chains)?;
    for record in atom_records {
        let key = (record.chain, record.resn, record.resi, record.insertion);
        let (sc_sasa, sap_score) = grouped.entry(key).or_insert((0.0, 0.0));
        *sc_sasa += record.sasa;
        // [ACCEPTED DEVIATION] Match Rosetta's positive-only per-residue SAP
        // accumulation while keeping the complete side-chain SASA above.
        *sap_score += record.sap_score.max(0.0);
    }

    let mut records = grouped
        .into_iter()
        .map(|((chain, resn, resi, insertion), (sc_sasa, sap_score))| {
            let max_sc_asa = get_sc_max_asa(&resn).unwrap();
            ResidueSapRecord {
                chain,
                resn,
                resi,
                insertion,
                sc_sasa,
                sap_score,
                max_sc_asa,
                relative_sc_sasa: (sc_sasa / max_sc_asa).clamp(0.0, 1.0),
            }
        })
        .collect::<Vec<_>>();
    records.sort_by(|a, b| (&a.chain, a.resi, &a.insertion).cmp(&(&b.chain, b.resi, &b.insertion)));
    Ok((records, unsupported))
}

fn residue_sap_records_to_dataframe(records: &[ResidueSapRecord]) -> DataFrame {
    df!(
        "chain" => records.iter().map(|r| r.chain.clone()).collect::<Vec<String>>(),
        "resn" => records.iter().map(|r| r.resn.clone()).collect::<Vec<String>>(),
        "resi" => records.iter().map(|r| r.resi).collect::<Vec<i64>>(),
        "insertion" => records.iter().map(|r| r.insertion.clone()).collect::<Vec<String>>(),
        "sc_sasa" => records.iter().map(|r| r.sc_sasa).collect::<Vec<f32>>(),
        "sap_score" => records.iter().map(|r| r.sap_score).collect::<Vec<f32>>(),
        "max_sc_asa" => records.iter().map(|r| r.max_sc_asa).collect::<Vec<f32>>(),
        "relative_sc_sasa" => records.iter().map(|r| r.relative_sc_sasa).collect::<Vec<f32>>(),
    )
    .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{WarningCode, load_model, run_with_threads};
    use std::io::BufReader;

    fn load_ubiquitin() -> PDB {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/1ubq.pdb");
        load_model(&path).unwrap().value
    }

    fn load_multi_chain() -> PDB {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/6bft.pdb");
        load_model(&path).unwrap().value
    }

    #[test]
    fn test_hydrophobicity_values() {
        // Glycine should be 0 (reference)
        assert!((get_hydrophobicity("GLY").unwrap() - 0.0).abs() < 1e-6);

        // Phenylalanine should be most hydrophobic (positive)
        let phe = get_hydrophobicity("PHE").unwrap();
        assert!(phe > 0.4, "PHE hydrophobicity should be high: {phe}");

        // Arginine should be most hydrophilic (negative)
        let arg = get_hydrophobicity("ARG").unwrap();
        assert!(arg < -0.4, "ARG hydrophobicity should be low: {arg}");

        // Unknown residues should return None
        assert!(get_hydrophobicity("XXX").is_none());

        for alias in ["HID", "HIE", "HIP", "HSD", "HSE", "HSP"] {
            assert_eq!(get_hydrophobicity(alias), get_hydrophobicity("HIS"));
            assert_eq!(get_sc_max_asa(alias), get_sc_max_asa("HIS"));
        }
        assert_eq!(get_hydrophobicity("MSE"), get_hydrophobicity("MET"));
        assert_eq!(get_sc_max_asa("MSE"), get_sc_max_asa("MET"));
    }

    #[test]
    fn test_max_asa_values() {
        // All standard amino acids should have values (using get_max_asa from sasa.rs)
        let amino_acids = [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        ];

        for aa in &amino_acids {
            let max_sasa = get_sc_max_asa(aa);
            assert!(max_sasa.is_some(), "Should have max SASA for {aa}");
            assert!(
                max_sasa.unwrap() > 0.0,
                "Max SASA for {aa} should be positive"
            );
        }

        // Larger residues should have larger max SASA
        let gly_sasa = get_sc_max_asa("GLY").unwrap();
        let trp_sasa = get_sc_max_asa("TRP").unwrap();
        assert!(
            trp_sasa > gly_sasa,
            "TRP should have larger max SASA than GLY"
        );
    }

    #[test]
    fn test_per_atom_sap_returns_data() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_per_atom_sap_score(&pdb, 1.4, 100, 0, 5.0, ""));
        let df = df.unwrap().value;

        // Check that we get results
        assert!(df.height() > 0, "SAP DataFrame should not be empty");

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
        assert_eq!(df.column("resi").unwrap().dtype(), &DataType::Int64);
        assert_eq!(df.column("atomi").unwrap().dtype(), &DataType::UInt64);
    }

    #[test]
    fn test_per_atom_sap_values_reasonable() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_per_atom_sap_score(&pdb, 1.4, 100, 0, 5.0, ""));
        let df = df.unwrap().value;

        // Get SAP scores
        let sap_values: Vec<f32> = df
            .column("sap_score")
            .unwrap()
            .f32()
            .unwrap()
            .iter()
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
    fn unsupported_sap_monomers_are_omitted_with_a_warning() {
        let input =
            b"ATOM      1  CB  SEC A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
END                                                                             \n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();

        let analysis = get_per_residue_sap_score(&pdb, 1.4, 100, 0, 5.0, "").unwrap();
        assert_eq!(analysis.value.height(), 0);
        assert!(
            analysis
                .warnings
                .iter()
                .any(|warning| warning.code == WarningCode::UnsupportedMonomer)
        );
    }

    #[test]
    fn test_per_residue_sap_returns_data() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_per_residue_sap_score(&pdb, 1.4, 100, 0, 5.0, ""));
        let df = df.unwrap().value;

        // Check that we get results
        assert!(df.height() > 0, "Residue SAP DataFrame should not be empty");

        // Check columns
        let columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(columns.contains(&"chain".to_string()));
        assert!(columns.contains(&"resn".to_string()));
        assert!(columns.contains(&"resi".to_string()));
        assert!(columns.contains(&"sc_sasa".to_string()));
        assert!(columns.contains(&"sap_score".to_string()));
    }

    #[test]
    fn residue_sap_keeps_complete_sidechain_sasa() {
        let pdb = load_ubiquitin();
        let atom = run_with_threads(1, || get_per_atom_sap_score(&pdb, 1.4, 100, 0, 5.0, ""));
        let atom = atom.unwrap().value;
        let residue = run_with_threads(1, || get_per_residue_sap_score(&pdb, 1.4, 100, 0, 5.0, ""));
        let residue = residue.unwrap().value;
        let mut expected = BTreeMap::<(String, String, i64, String), (f32, f32)>::new();

        for row in 0..atom.height() {
            let key = (
                atom.column("chain")
                    .unwrap()
                    .str()
                    .unwrap()
                    .get(row)
                    .unwrap()
                    .to_string(),
                atom.column("resn")
                    .unwrap()
                    .str()
                    .unwrap()
                    .get(row)
                    .unwrap()
                    .to_string(),
                atom.column("resi")
                    .unwrap()
                    .i64()
                    .unwrap()
                    .get(row)
                    .unwrap(),
                atom.column("insertion")
                    .unwrap()
                    .str()
                    .unwrap()
                    .get(row)
                    .unwrap()
                    .to_string(),
            );
            let entry = expected.entry(key).or_default();
            entry.0 += atom
                .column("sasa")
                .unwrap()
                .f32()
                .unwrap()
                .get(row)
                .unwrap();
            entry.1 += atom
                .column("sap_score")
                .unwrap()
                .f32()
                .unwrap()
                .get(row)
                .unwrap()
                .max(0.0);
        }

        assert_eq!(residue.height(), expected.len());
        for row in 0..residue.height() {
            let key = (
                residue
                    .column("chain")
                    .unwrap()
                    .str()
                    .unwrap()
                    .get(row)
                    .unwrap()
                    .to_string(),
                residue
                    .column("resn")
                    .unwrap()
                    .str()
                    .unwrap()
                    .get(row)
                    .unwrap()
                    .to_string(),
                residue
                    .column("resi")
                    .unwrap()
                    .i64()
                    .unwrap()
                    .get(row)
                    .unwrap(),
                residue
                    .column("insertion")
                    .unwrap()
                    .str()
                    .unwrap()
                    .get(row)
                    .unwrap()
                    .to_string(),
            );
            let (expected_sasa, expected_score) = expected[&key];
            let actual_sasa = residue
                .column("sc_sasa")
                .unwrap()
                .f32()
                .unwrap()
                .get(row)
                .unwrap();
            let actual_score = residue
                .column("sap_score")
                .unwrap()
                .f32()
                .unwrap()
                .get(row)
                .unwrap();
            assert!(
                (actual_sasa - expected_sasa).abs() < 1e-4,
                "wrong SASA for {key:?}"
            );
            assert!(
                (actual_score - expected_score).abs() < 1e-4,
                "wrong SAP for {key:?}"
            );
        }
    }

    #[test]
    fn test_sap_multi_chain() {
        let pdb = load_multi_chain();
        let df = run_with_threads(1, || get_per_residue_sap_score(&pdb, 1.4, 100, 0, 5.0, ""));
        let df = df.unwrap().value;

        // Should have results
        assert!(df.height() > 0, "Multi-chain SAP should not be empty");

        // Should have multiple chains
        let chain_count = df.column("chain").unwrap().unique().unwrap().len();
        assert!(
            chain_count > 1,
            "Should have multiple chains: {chain_count}"
        );
    }
}
