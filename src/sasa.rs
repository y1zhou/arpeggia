//! Solvent Accessible Surface Area (SASA) calculations.
//!
//! This module provides functions for calculating SASA at different levels
//! of granularity (atom, residue, chain) and related metrics like dSASA
//! (buried surface area) and relative SASA.

use crate::contacts::InteractingEntity;
use crate::structure::{StructurePreparation, parse_groups, prepare_structure};
use pdbtbx::*;
use polars::prelude::*;
use std::collections::{BTreeMap, HashSet};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum AtomPolarity {
    Polar,
    Hydrophobic,
    Unclassified,
}

impl AtomPolarity {
    fn as_str(self) -> &'static str {
        match self {
            Self::Polar => "polar",
            Self::Hydrophobic => "hydrophobic",
            Self::Unclassified => "unclassified",
        }
    }
}

#[derive(Clone, Copy)]
enum RadiusScheme {
    StandardProtor,
    SapReduce,
}

#[derive(Clone, Debug)]
pub(crate) struct AtomSasaRecord {
    pub(crate) chain: String,
    pub(crate) resn: String,
    pub(crate) resi: i64,
    pub(crate) insertion: String,
    pub(crate) altloc: String,
    pub(crate) atomn: String,
    pub(crate) atomi: u64,
    pub(crate) sasa: f32,
    pub(crate) polarity: AtomPolarity,
}

#[derive(Clone, Debug)]
pub(crate) struct ResidueSasaRecord {
    pub(crate) chain: String,
    pub(crate) resn: String,
    pub(crate) resi: i64,
    pub(crate) insertion: String,
    pub(crate) sasa: f32,
    pub(crate) polar_sasa: f32,
    pub(crate) hydrophobic_sasa: f32,
    pub(crate) unclassified_sasa: f32,
}

#[derive(Clone, Debug)]
pub(crate) struct ChainSasaRecord {
    pub(crate) chain: String,
    pub(crate) sasa: f32,
    pub(crate) polar_sasa: f32,
    pub(crate) hydrophobic_sasa: f32,
    pub(crate) unclassified_sasa: f32,
}

#[derive(Clone, Debug)]
pub(crate) struct RelativeSasaRecord {
    pub(crate) chain: String,
    pub(crate) resn: String,
    pub(crate) resi: i64,
    pub(crate) insertion: String,
    pub(crate) sasa: f32,
    pub(crate) polar_sasa: f32,
    pub(crate) hydrophobic_sasa: f32,
    pub(crate) unclassified_sasa: f32,
    pub(crate) relative_sasa: Option<f32>,
}

/// Two-sided buried area and its Rosetta `SasaFilter` polarity partition.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DsasaResult {
    /// Two-sided buried surface area.
    pub dsasa: f32,
    /// Buried area assigned to polar atoms.
    pub polar_dsasa: f32,
    /// Buried area assigned to hydrophobic atoms.
    pub hydrophobic_dsasa: f32,
    /// Buried area whose atom polarity is unsupported.
    pub unclassified_dsasa: f32,
}

/// Calculate solvent accessible surface area (SASA) for each atom in a PDB structure.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `remove_hydrogens` - Whether to remove hydrogen atoms before calculation
/// * `chains` - Comma-separated chain IDs to include (e.g., "A,B,C"). Empty string includes all chains.
///
/// # Returns
///
/// A Polars `DataFrame` with columns:
/// - `atomi`, `sasa`
/// - `chain`, `resn`, `resi`, `insertion`, `altloc`, `atomn`
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_atom_sasa};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let pdb = load_model(&input_file).unwrap().value;
///
/// // Calculate SASA for all chains
/// let sasa_df = get_atom_sasa(&pdb, 1.4, 100, 0, true, "").unwrap().value;
/// println!("Calculated SASA for {} atoms", sasa_df.height());
///
/// // Calculate SASA for only chains A and B
/// let sasa_ab = get_atom_sasa(&pdb, 1.4, 100, 0, true, "A,B").unwrap().value;
/// ```
pub fn get_atom_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    remove_hydrogens: bool,
    chains: &str,
) -> crate::ArpeggiaResult<crate::Analysis<DataFrame>> {
    let analysis = validate_sasa_input(pdb, probe_radius, n_points, model_num, chains)?;
    let records = calculate_atom_sasa_records(
        &analysis.value,
        probe_radius,
        n_points,
        model_num,
        remove_hydrogens,
        chains,
    )?;
    let mut warnings = analysis.warnings;
    append_polarity_warning(&records, &mut warnings);
    Ok(crate::Analysis::new(
        atom_sasa_records_to_dataframe(&records),
        warnings,
    ))
}

fn append_polarity_warning(records: &[AtomSasaRecord], warnings: &mut Vec<crate::AnalysisWarning>) {
    let count = records
        .iter()
        .filter(|record| record.polarity == AtomPolarity::Unclassified)
        .count();
    if count > 0 {
        warnings.push(crate::AnalysisWarning::new(
            crate::WarningCode::UnsupportedPolarity,
            format!(
                "{count} atoms are outside the canonical Rosetta SasaFilter mapping; their area is reported as unclassified"
            ),
        ));
    }
}

pub(crate) fn calculate_atom_sasa_records(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    remove_hydrogens: bool,
    chains: &str,
) -> crate::ArpeggiaResult<Vec<AtomSasaRecord>> {
    calculate_atom_sasa_records_with_scheme(
        pdb,
        probe_radius,
        n_points,
        model_num,
        remove_hydrogens,
        chains,
        RadiusScheme::StandardProtor,
    )
}

pub(crate) fn calculate_sap_atom_sasa_records(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
) -> crate::ArpeggiaResult<Vec<AtomSasaRecord>> {
    calculate_atom_sasa_records_with_scheme(
        pdb,
        probe_radius,
        n_points,
        model_num,
        false,
        chains,
        RadiusScheme::SapReduce,
    )
}

#[allow(clippy::too_many_arguments)]
fn calculate_atom_sasa_records_with_scheme(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    remove_hydrogens: bool,
    chains: &str,
    radius_scheme: RadiusScheme,
) -> crate::ArpeggiaResult<Vec<AtomSasaRecord>> {
    let pdb_filtered = prepare_structure(pdb, model_num, remove_hydrogens, chains);
    calculate_prepared_atom_sasa_records(&pdb_filtered, probe_radius, n_points, None, radius_scheme)
}

fn calculate_prepared_atom_sasa_records(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    chains: Option<&HashSet<String>>,
    radius_scheme: RadiusScheme,
) -> crate::ArpeggiaResult<Vec<AtomSasaRecord>> {
    use rust_sasa::{Atom as SASAAtom, calculate_sasa_internal};

    let atom_hierarchy = pdb
        .atoms_with_hierarchy()
        .filter(|entity| chains.is_none_or(|chains| chains.contains(entity.chain().id())))
        .collect::<Vec<_>>();
    let atoms = atom_hierarchy
        .iter()
        .map(|x| {
            let radius = match radius_scheme {
                RadiusScheme::StandardProtor => rust_sasa::utils::get_protor_radius(
                    x.residue().name().unwrap_or(""),
                    x.atom().name(),
                )
                .or_else(|| {
                    x.atom()
                        .element()?
                        .atomic_radius()
                        .van_der_waals
                        .map(|radius| radius as f32)
                }),
                RadiusScheme::SapReduce => sap_reduce_radius(x),
            }
            .ok_or_else(|| {
                crate::ArpeggiaError::InvalidArgument(format!(
                    "atom {} ({}) has no usable van der Waals radius",
                    x.atom().serial_number(),
                    x.atom().name()
                ))
            })?;
            Ok(SASAAtom {
                position: [
                    x.atom().pos().0 as f32,
                    x.atom().pos().1 as f32,
                    x.atom().pos().2 as f32,
                ],
                radius,
                id: x.atom().serial_number(),
                parent_id: None,
            })
        })
        .collect::<crate::ArpeggiaResult<Vec<_>>>()?;
    let atom_sasa = calculate_sasa_internal(
        &atoms,
        probe_radius,
        n_points,
        rayon::current_num_threads() as isize,
    );

    atom_hierarchy
        .iter()
        .zip(atom_sasa)
        .map(|(hier, sasa)| {
            let entity = InteractingEntity::from_hier(hier);
            Ok(AtomSasaRecord {
                chain: entity.chain,
                resn: entity.resn,
                resi: i64::try_from(entity.resi).map_err(|_| {
                    crate::ArpeggiaError::InvalidArgument("residue identifier exceeds Int64".into())
                })?,
                insertion: entity.insertion,
                altloc: entity.altloc,
                atomn: entity.atomn,
                atomi: u64::try_from(entity.atomi).map_err(|_| {
                    crate::ArpeggiaError::InvalidArgument("atom identifier exceeds UInt64".into())
                })?,
                sasa,
                polarity: rosetta_sasa_filter_polarity(hier),
            })
        })
        .collect()
}

fn sap_reduce_radius(entity: &AtomConformerResidueChainModel) -> Option<f32> {
    let residue = entity.conformer().name();
    let atom = entity.atom();
    let atom_name = atom.name();
    Some(match atom.element()? {
        Element::C if atom_name == "C" => {
            if entity.residue().atoms().any(|atom| atom.name() == "OXT") {
                1.75
            } else {
                1.65
            }
        }
        Element::C => 1.75,
        Element::N => 1.55,
        Element::O => 1.4,
        Element::S => 1.8,
        Element::Se if residue == "MSE" => 1.8,
        Element::H if reduce_unit_hydrogen(residue, atom_name) => 1.0,
        Element::H => 1.17,
        element => element.atomic_radius().van_der_waals? as f32,
    })
}

fn reduce_unit_hydrogen(residue: &str, atom: &str) -> bool {
    if matches!(atom, "H" | "HN" | "H1" | "H2" | "H3" | "1H" | "2H" | "3H") {
        return true;
    }
    match residue {
        "ARG" => atom == "HE" || atom.contains("HH"),
        "ASN" => atom.contains("HD2"),
        "CYS" => atom == "HG",
        "GLN" => atom.contains("HE2"),
        "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP" => {
            matches!(atom, "HD1" | "HE2")
        }
        "LYS" => atom.contains("HZ"),
        "PHE" | "TRP" | "TYR" => atom != "HA" && !atom.contains("HB"),
        "SER" => atom == "HG",
        "THR" => atom == "HG1",
        _ => false,
    }
}

fn rosetta_sasa_filter_polarity(entity: &AtomConformerResidueChainModel) -> AtomPolarity {
    // [ACCEPTED DEVIATION] This reproduces Rosetta SasaFilter's atom partition
    // while areas remain Shrake-Rupley/ProtOr rather than LeGrand atom-type SASA.
    // [WARNING] Cross-tool benchmarking shows these Shrake-Rupley areas are not
    // numerically equivalent to Rosetta's full-atom LeGrand areas.
    let residue = entity.conformer().name();
    let atom = entity.atom().name();
    if !matches!(
        residue,
        "ALA"
            | "ARG"
            | "ASN"
            | "ASP"
            | "CYS"
            | "GLN"
            | "GLU"
            | "GLY"
            | "HIS"
            | "HID"
            | "HIE"
            | "HIP"
            | "HSD"
            | "HSE"
            | "HSP"
            | "ILE"
            | "LEU"
            | "LYS"
            | "MET"
            | "PHE"
            | "PRO"
            | "SER"
            | "THR"
            | "TRP"
            | "TYR"
            | "VAL"
    ) {
        return AtomPolarity::Unclassified;
    }

    let polar = if entity.atom().element() == Some(&Element::H) {
        matches!(atom, "H" | "H1" | "H2" | "H3" | "HN" | "1H" | "2H" | "3H")
            || matches!(
                (residue, atom),
                ("ARG", "HE" | "HH11" | "HH12" | "HH21" | "HH22")
                    | ("ASN", "HD21" | "HD22")
                    | ("GLN", "HE21" | "HE22")
                    | (
                        "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP",
                        "HD1" | "HE2"
                    )
                    | ("LYS", "HZ1" | "HZ2" | "HZ3")
                    | ("SER" | "CYS", "HG")
                    | ("THR", "HG1")
                    | ("TRP", "HE1")
                    | ("TYR", "HH")
            )
    } else {
        matches!(atom, "O" | "OXT")
            || (atom == "N" && residue != "PRO")
            || matches!(
                (residue, atom),
                ("ARG", "NE" | "NH1" | "NH2")
                    | ("ASN", "OD1" | "ND2")
                    | ("ASP", "OD1" | "OD2")
                    | ("GLN", "OE1" | "NE2")
                    | ("GLU", "OE1" | "OE2")
                    | (
                        "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP",
                        "ND1" | "NE2"
                    )
                    | ("LYS", "NZ")
                    | ("SER", "OG")
                    | ("THR", "OG1")
                    | ("TRP", "NE1")
                    | ("TYR", "OH")
            )
    };
    if polar {
        AtomPolarity::Polar
    } else {
        AtomPolarity::Hydrophobic
    }
}

pub(crate) fn atom_sasa_records_to_dataframe(records: &[AtomSasaRecord]) -> DataFrame {
    df!(
        "atomi" => records.iter().map(|x| x.atomi).collect::<Vec<u64>>(),
        "sasa" => records.iter().map(|x| x.sasa).collect::<Vec<f32>>(),
        "polarity" => records.iter().map(|x| x.polarity.as_str()).collect::<Vec<&str>>(),
        "chain" => records.iter().map(|x| x.chain.clone()).collect::<Vec<String>>(),
        "resn" => records.iter().map(|x| x.resn.clone()).collect::<Vec<String>>(),
        "resi" => records.iter().map(|x| x.resi).collect::<Vec<i64>>(),
        "insertion" => records.iter().map(|x| x.insertion.clone()).collect::<Vec<String>>(),
        "altloc" => records.iter().map(|x| x.altloc.clone()).collect::<Vec<String>>(),
        "atomn" => records.iter().map(|x| x.atomn.clone()).collect::<Vec<String>>(),
    )
    .unwrap()
    .sort(["atomi"], Default::default())
    .unwrap()
}

/// Calculate solvent accessible surface area (SASA) aggregated by residue.
///
/// Aggregates the same selected atom records used by atom-level SASA.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `chains` - Comma-separated chain IDs to include (e.g., "A,B,C"). Empty string includes all chains.
///
/// # Returns
///
/// A Polars `DataFrame` with columns:
/// - `chain`, `resn`, `resi`, `insertion`, `sasa`, `polar_sasa`,
///   `hydrophobic_sasa`, `unclassified_sasa`
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_residue_sasa};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let pdb = load_model(&input_file).unwrap().value;
///
/// // Calculate SASA for all chains
/// let sasa_df = get_residue_sasa(&pdb, 1.4, 100, 0, "").unwrap().value;
/// println!("Calculated SASA for {} residues", sasa_df.height());
///
/// // Calculate SASA for only chain A
/// let sasa_a = get_residue_sasa(&pdb, 1.4, 100, 0, "A").unwrap().value;
/// ```
pub fn get_residue_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
) -> crate::ArpeggiaResult<crate::Analysis<DataFrame>> {
    let analysis = validate_sasa_input(pdb, probe_radius, n_points, model_num, chains)?;
    let atom_records = calculate_atom_sasa_records(
        &analysis.value,
        probe_radius,
        n_points,
        model_num,
        true,
        chains,
    )?;
    let mut warnings = analysis.warnings;
    append_polarity_warning(&atom_records, &mut warnings);
    Ok(crate::Analysis::new(
        residue_sasa_records_to_dataframe(&aggregate_residue_records(atom_records)),
        warnings,
    ))
}

fn aggregate_residue_records(records: Vec<AtomSasaRecord>) -> Vec<ResidueSasaRecord> {
    let mut totals = BTreeMap::<(String, String, i64, String), [f32; 4]>::new();
    for atom in records {
        let values = totals
            .entry((atom.chain, atom.resn, atom.resi, atom.insertion))
            .or_default();
        values[0] += atom.sasa;
        values[match atom.polarity {
            AtomPolarity::Polar => 1,
            AtomPolarity::Hydrophobic => 2,
            AtomPolarity::Unclassified => 3,
        }] += atom.sasa;
    }
    totals
        .into_iter()
        .map(
            |((chain, resn, resi, insertion), values)| ResidueSasaRecord {
                chain,
                resn,
                resi,
                insertion,
                sasa: values[0],
                polar_sasa: values[1],
                hydrophobic_sasa: values[2],
                unclassified_sasa: values[3],
            },
        )
        .collect()
}

pub(crate) fn residue_sasa_records_to_dataframe(records: &[ResidueSasaRecord]) -> DataFrame {
    df!(
        "chain" => records.iter().map(|r| r.chain.clone()).collect::<Vec<String>>(),
        "resn" => records.iter().map(|r| r.resn.clone()).collect::<Vec<String>>(),
        "resi" => records.iter().map(|r| r.resi).collect::<Vec<i64>>(),
        "insertion" => records.iter().map(|r| r.insertion.clone()).collect::<Vec<String>>(),
        "sasa" => records.iter().map(|r| r.sasa).collect::<Vec<f32>>(),
        "polar_sasa" => records.iter().map(|r| r.polar_sasa).collect::<Vec<f32>>(),
        "hydrophobic_sasa" => records.iter().map(|r| r.hydrophobic_sasa).collect::<Vec<f32>>(),
        "unclassified_sasa" => records.iter().map(|r| r.unclassified_sasa).collect::<Vec<f32>>(),
    )
    .unwrap()
    .sort(["chain", "resi", "insertion"], Default::default())
    .unwrap()
}

/// Calculate solvent accessible surface area (SASA) aggregated by chain.
///
/// Aggregates the same selected atom records used by atom-level SASA.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `chains` - Comma-separated chain IDs to include (e.g., "A,B,C"). Empty string includes all chains.
///
/// # Returns
///
/// A Polars `DataFrame` with columns:
/// - `chain`, `sasa`, `polar_sasa`, `hydrophobic_sasa`, `unclassified_sasa`
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_chain_sasa};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let pdb = load_model(&input_file).unwrap().value;
///
/// // Calculate SASA for all chains
/// let sasa_df = get_chain_sasa(&pdb, 1.4, 100, 0, "").unwrap().value;
/// println!("Calculated SASA for {} chains", sasa_df.height());
///
/// // Calculate SASA for only chains A and B
/// let sasa_ab = get_chain_sasa(&pdb, 1.4, 100, 0, "A,B").unwrap().value;
/// ```
pub fn get_chain_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
) -> crate::ArpeggiaResult<crate::Analysis<DataFrame>> {
    let analysis = validate_sasa_input(pdb, probe_radius, n_points, model_num, chains)?;
    let atom_records = calculate_atom_sasa_records(
        &analysis.value,
        probe_radius,
        n_points,
        model_num,
        true,
        chains,
    )?;
    let mut warnings = analysis.warnings;
    append_polarity_warning(&atom_records, &mut warnings);
    Ok(crate::Analysis::new(
        chain_sasa_records_to_dataframe(&aggregate_chain_records(atom_records)),
        warnings,
    ))
}

fn aggregate_chain_records(records: Vec<AtomSasaRecord>) -> Vec<ChainSasaRecord> {
    let mut totals = BTreeMap::<String, [f32; 4]>::new();
    for atom in records {
        let values = totals.entry(atom.chain).or_default();
        values[0] += atom.sasa;
        values[match atom.polarity {
            AtomPolarity::Polar => 1,
            AtomPolarity::Hydrophobic => 2,
            AtomPolarity::Unclassified => 3,
        }] += atom.sasa;
    }
    totals
        .into_iter()
        .map(|(chain, values)| ChainSasaRecord {
            chain,
            sasa: values[0],
            polar_sasa: values[1],
            hydrophobic_sasa: values[2],
            unclassified_sasa: values[3],
        })
        .collect()
}

pub(crate) fn chain_sasa_records_to_dataframe(records: &[ChainSasaRecord]) -> DataFrame {
    df!(
        "chain" => records.iter().map(|r| r.chain.clone()).collect::<Vec<String>>(),
        "sasa" => records.iter().map(|r| r.sasa).collect::<Vec<f32>>(),
        "polar_sasa" => records.iter().map(|r| r.polar_sasa).collect::<Vec<f32>>(),
        "hydrophobic_sasa" => records.iter().map(|r| r.hydrophobic_sasa).collect::<Vec<f32>>(),
        "unclassified_sasa" => records.iter().map(|r| r.unclassified_sasa).collect::<Vec<f32>>(),
    )
    .unwrap()
    .sort(["chain"], Default::default())
    .unwrap()
}

/// Calculate the buried surface area (dSASA) at the interface between two chain groups.
///
/// The buried surface area is calculated as:
/// dSASA = SASA_group1 + SASA_group2 - SASA_complex
///
/// This is the two-sided buried area. Divide by two only when a one-sided
/// interface-area convention is explicitly desired.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `groups` - Chain groups specification (e.g., "A,B/C,D")
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
///
/// # Returns
///
/// The buried surface area at the interface in square Ångströms.
pub fn get_dsasa(
    pdb: &PDB,
    groups: &str,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
) -> crate::ArpeggiaResult<crate::Analysis<f32>> {
    Ok(
        get_dsasa_components(pdb, groups, probe_radius, n_points, model_num)?
            .map(|result| result.dsasa),
    )
}

/// Calculate two-sided dSASA and its additive polarity partition.
pub fn get_dsasa_components(
    pdb: &PDB,
    groups: &str,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
) -> crate::ArpeggiaResult<crate::Analysis<DsasaResult>> {
    let analysis = validate_sasa_input(pdb, probe_radius, n_points, model_num, "")?;
    let pdb = &analysis.value;
    // Get all chains in the PDB
    let all_chains: HashSet<String> = pdb.chains().map(|c| c.id().to_string()).collect();

    // Parse groups
    let (group1_chains, group2_chains) = parse_groups(&all_chains, groups)?;
    if !group1_chains.is_disjoint(&group2_chains) {
        return Err(crate::ArpeggiaError::InvalidArgument(
            "dSASA chain groups must be disjoint".into(),
        ));
    }

    // Get combined chains (union of both groups)
    let combined_group_chains: HashSet<String> =
        group1_chains.union(&group2_chains).cloned().collect();

    let prepared = StructurePreparation::new(model_num)
        .remove_hydrogens(true)
        .remove_solvent_and_ions(true)
        .chain_set(&combined_group_chains)
        .prepare(pdb);

    let combined_atom_sasa = calculate_prepared_atom_sasa_records(
        &prepared,
        probe_radius,
        n_points,
        None,
        RadiusScheme::StandardProtor,
    )?;
    let mut warnings = analysis.warnings;
    append_polarity_warning(&combined_atom_sasa, &mut warnings);
    let combined_sasa = aggregate_chain_records(combined_atom_sasa);
    let combined_total = sum_chain_sasa(&combined_sasa);

    let group1_sasa = aggregate_chain_records(calculate_prepared_atom_sasa_records(
        &prepared,
        probe_radius,
        n_points,
        Some(&group1_chains),
        RadiusScheme::StandardProtor,
    )?);
    let group1_total = sum_chain_sasa(&group1_sasa);

    let group2_sasa = aggregate_chain_records(calculate_prepared_atom_sasa_records(
        &prepared,
        probe_radius,
        n_points,
        Some(&group2_chains),
        RadiusScheme::StandardProtor,
    )?);
    let group2_total = sum_chain_sasa(&group2_sasa);

    // [ACCEPTED DEVIATION] Arpeggia reports the reference-compatible two-sided
    // buried area, without the optional interface-area division by two.
    Ok(crate::Analysis::new(
        DsasaResult {
            dsasa: group1_total.dsasa + group2_total.dsasa - combined_total.dsasa,
            polar_dsasa: group1_total.polar_dsasa + group2_total.polar_dsasa
                - combined_total.polar_dsasa,
            hydrophobic_dsasa: group1_total.hydrophobic_dsasa + group2_total.hydrophobic_dsasa
                - combined_total.hydrophobic_dsasa,
            unclassified_dsasa: group1_total.unclassified_dsasa + group2_total.unclassified_dsasa
                - combined_total.unclassified_dsasa,
        },
        warnings,
    ))
}

fn sum_chain_sasa(records: &[ChainSasaRecord]) -> DsasaResult {
    records
        .iter()
        .fold(DsasaResult::default(), |mut total, record| {
            total.dsasa += record.sasa;
            total.polar_dsasa += record.polar_sasa;
            total.hydrophobic_dsasa += record.hydrophobic_sasa;
            total.unclassified_dsasa += record.unclassified_sasa;
            total
        })
}

/// Maximum solvent accessible surface area (`MaxASA`) values for amino acids.
///
/// Values are from Tien et al. (2013) "Maximum Allowed Solvent Accessibilities of Residues in Proteins"
/// PLOS ONE. These theoretical values represent the maximum possible SASA for each amino acid
/// in a Gly-X-Gly tripeptide.
///
/// Returns the `MaxASA` value in Å² for a given 3-letter amino acid code, or None if unknown.
pub fn get_max_asa(resn: &str) -> Option<f32> {
    match resn.to_uppercase().as_str() {
        "ALA" => Some(129.0),
        "ARG" => Some(274.0),
        "ASN" => Some(195.0),
        "ASP" => Some(193.0),
        "CYS" => Some(167.0),
        "GLU" => Some(223.0),
        "GLN" => Some(225.0),
        "GLY" => Some(104.0),
        "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP" | "MET" | "MSE" => Some(224.0),
        "ILE" => Some(197.0),
        "LEU" => Some(201.0),
        "LYS" => Some(236.0),
        "PHE" => Some(240.0),
        "PRO" => Some(159.0),
        "SER" => Some(155.0),
        "THR" => Some(172.0),
        "TRP" => Some(285.0),
        "TYR" => Some(263.0),
        "VAL" => Some(174.0),
        _ => None,
    }
}

/// Calculate relative solvent accessible surface area (RSA) for each residue.
///
/// RSA is calculated as the ratio of observed SASA to the maximum possible SASA
/// for each amino acid type, based on Tien et al. (2013) theoretical values.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `chains` - Comma-separated chain IDs to include (e.g., "A,B,C"). Empty string includes all chains.
///
/// # Returns
///
/// A Polars `DataFrame` with columns:
/// - chain, resn, resi, insertion, sasa, `polar_sasa`, `hydrophobic_sasa`,
///   `unclassified_sasa`, `relative_sasa`
///
/// The `relative_sasa` column contains values between 0 and ~1 (can slightly exceed 1 due to
/// structural context), or null for non-standard amino acids.
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_relative_sasa};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let pdb = load_model(&input_file).unwrap().value;
///
/// // Calculate RSA for all chains
/// let rsa_df = get_relative_sasa(&pdb, 1.4, 100, 0, "").unwrap().value;
///
/// // Calculate RSA for only chain A
/// let rsa_a = get_relative_sasa(&pdb, 1.4, 100, 0, "A").unwrap().value;
/// ```
pub fn get_relative_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
) -> crate::ArpeggiaResult<crate::Analysis<DataFrame>> {
    let analysis = validate_sasa_input(pdb, probe_radius, n_points, model_num, chains)?;
    let atom_records = calculate_atom_sasa_records(
        &analysis.value,
        probe_radius,
        n_points,
        model_num,
        true,
        chains,
    )?;
    let mut warnings = analysis.warnings;
    append_polarity_warning(&atom_records, &mut warnings);
    let residue_records = aggregate_residue_records(atom_records);
    Ok(crate::Analysis::new(
        relative_sasa_records_to_dataframe(&relative_sasa_from_residue_records(residue_records)),
        warnings,
    ))
}

fn relative_sasa_from_residue_records(records: Vec<ResidueSasaRecord>) -> Vec<RelativeSasaRecord> {
    records
        .into_iter()
        .map(|record| {
            let relative_sasa = get_max_asa(&record.resn)
                .filter(|max_sasa| *max_sasa > 0.0)
                .map(|max_sasa| record.sasa / max_sasa);

            RelativeSasaRecord {
                chain: record.chain,
                resn: record.resn,
                resi: record.resi,
                insertion: record.insertion,
                sasa: record.sasa,
                polar_sasa: record.polar_sasa,
                hydrophobic_sasa: record.hydrophobic_sasa,
                unclassified_sasa: record.unclassified_sasa,
                relative_sasa,
            }
        })
        .collect()
}

pub(crate) fn validate_sasa_input(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
) -> crate::ArpeggiaResult<crate::Analysis<PDB>> {
    if !probe_radius.is_finite() || probe_radius <= 0.0 {
        return Err(crate::ArpeggiaError::InvalidArgument(
            "probe_radius must be finite and positive".into(),
        ));
    }
    if n_points == 0 {
        return Err(crate::ArpeggiaError::InvalidArgument(
            "n_points must be positive".into(),
        ));
    }
    let models = pdb
        .models()
        .map(|model| model.serial_number())
        .collect::<Vec<_>>();
    if models.is_empty() || (model_num != 0 && !models.contains(&model_num)) {
        return Err(crate::ArpeggiaError::InvalidArgument(format!(
            "model {model_num} does not exist"
        )));
    }
    let available = pdb
        .chains()
        .map(|chain| chain.id().to_string())
        .collect::<HashSet<_>>();
    let requested = chains
        .split(',')
        .map(str::trim)
        .filter(|chain| !chain.is_empty())
        .collect::<Vec<_>>();
    if let Some(chain) = requested.iter().find(|chain| !available.contains(**chain)) {
        return Err(crate::ArpeggiaError::InvalidArgument(format!(
            "chain {chain} does not exist"
        )));
    }
    let mut prepared = pdb.clone();
    let warnings = crate::structure::select_conformers(&mut prepared);
    Ok(crate::Analysis::new(prepared, warnings))
}

fn relative_sasa_records_to_dataframe(records: &[RelativeSasaRecord]) -> DataFrame {
    df!(
        "chain" => records.iter().map(|r| r.chain.clone()).collect::<Vec<String>>(),
        "resn" => records.iter().map(|r| r.resn.clone()).collect::<Vec<String>>(),
        "resi" => records.iter().map(|r| r.resi).collect::<Vec<i64>>(),
        "insertion" => records.iter().map(|r| r.insertion.clone()).collect::<Vec<String>>(),
        "sasa" => records.iter().map(|r| r.sasa).collect::<Vec<f32>>(),
        "polar_sasa" => records.iter().map(|r| r.polar_sasa).collect::<Vec<f32>>(),
        "hydrophobic_sasa" => records.iter().map(|r| r.hydrophobic_sasa).collect::<Vec<f32>>(),
        "unclassified_sasa" => records.iter().map(|r| r.unclassified_sasa).collect::<Vec<f32>>(),
        "relative_sasa" => records.iter().map(|r| r.relative_sasa).collect::<Vec<Option<f32>>>(),
    )
    .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{WarningCode, load_model, run_with_threads};

    fn sum_float_col(df: &DataFrame, colname: &str) -> f32 {
        df.column(colname)
            .unwrap()
            .f32()
            .unwrap()
            .sum()
            .unwrap_or(0.0)
    }

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
    fn test_get_atom_sasa_returns_data() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_atom_sasa(&pdb, 1.4, 100, 0, true, ""));
        let df = df.unwrap().value;

        // Check that we get results
        assert!(df.height() > 0, "SASA DataFrame should not be empty");

        // Check that the expected columns exist
        let columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(
            columns.contains(&"atomi".to_string()),
            "Should have 'atomi' column"
        );
        assert!(
            columns.contains(&"sasa".to_string()),
            "Should have 'sasa' column"
        );
        assert!(
            columns.contains(&"chain".to_string()),
            "Should have 'chain' column"
        );
        assert!(
            columns.contains(&"resn".to_string()),
            "Should have 'resn' column"
        );
        assert!(
            columns.contains(&"resi".to_string()),
            "Should have 'resi' column"
        );
        assert!(
            columns.contains(&"atomn".to_string()),
            "Should have 'atomn' column"
        );
        assert!(columns.contains(&"polarity".to_string()));
        assert_eq!(df.column("resi").unwrap().dtype(), &DataType::Int64);
        assert_eq!(df.column("atomi").unwrap().dtype(), &DataType::UInt64);
    }

    #[test]
    fn test_get_atom_sasa_values_reasonable() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_atom_sasa(&pdb, 1.4, 100, 0, true, ""));
        let df = df.unwrap().value;

        // Get the SASA column and check values are non-negative
        let sasa_col = df.column("sasa").unwrap();
        let sasa_values: Vec<f32> = sasa_col.f32().unwrap().iter().flatten().collect();

        assert!(
            sasa_values.iter().all(|&v| v >= 0.0),
            "All SASA values should be non-negative"
        );

        // Check that some atoms have non-zero SASA (surface exposed)
        let non_zero_count = sasa_values.iter().filter(|&&v| v > 0.0).count();
        assert!(non_zero_count > 0, "Some atoms should have non-zero SASA");
    }

    #[test]
    fn unsupported_polarity_is_reported_without_dropping_area() {
        let input = b"HETATM    1 SE   MSE A   1       0.000   0.000   0.000  1.00 20.00          SE  \nEND\n";
        let pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(std::io::BufReader::new(input.as_slice()))
            .unwrap()
            .0;
        let analysis = get_atom_sasa(&pdb, 1.4, 100, 0, true, "").unwrap();
        assert_eq!(analysis.value.height(), 1);
        assert_eq!(
            analysis
                .value
                .column("polarity")
                .unwrap()
                .str()
                .unwrap()
                .get(0),
            Some("unclassified")
        );
        assert!(
            analysis
                .warnings
                .iter()
                .any(|warning| warning.code == WarningCode::UnsupportedPolarity)
        );
    }

    #[test]
    fn missing_radius_returns_a_typed_error() {
        let input = b"ATOM      1  QQ  ALA A   1       0.000   0.000   0.000  1.00 20.00              \nEND\n";
        let pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(std::io::BufReader::new(input.as_slice()))
            .unwrap()
            .0;
        assert!(matches!(
            get_atom_sasa(&pdb, 1.4, 100, 0, true, ""),
            Err(crate::ArpeggiaError::InvalidArgument(message))
                if message.contains("van der Waals radius")
        ));
    }

    #[test]
    fn test_get_residue_sasa_returns_data() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_residue_sasa(&pdb, 1.4, 100, 0, ""));
        let df = df.unwrap().value;

        // Check that we get results
        assert!(
            df.height() > 0,
            "Residue SASA DataFrame should not be empty"
        );

        // Check that the expected columns exist
        let columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(
            columns.contains(&"chain".to_string()),
            "Should have 'chain' column"
        );
        assert!(
            columns.contains(&"resn".to_string()),
            "Should have 'resn' column"
        );
        assert!(
            columns.contains(&"resi".to_string()),
            "Should have 'resi' column"
        );
        assert!(
            columns.contains(&"insertion".to_string()),
            "Should have 'insertion' column"
        );
        assert!(
            columns.contains(&"sasa".to_string()),
            "Should have 'sasa' column"
        );
        assert!(
            columns.contains(&"polar_sasa".to_string()),
            "Should have 'polar_sasa' column"
        );
        assert!(columns.contains(&"hydrophobic_sasa".to_string()));
        assert!(columns.contains(&"unclassified_sasa".to_string()));
    }

    #[test]
    fn test_get_residue_sasa_aggregation() {
        let pdb = load_ubiquitin();
        let atom_df = run_with_threads(1, || get_atom_sasa(&pdb, 1.4, 100, 0, true, ""));
        let atom_df = atom_df.unwrap().value;
        let residue_df = run_with_threads(1, || get_residue_sasa(&pdb, 1.4, 100, 0, ""));
        let residue_df = residue_df.unwrap().value;

        // There should be fewer rows in residue-level than atom-level
        assert!(
            residue_df.height() < atom_df.height(),
            "Residue-level should have fewer rows than atom-level: {} vs {}",
            residue_df.height(),
            atom_df.height()
        );

        // All levels aggregate the exact same per-atom calculation.
        let atom_total: f32 = sum_float_col(&atom_df, "sasa");
        let residue_total: f32 = sum_float_col(&residue_df, "sasa");
        assert!(
            (residue_total - atom_total).abs() < 0.001,
            "atom={atom_total}, residue={residue_total}"
        );
        let partition_total = sum_float_col(&residue_df, "polar_sasa")
            + sum_float_col(&residue_df, "hydrophobic_sasa")
            + sum_float_col(&residue_df, "unclassified_sasa");
        assert!((partition_total - residue_total).abs() < 0.001);

        let chain_df = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));
        let chain_df = chain_df.unwrap().value;
        assert!((sum_float_col(&chain_df, "sasa") - atom_total).abs() < 0.001);
    }

    #[test]
    fn test_get_chain_sasa_returns_data() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));
        let df = df.unwrap().value;

        // Check that we get results
        assert!(df.height() > 0, "Chain SASA DataFrame should not be empty");

        // Check that the expected columns exist
        let columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(
            columns.contains(&"chain".to_string()),
            "Should have 'chain' column"
        );
        assert!(
            columns.contains(&"sasa".to_string()),
            "Should have 'sasa' column"
        );
    }

    #[test]
    fn test_get_chain_sasa_single_chain() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));
        let df = df.unwrap().value;

        // Ubiquitin (1ubq) has a single chain A
        assert_eq!(df.height(), 1, "1ubq should have 1 chain");

        // Check that the chain is A
        let chain_col = df.column("chain").unwrap();
        let chain_id = chain_col.str().unwrap().get(0).unwrap();
        assert_eq!(chain_id, "A", "Chain should be A");
    }

    #[test]
    fn test_get_chain_sasa_multi_chain() {
        let pdb = load_multi_chain();
        let df = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));
        let df = df.unwrap().value;

        // 6bft should have multiple chains
        assert!(df.height() > 1, "6bft should have multiple chains");

        // Check that all SASA values are non-negative
        let sasa_col = df.column("sasa").unwrap();
        let sasa_values: Vec<f32> = sasa_col.f32().unwrap().iter().flatten().collect();

        assert!(
            sasa_values.iter().all(|&v| v >= 0.0),
            "All chain SASA values should be non-negative"
        );
    }

    #[test]
    fn test_sasa_probe_radius_effect() {
        let pdb = load_ubiquitin();

        // Smaller probe radius should result in larger SASA
        let small_probe = run_with_threads(1, || get_chain_sasa(&pdb, 1.0, 100, 0, ""));
        let small_probe = small_probe.unwrap().value;
        let large_probe = run_with_threads(1, || get_chain_sasa(&pdb, 2.0, 100, 0, ""));
        let large_probe = large_probe.unwrap().value;

        let small_sasa: f32 = small_probe
            .column("sasa")
            .unwrap()
            .f32()
            .unwrap()
            .get(0)
            .unwrap();
        let large_sasa: f32 = large_probe
            .column("sasa")
            .unwrap()
            .f32()
            .unwrap()
            .get(0)
            .unwrap();

        assert!(
            small_sasa > large_sasa,
            "Smaller probe radius should give larger SASA: {small_sasa} vs {large_sasa}"
        );
    }

    #[test]
    fn test_sasa_regression_ubiquitin() {
        // Regression test to ensure SASA values remain consistent
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));
        let df = df.unwrap().value;

        let total_sasa: f32 = df.column("sasa").unwrap().f32().unwrap().get(0).unwrap();

        // Expected value from rust-sasa 0.9.0 with default parameters
        // This should be around 4813 Å² for ubiquitin
        let expected_sasa = 4813.0;
        let tolerance = 100.0; // Allow some tolerance for minor differences

        assert!(
            (total_sasa - expected_sasa).abs() < tolerance,
            "Ubiquitin total SASA should be around {expected_sasa} Å², got {total_sasa} Å²"
        );
    }

    #[test]
    fn test_get_dsasa_returns_positive() {
        // 6bft has multiple chains with interfaces
        let pdb = load_multi_chain();

        // Calculate dSASA between groups A,B,C and G,H,L
        let dsasa = run_with_threads(1, || get_dsasa(&pdb, "A,B,C/G,H,L", 1.4, 100, 0));
        let dsasa = dsasa.unwrap().value;

        // dSASA should be positive for an interface
        assert!(dsasa > 0.0, "dSASA should be positive, got {dsasa}");
    }

    #[test]
    fn test_get_dsasa_interface_value() {
        // 6bft has multiple chains with interfaces
        let pdb = load_multi_chain();

        // Calculate dSASA between groups A,B,C and G,H,L
        let dsasa = run_with_threads(1, || get_dsasa(&pdb, "C/H,L", 1.4, 100, 0));
        let dsasa = dsasa.unwrap().value;

        // Regression test: the dSASA should be around 1644-1665 Å²
        // as calculated from PyMOL and Rosetta InterfaceAnalyzer
        let expected_dsasa = 1650.0;
        let tolerance = 50.0; // Allow some tolerance

        assert!(
            (dsasa - expected_dsasa).abs() < tolerance,
            "6bft dSASA should be around {expected_dsasa} Å², got {dsasa} Å²"
        );
    }

    #[test]
    fn test_get_dsasa_symmetric() {
        // dSASA should be the same regardless of which group is first
        let pdb = load_multi_chain();

        let dsasa1 = run_with_threads(1, || get_dsasa(&pdb, "A,B,C/G,H,L", 1.4, 100, 0));
        let dsasa2 = run_with_threads(1, || get_dsasa(&pdb, "G,H,L/A,B,C", 1.4, 100, 0));
        let dsasa1 = dsasa1.unwrap().value;
        let dsasa2 = dsasa2.unwrap().value;

        let diff = (dsasa1 - dsasa2).abs();
        assert!(
            diff < 1.0,
            "dSASA should be symmetric: {dsasa1} vs {dsasa2}"
        );
    }

    #[test]
    fn dsasa_is_two_sided_partitioned_and_requires_disjoint_groups() {
        let pdb = load_multi_chain();
        let result = run_with_threads(1, || {
            get_dsasa_components(&pdb, "C/H,L", 1.4, 100, 0)
                .unwrap()
                .value
        });
        assert!(result.dsasa > 0.0);
        assert!(
            (result.polar_dsasa + result.hydrophobic_dsasa + result.unclassified_dsasa
                - result.dsasa)
                .abs()
                < 0.01
        );
        assert!(get_dsasa_components(&pdb, "/", 1.4, 100, 0).is_err());
        assert!(get_dsasa_components(&pdb, "C/H,L", f32::NAN, 100, 0).is_err());
        assert!(get_dsasa_components(&pdb, "C/H,L", 1.4, 0, 0).is_err());
    }

    #[test]
    fn test_get_relative_sasa_returns_data() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_relative_sasa(&pdb, 1.4, 100, 0, ""));
        let df = df.unwrap().value;

        // Check that we get results
        assert!(
            df.height() > 0,
            "Relative SASA DataFrame should not be empty"
        );

        // Check that the expected columns exist
        let columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(
            columns.contains(&"chain".to_string()),
            "Should have 'chain' column"
        );
        assert!(
            columns.contains(&"resn".to_string()),
            "Should have 'resn' column"
        );
        assert!(
            columns.contains(&"resi".to_string()),
            "Should have 'resi' column"
        );
        assert!(
            columns.contains(&"sasa".to_string()),
            "Should have 'sasa' column"
        );
        assert!(
            columns.contains(&"relative_sasa".to_string()),
            "Should have 'relative_sasa' column"
        );
    }

    #[test]
    fn test_get_relative_sasa_values_bounded() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_relative_sasa(&pdb, 1.4, 100, 0, ""));
        let df = df.unwrap().value;

        // Get relative_sasa values
        let rsa_values: Vec<f32> = df
            .column("relative_sasa")
            .unwrap()
            .f32()
            .unwrap()
            .iter()
            .flatten()
            .collect();

        // All values should be non-negative
        assert!(
            rsa_values.iter().all(|&v| v >= 0.0),
            "All relative_sasa values should be non-negative"
        );

        // Most values should be <= 1.0 (some may slightly exceed due to structural context)
        let below_threshold = rsa_values.iter().filter(|&&v| v <= 1.5).count();
        let ratio = below_threshold as f64 / rsa_values.len() as f64;
        assert!(
            ratio > 0.95,
            "Most relative_sasa values should be <= 1.5: ratio={ratio}"
        );
    }

    #[test]
    fn test_get_max_asa_standard_amino_acids() {
        // Test that all standard amino acids have MaxASA values
        let amino_acids = [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        ];

        for aa in &amino_acids {
            let max_asa = get_max_asa(aa);
            assert!(max_asa.is_some(), "Should have MaxASA value for {aa}");
            assert!(max_asa.unwrap() > 0.0, "MaxASA for {aa} should be positive");
        }
    }

    #[test]
    fn test_get_max_asa_unknown_residue() {
        // Unknown residues should return None
        assert!(get_max_asa("XXX").is_none());
        assert_eq!(get_max_asa("HID"), get_max_asa("HIS"));
        assert_eq!(get_max_asa("MSE"), get_max_asa("MET"));
        assert!(get_max_asa("HOH").is_none());
        assert!(get_max_asa("").is_none());
    }

    #[test]
    fn test_chain_filter_empty_keeps_all() {
        let pdb = load_multi_chain();

        // Empty chain filter should keep all chains
        let df_all = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));
        let df_all = df_all.unwrap().value;
        let chain_count_all = df_all.height();

        assert!(
            chain_count_all > 1,
            "Multi-chain structure should have multiple chains: {chain_count_all}"
        );
    }

    #[test]
    fn test_chain_filter_single_chain() {
        let pdb = load_multi_chain();

        // Filter to only chain A
        let df_a = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, "A"));
        let df_a = df_a.unwrap().value;

        assert_eq!(
            df_a.height(),
            1,
            "Should have only one chain when filtering to A"
        );

        // Check that the chain is A
        let chain_col = df_a.column("chain").unwrap();
        let chain_id = chain_col.str().unwrap().get(0).unwrap();
        assert_eq!(chain_id, "A", "Chain should be A");
    }

    #[test]
    fn test_chain_filter_multiple_chains() {
        let pdb = load_multi_chain();

        // Filter to chains A, B
        let df_ab = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, "A,B"));
        let df_ab = df_ab.unwrap().value;
        let df_all = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));
        let df_all = df_all.unwrap().value;

        assert!(
            df_ab.height() <= df_all.height(),
            "Filtered results should have equal or fewer chains: {} vs {}",
            df_ab.height(),
            df_all.height()
        );

        // Check that we only have A and B chains
        let chain_col = df_ab.column("chain").unwrap();
        let chain_ids: Vec<&str> = chain_col.str().unwrap().iter().flatten().collect();
        for chain_id in &chain_ids {
            assert!(
                *chain_id == "A" || *chain_id == "B",
                "Only A and B chains should be present, got: {chain_id}"
            );
        }
    }

    #[test]
    fn reduce_histidine_hydrogen_radii_match_rosetta_types() {
        assert!(reduce_unit_hydrogen("HIS", "HE2"));
        assert!(!reduce_unit_hydrogen("HIS", "HD2"));
    }
}
