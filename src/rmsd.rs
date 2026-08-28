//! Rigid protein-structure superposition and atom correspondence.

use crate::contacts::residues::ResidueExt;
use crate::structure::selected_model;
use crate::{Analysis, AnalysisWarning, ArpeggiaError, ArpeggiaResult, WarningCode};
use nalgebra as na;
use pdbtbx::{ContainsAtomConformer, Element, PDB};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};

/// Atom population used for structural superposition.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, clap::ValueEnum)]
pub enum AtomSubset {
    /// Alpha-carbon atoms only.
    #[default]
    Ca,
    /// Protein N, CA, C, O, and OXT atoms.
    Backbone,
    /// All non-hydrogen atoms in protein residues and terminal caps.
    Heavy,
    /// All atoms in protein residues and terminal caps.
    All,
}

/// Calculate optimally superposed RMSD between equal-sized coordinate arrays.
///
/// Coordinates are uniformly weighted. Reflection and scaling are forbidden.
pub fn kabsch_rmsd(reference: &[[f64; 3]], mobile: &[[f64; 3]]) -> ArpeggiaResult<f64> {
    validate_coordinate_arrays(reference, mobile)?;
    let reference = center(reference);
    let mobile = center(mobile);
    kabsch_centered_rmsd(&reference, &mobile)
}

/// Calculate RMSD between two parsed structures under an exact atom selection.
pub fn get_rmsd(
    reference: &PDB,
    mobile: &PDB,
    model_num: usize,
    residues: &str,
    atoms: AtomSubset,
) -> ArpeggiaResult<Analysis<f64>> {
    let selector = ResidueSelector::parse(residues)?;
    let reference = select_coordinates(reference, model_num, &selector, atoms)?;
    let mobile = select_coordinates(mobile, model_num, &selector, atoms)?;

    // TODO: Future correspondence may use sequence and/or structural alignment
    // before creating these equal-shaped coordinate arrays, similar to PyMOL.
    validate_correspondence(&reference.keys, &mobile.keys)?;
    let rmsd = kabsch_rmsd(&reference.coordinates, &mobile.coordinates)?;
    let mut warnings = reference.warnings;
    warnings.extend(mobile.warnings);
    Ok(Analysis::new(rmsd, warnings))
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct AtomIdentity {
    chain: String,
    residue_number: isize,
    insertion: String,
    residue_name: String,
    atom_name: String,
}

#[derive(Clone, Debug)]
pub(crate) struct SelectedCoordinates {
    pub(crate) keys: Vec<AtomIdentity>,
    pub(crate) coordinates: Vec<[f64; 3]>,
    pub(crate) warnings: Vec<AnalysisWarning>,
}

pub(crate) fn select_coordinates(
    pdb: &PDB,
    model_num: usize,
    selector: &ResidueSelector,
    atoms: AtomSubset,
) -> ArpeggiaResult<SelectedCoordinates> {
    let model = selected_model(pdb, model_num)?;
    selector.validate_chains(model.chains().map(|chain| chain.id()))?;

    let mut selected = Vec::new();
    for chain in model.chains() {
        for residue in chain.residues() {
            if !selector.matches(
                chain.id(),
                residue.serial_number(),
                residue.insertion_code(),
            ) {
                continue;
            }
            let residue_name = residue.name().unwrap_or("");
            let is_protein = residue.resn().is_some();
            let is_cap = matches!(residue_name, "ACE" | "NH2");
            if !(is_protein || is_cap && matches!(atoms, AtomSubset::Heavy | AtomSubset::All)) {
                continue;
            }
            for hierarchy in residue.atoms_with_hierarchy() {
                let atom = hierarchy.atom();
                if !atom_in_subset(atom.name(), atom.element(), atoms) {
                    continue;
                }
                let coordinate = [atom.x(), atom.y(), atom.z()];
                if !coordinate.iter().all(|value| value.is_finite()) {
                    return Err(ArpeggiaError::Calculation(format!(
                        "atom {} {}:{}{} {} has non-finite coordinates",
                        atom.name(),
                        chain.id(),
                        residue.serial_number(),
                        residue.insertion_code().unwrap_or(""),
                        residue_name
                    )));
                }
                selected.push((
                    AtomIdentity {
                        chain: chain.id().to_string(),
                        residue_number: residue.serial_number(),
                        insertion: residue.insertion_code().unwrap_or("").to_string(),
                        residue_name: residue_name.to_string(),
                        atom_name: atom.name().to_string(),
                    },
                    coordinate,
                ));
            }
        }
    }
    selected.sort_unstable_by(|left, right| left.0.cmp(&right.0));
    if let Some(duplicate) = selected
        .windows(2)
        .find(|pair| pair[0].0 == pair[1].0)
        .map(|pair| &pair[0].0)
    {
        return Err(ArpeggiaError::Calculation(format!(
            "selected atom identity occurs more than once: {duplicate:?}"
        )));
    }
    if selected.is_empty() {
        return Err(ArpeggiaError::InvalidArgument(
            "structure selection contains no atoms".into(),
        ));
    }

    let mut warnings = Vec::new();
    if model_num == 0 && pdb.model_count() > 1 {
        warnings.push(AnalysisWarning::new(
            WarningCode::ModelSelected,
            format!(
                "selected first model {} from {} coordinate models",
                model.serial_number(),
                pdb.model_count()
            ),
        ));
    }
    let (keys, coordinates) = selected.into_iter().unzip();
    Ok(SelectedCoordinates {
        keys,
        coordinates,
        warnings,
    })
}

fn atom_in_subset(name: &str, element: Option<&Element>, subset: AtomSubset) -> bool {
    match subset {
        AtomSubset::Ca => name == "CA",
        AtomSubset::Backbone => matches!(name, "N" | "CA" | "C" | "O" | "OXT"),
        AtomSubset::Heavy => element != Some(&Element::H),
        AtomSubset::All => true,
    }
}

fn validate_correspondence(
    reference: &[AtomIdentity],
    mobile: &[AtomIdentity],
) -> ArpeggiaResult<()> {
    if reference.len() != mobile.len() {
        return Err(ArpeggiaError::Calculation(format!(
            "atom correspondence mismatch: reference has {} selected atoms but mobile has {}",
            reference.len(),
            mobile.len()
        )));
    }
    if let Some((index, (left, right))) = reference
        .iter()
        .zip(mobile)
        .enumerate()
        .find(|(_, (left, right))| left != right)
    {
        return Err(ArpeggiaError::Calculation(format!(
            "atom correspondence mismatch at position {index}: {left:?} != {right:?}"
        )));
    }
    Ok(())
}

fn validate_coordinate_arrays(reference: &[[f64; 3]], mobile: &[[f64; 3]]) -> ArpeggiaResult<()> {
    if reference.len() != mobile.len() {
        return Err(ArpeggiaError::InvalidArgument(format!(
            "coordinate arrays must have equal lengths, got {} and {}",
            reference.len(),
            mobile.len()
        )));
    }
    if reference.len() < 3 {
        return Err(ArpeggiaError::InvalidArgument(
            "RMSD requires at least three coordinate pairs".into(),
        ));
    }
    if reference
        .iter()
        .chain(mobile)
        .flatten()
        .any(|value| !value.is_finite())
    {
        return Err(ArpeggiaError::InvalidArgument(
            "coordinates must be finite".into(),
        ));
    }
    if !is_non_collinear(reference) || !is_non_collinear(mobile) {
        return Err(ArpeggiaError::InvalidArgument(
            "RMSD requires at least three non-collinear points in each structure".into(),
        ));
    }
    Ok(())
}

fn is_non_collinear(points: &[[f64; 3]]) -> bool {
    let origin = na::Vector3::from(points[0]);
    let vectors = points[1..]
        .iter()
        .map(|point| na::Vector3::from(*point) - origin)
        .collect::<Vec<_>>();
    let scale = vectors.iter().map(na::Vector3::norm).fold(0.0, f64::max);
    if scale == 0.0 {
        return false;
    }
    let tolerance = scale * scale * f64::EPSILON * points.len() as f64 * 32.0;
    vectors.iter().enumerate().any(|(index, left)| {
        vectors[index + 1..]
            .iter()
            .any(|right| left.cross(right).norm() > tolerance)
    })
}

fn center(points: &[[f64; 3]]) -> Vec<[f64; 3]> {
    let n = points.len() as f64;
    let centroid = points.iter().fold([0.0; 3], |mut sum, point| {
        for axis in 0..3 {
            sum[axis] += point[axis];
        }
        sum
    });
    points
        .iter()
        .map(|point| {
            [
                point[0] - centroid[0] / n,
                point[1] - centroid[1] / n,
                point[2] - centroid[2] / n,
            ]
        })
        .collect()
}

pub(crate) fn kabsch_centered_rmsd(
    reference: &[[f64; 3]],
    mobile: &[[f64; 3]],
) -> ArpeggiaResult<f64> {
    let mut covariance = na::Matrix3::zeros();
    for (reference, mobile) in reference.iter().zip(mobile) {
        let reference = na::Vector3::from(*reference);
        let mobile = na::Vector3::from(*mobile);
        covariance += mobile * reference.transpose();
    }
    let svd = covariance.svd(true, true);
    let u = svd.u.ok_or_else(|| {
        ArpeggiaError::Calculation("Kabsch SVD did not return left singular vectors".into())
    })?;
    let v_t = svd.v_t.ok_or_else(|| {
        ArpeggiaError::Calculation("Kabsch SVD did not return right singular vectors".into())
    })?;
    let v = v_t.transpose();
    let mut correction = na::Matrix3::identity();
    correction[(2, 2)] = (v * u.transpose()).determinant().signum();
    let rotation = v * correction * u.transpose();

    let squared_error = reference
        .iter()
        .zip(mobile)
        .map(|(reference, mobile)| {
            let delta = na::Vector3::from(*reference) - rotation * na::Vector3::from(*mobile);
            delta.norm_squared()
        })
        .sum::<f64>();
    Ok((squared_error / reference.len() as f64).sqrt())
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ResidueBound {
    number: isize,
    insertion: Option<String>,
}

impl ResidueBound {
    fn parse(value: &str) -> Option<Self> {
        let digit_start = usize::from(value.starts_with('-'));
        let digit_count = value[digit_start..]
            .bytes()
            .take_while(u8::is_ascii_digit)
            .count();
        if digit_count == 0 {
            return None;
        }
        let number_end = digit_start + digit_count;
        let number = value[..number_end].parse().ok()?;
        let insertion = &value[number_end..];
        if insertion.is_empty() {
            Some(Self {
                number,
                insertion: None,
            })
        } else if insertion.len() == 1 && insertion.bytes().all(|byte| byte.is_ascii_alphanumeric())
        {
            Some(Self {
                number,
                insertion: Some(insertion.to_ascii_uppercase()),
            })
        } else {
            None
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum ResidueSpan {
    One(ResidueBound),
    Range(ResidueBound, ResidueBound),
}

impl ResidueSpan {
    fn parse(value: &str) -> ArpeggiaResult<Self> {
        if let Some(bound) = ResidueBound::parse(value) {
            return Ok(Self::One(bound));
        }
        for (index, byte) in value.bytes().enumerate().skip(1) {
            if byte != b'-' {
                continue;
            }
            let Some(start) = ResidueBound::parse(&value[..index]) else {
                continue;
            };
            let Some(end) = ResidueBound::parse(&value[index + 1..]) else {
                continue;
            };
            if compare_bounds(&start, &end) == Ordering::Greater {
                return Err(ArpeggiaError::InvalidArgument(format!(
                    "residue range starts after it ends: {value}"
                )));
            }
            return Ok(Self::Range(start, end));
        }
        Err(ArpeggiaError::InvalidArgument(format!(
            "invalid residue selection: {value}"
        )))
    }

    fn matches(&self, number: isize, insertion: Option<&str>) -> bool {
        match self {
            Self::One(bound) => {
                number == bound.number
                    && bound
                        .insertion
                        .as_deref()
                        .is_none_or(|expected| insertion == Some(expected))
            }
            Self::Range(start, end) => {
                lower_bound_matches(number, insertion, start)
                    && upper_bound_matches(number, insertion, end)
            }
        }
    }
}

fn compare_bounds(left: &ResidueBound, right: &ResidueBound) -> Ordering {
    left.number
        .cmp(&right.number)
        .then_with(|| left.insertion.cmp(&right.insertion))
}

fn lower_bound_matches(number: isize, insertion: Option<&str>, bound: &ResidueBound) -> bool {
    number > bound.number
        || (number == bound.number
            && bound
                .insertion
                .as_deref()
                .is_none_or(|start| insertion.unwrap_or("") >= start))
}

fn upper_bound_matches(number: isize, insertion: Option<&str>, bound: &ResidueBound) -> bool {
    number < bound.number
        || (number == bound.number
            && bound
                .insertion
                .as_deref()
                .is_none_or(|end| insertion.unwrap_or("") <= end))
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct ResidueSelector {
    chains: BTreeMap<String, Option<Vec<ResidueSpan>>>,
}

impl ResidueSelector {
    pub(crate) fn parse(value: &str) -> ArpeggiaResult<Self> {
        if value.trim().is_empty() {
            return Ok(Self::default());
        }
        let mut chains: BTreeMap<String, Option<Vec<ResidueSpan>>> = BTreeMap::new();
        for clause in value.split(',').map(str::trim) {
            if clause.is_empty() {
                return Err(ArpeggiaError::InvalidArgument(
                    "residue selection contains an empty clause".into(),
                ));
            }
            let (chain, span) = clause
                .split_once(':')
                .map_or((clause, None), |(chain, span)| (chain, Some(span)));
            if chain.is_empty() {
                return Err(ArpeggiaError::InvalidArgument(
                    "residue selection contains an empty chain".into(),
                ));
            }
            match span {
                None => {
                    chains.insert(chain.to_string(), None);
                }
                Some("") => {
                    return Err(ArpeggiaError::InvalidArgument(format!(
                        "chain {chain} has an empty residue selection"
                    )));
                }
                Some(span) => {
                    if chains.get(chain).is_some_and(Option::is_none) {
                        continue;
                    }
                    chains
                        .entry(chain.to_string())
                        .or_insert_with(|| Some(Vec::new()))
                        .as_mut()
                        .expect("entry was initialized with a residue list")
                        .push(ResidueSpan::parse(span)?);
                }
            }
        }
        Ok(Self { chains })
    }

    fn matches(&self, chain: &str, number: isize, insertion: Option<&str>) -> bool {
        if self.chains.is_empty() {
            return true;
        }
        match self.chains.get(chain) {
            Some(None) => true,
            Some(Some(spans)) => spans.iter().any(|span| span.matches(number, insertion)),
            None => false,
        }
    }

    fn validate_chains<'a>(&self, chains: impl Iterator<Item = &'a str>) -> ArpeggiaResult<()> {
        if self.chains.is_empty() {
            return Ok(());
        }
        let present = chains.collect::<BTreeSet<_>>();
        let unknown = self
            .chains
            .keys()
            .filter(|chain| !present.contains(chain.as_str()))
            .cloned()
            .collect::<Vec<_>>();
        if unknown.is_empty() {
            Ok(())
        } else {
            Err(ArpeggiaError::InvalidArgument(format!(
                "unknown chain identifiers: {}",
                unknown.join(",")
            )))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kabsch_removes_translation_and_rotation() {
        let reference = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ];
        let mobile = reference.map(|[x, y, z]| [-y + 10.0, x - 4.0, z + 2.0]);
        let forward = kabsch_rmsd(&reference, &mobile).unwrap();
        let reverse = kabsch_rmsd(&mobile, &reference).unwrap();
        assert!(forward < 1e-12);
        assert!((forward - reverse).abs() < 1e-12);
    }

    #[test]
    fn kabsch_does_not_reflect() {
        let reference = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ];
        let reflected = reference.map(|[x, y, z]| [-x, y, z]);
        assert!(kabsch_rmsd(&reference, &reflected).unwrap() > 0.1);
    }

    #[test]
    fn kabsch_rejects_invalid_arrays() {
        assert!(kabsch_rmsd(&[[0.0; 3]; 2], &[[0.0; 3]; 2]).is_err());
        assert!(
            kabsch_rmsd(
                &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
                &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            )
            .is_err()
        );
    }

    #[test]
    fn residue_selector_supports_multichain_ranges() {
        let selector = ResidueSelector::parse("A:1-100,A:110-120,B,C:1,C:3,C:9-20").unwrap();
        assert!(selector.matches("A", 50, None));
        assert!(!selector.matches("A", 105, None));
        assert!(selector.matches("B", -999, None));
        assert!(selector.matches("C", 3, None));
        assert!(!selector.matches("C", 4, None));
        assert!(selector.matches("C", 12, None));
    }

    #[test]
    fn residue_selector_supports_negative_numbers_and_insertions() {
        let selector = ResidueSelector::parse("A:-5--1,B:10A-20").unwrap();
        assert!(selector.matches("A", -3, None));
        assert!(!selector.matches("A", 1, None));
        assert!(!selector.matches("B", 10, None));
        assert!(selector.matches("B", 10, Some("A")));
        assert!(selector.matches("B", 20, Some("B")));
    }

    #[test]
    fn identical_structure_has_zero_rmsd() {
        let input = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
        let pdb = crate::load_model(&input).unwrap().value;
        let analysis = get_rmsd(&pdb, &pdb, 0, "A:1-20", AtomSubset::Ca).unwrap();
        assert!(analysis.value < 1e-12);
    }
}
