//! Rigid protein-structure superposition and atom correspondence.

use crate::contacts::residues::ResidueExt;
use crate::structure::{select_conformers, selected_model};
use crate::{Analysis, AnalysisWarning, ArpeggiaError, ArpeggiaResult, WarningCode};
use nalgebra as na;
use pdbtbx::{ContainsAtomConformer, Element, PDB};
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
    if reference.len() != mobile.len() {
        return Err(ArpeggiaError::InvalidArgument(format!(
            "coordinate arrays must have equal lengths, got {} and {}",
            reference.len(),
            mobile.len()
        )));
    }
    let reference = prepare_coordinates(reference)?;
    let mobile = prepare_coordinates(mobile)?;
    kabsch_prepared_rmsd(&reference, &mobile)
}

/// Calculate RMSD between two parsed structures under an exact atom selection.
///
/// The structures are consumed so alternate conformers can be selected without cloning.
pub fn get_rmsd(
    mut reference: PDB,
    mut mobile: PDB,
    model_num: usize,
    residues: &str,
    atoms: AtomSubset,
) -> ArpeggiaResult<Analysis<f64>> {
    let mut warnings = select_conformers(&mut reference);
    warnings.extend(select_conformers(&mut mobile));
    let selector = ResidueSelector::parse(residues)?;
    let reference = select_coordinates(&reference, model_num, &selector, atoms)?;
    let mobile = select_coordinates(&mobile, model_num, &selector, atoms)?;

    // TODO: Future correspondence may use sequence and/or structural alignment
    // before creating these equal-shaped coordinate arrays, similar to PyMOL.
    validate_correspondence(&reference.keys, &mobile.keys)?;
    let rmsd = kabsch_rmsd(&reference.coordinates, &mobile.coordinates)?;
    warnings.extend(reference.warnings);
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
            // TODO: Add an exact-correspondence all-atom mode that retains arbitrary
            // polymers, ligands, and modified residues when both inputs match.
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
        AtomSubset::Heavy => {
            element != Some(&Element::H) && (element.is_some() || !is_hydrogen_atom_name(name))
        }
        AtomSubset::All => true,
    }
}

fn is_hydrogen_atom_name(name: &str) -> bool {
    matches!(
        name.trim_start_matches(|character: char| character.is_ascii_digit())
            .as_bytes()
            .first(),
        Some(b'H' | b'D' | b'h' | b'd')
    )
}

pub(crate) fn validate_correspondence(
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

pub(crate) struct PreparedCoordinates {
    points: Vec<[f64; 3]>,
    centroid: [f64; 3],
    scale: f64,
}

impl PreparedCoordinates {
    pub(crate) fn len(&self) -> usize {
        self.points.len()
    }
}

pub(crate) fn prepare_coordinates(coordinates: &[[f64; 3]]) -> ArpeggiaResult<PreparedCoordinates> {
    if coordinates.len() < 3 {
        return Err(ArpeggiaError::InvalidArgument(
            "RMSD requires at least three coordinate pairs".into(),
        ));
    }
    if coordinates.iter().flatten().any(|value| !value.is_finite()) {
        return Err(ArpeggiaError::InvalidArgument(
            "coordinates must be finite".into(),
        ));
    }
    let (points, scale) = normalized_displacements(coordinates)
        .ok_or_else(|| ArpeggiaError::Calculation("coordinates have no finite scale".into()))?;
    if !is_non_collinear(&points) {
        return Err(ArpeggiaError::InvalidArgument(
            "RMSD requires at least three non-collinear points in each structure".into(),
        ));
    }
    let n = coordinates.len() as f64;
    let centroid = points.iter().fold([0.0; 3], |mut sum, point| {
        for axis in 0..3 {
            sum[axis] += point[axis] / n;
        }
        sum
    });
    if centroid.iter().any(|value| !value.is_finite()) {
        return Err(ArpeggiaError::Calculation(
            "coordinate centering produced a non-finite value".into(),
        ));
    }
    Ok(PreparedCoordinates {
        points,
        centroid,
        scale,
    })
}

fn is_non_collinear(displacements: &[[f64; 3]]) -> bool {
    let (baseline, scale) = displacements[1..]
        .iter()
        .map(|point| {
            let vector = na::Vector3::from(*point);
            let norm = vector.norm();
            (vector, norm)
        })
        .max_by(|left, right| left.1.total_cmp(&right.1))
        .expect("coordinate arrays contain at least three points");
    if scale == 0.0 {
        return false;
    }
    let tolerance = scale * scale * f64::EPSILON * displacements.len() as f64 * 32.0;
    displacements[1..]
        .iter()
        .any(|point| baseline.cross(&na::Vector3::from(*point)).norm() > tolerance)
}

fn normalized_displacements(points: &[[f64; 3]]) -> Option<(Vec<[f64; 3]>, f64)> {
    let anchor = points[0];
    let mut displacements = points
        .iter()
        .map(|point| {
            [
                point[0] - anchor[0],
                point[1] - anchor[1],
                point[2] - anchor[2],
            ]
        })
        .collect::<Vec<_>>();
    let mut scale = displacements
        .iter()
        .flatten()
        .map(|value| value.abs())
        .max_by(f64::total_cmp)?;
    if !scale.is_finite() {
        scale = points
            .iter()
            .flatten()
            .map(|value| value.abs())
            .max_by(f64::total_cmp)?;
        if scale == 0.0 {
            return None;
        }
        for (displacement, point) in displacements.iter_mut().zip(points) {
            for axis in 0..3 {
                displacement[axis] = point[axis] / scale - anchor[axis] / scale;
            }
        }
    } else if scale > 0.0 {
        for displacement in &mut displacements {
            for value in displacement {
                *value /= scale;
            }
        }
    } else {
        return None;
    }
    Some((displacements, scale))
}

pub(crate) fn kabsch_prepared_rmsd(
    reference: &PreparedCoordinates,
    mobile: &PreparedCoordinates,
) -> ArpeggiaResult<f64> {
    if reference.scale == mobile.scale && reference.points == mobile.points {
        return Ok(0.0);
    }
    let scale = reference.scale.max(mobile.scale);
    let reference_factor = reference.scale / scale;
    let mobile_factor = mobile.scale / scale;
    let reference_centroid = na::Vector3::from(reference.centroid) * reference_factor;
    let mobile_centroid = na::Vector3::from(mobile.centroid) * mobile_factor;
    let mut covariance = na::Matrix3::zeros();
    for (reference, mobile) in reference.points.iter().zip(&mobile.points) {
        let reference = na::Vector3::from(*reference) * reference_factor - reference_centroid;
        let mobile = na::Vector3::from(*mobile) * mobile_factor - mobile_centroid;
        covariance += mobile * reference.transpose();
    }
    let svd = covariance
        .try_svd(true, true, f64::EPSILON * 5.0, 100)
        .ok_or_else(|| ArpeggiaError::Calculation("Kabsch SVD failed to converge".into()))?;
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

    let fitted_residual = aligned_residual_norm(
        reference,
        reference_factor,
        mobile,
        mobile_factor,
        &rotation,
    );
    let identity_residual = aligned_residual_norm(
        reference,
        reference_factor,
        mobile,
        mobile_factor,
        &na::Matrix3::identity(),
    );
    let (residual_norm, fitted) = if fitted_residual < identity_residual {
        (fitted_residual, true)
    } else {
        (identity_residual, false)
    };
    let rmsd = residual_norm * (scale / (reference.points.len() as f64).sqrt());
    let solver_tolerance = f64::EPSILON * (reference.points.len() as f64).sqrt() * 64.0;
    let unreliable_fitted_residual = fitted && residual_norm <= solver_tolerance;
    let unreliable_identity_rotation =
        !fitted && residual_norm <= solver_tolerance && rmsd > 1e-6 && {
            let rotation = identity_rotational_residual_norm(
                reference,
                reference_factor,
                mobile,
                mobile_factor,
            );
            rotation <= solver_tolerance
                && rotation * (scale / (reference.points.len() as f64).sqrt()) > 1e-6
        };
    if (unreliable_fitted_residual && rmsd > 1e-6) || unreliable_identity_rotation {
        return Err(ArpeggiaError::Calculation(
            "coordinate scale prevents a reliable Kabsch residual in Angstroms".into(),
        ));
    }
    if rmsd.is_finite() {
        Ok(rmsd)
    } else {
        Err(ArpeggiaError::Calculation(
            "Kabsch RMSD produced a non-finite result".into(),
        ))
    }
}

fn aligned_residual_norm(
    reference: &PreparedCoordinates,
    reference_factor: f64,
    mobile: &PreparedCoordinates,
    mobile_factor: f64,
    rotation: &na::Matrix3<f64>,
) -> f64 {
    let n = reference.points.len() as f64;
    let residual_centroid = reference.points.iter().zip(&mobile.points).fold(
        na::Vector3::zeros(),
        |sum, (reference, mobile)| {
            let reference = na::Vector3::from(*reference) * reference_factor;
            let mobile = na::Vector3::from(*mobile) * mobile_factor;
            sum + (reference - rotation * mobile) / n
        },
    );
    reference
        .points
        .iter()
        .zip(&mobile.points)
        .fold(0.0_f64, |norm, (reference, mobile)| {
            let reference = na::Vector3::from(*reference) * reference_factor;
            let mobile = na::Vector3::from(*mobile) * mobile_factor;
            let delta = reference - rotation * mobile - residual_centroid;
            delta.iter().fold(norm, |norm, value| norm.hypot(*value))
        })
}

fn identity_rotational_residual_norm(
    reference: &PreparedCoordinates,
    reference_factor: f64,
    mobile: &PreparedCoordinates,
    mobile_factor: f64,
) -> f64 {
    let n = reference.points.len() as f64;
    let residual_centroid = reference.points.iter().zip(&mobile.points).fold(
        na::Vector3::zeros(),
        |sum, (reference, mobile)| {
            let reference = na::Vector3::from(*reference) * reference_factor;
            let mobile = na::Vector3::from(*mobile) * mobile_factor;
            sum + (reference - mobile) / n
        },
    );
    let mobile_centroid = na::Vector3::from(mobile.centroid) * mobile_factor;
    let (inertia, torque) = reference.points.iter().zip(&mobile.points).fold(
        (na::Matrix3::zeros(), na::Vector3::zeros()),
        |(inertia, torque), (reference, mobile)| {
            let reference = na::Vector3::from(*reference) * reference_factor;
            let mobile = na::Vector3::from(*mobile) * mobile_factor;
            let centered_mobile = mobile - mobile_centroid;
            let residual = reference - mobile - residual_centroid;
            (
                inertia + na::Matrix3::identity() * centered_mobile.norm_squared()
                    - centered_mobile * centered_mobile.transpose(),
                torque + centered_mobile.cross(&residual),
            )
        },
    );
    let Some(rotation) = inertia.lu().solve(&torque) else {
        return 0.0;
    };
    mobile.points.iter().fold(0.0_f64, |norm, mobile| {
        let mobile = na::Vector3::from(*mobile) * mobile_factor - mobile_centroid;
        rotation
            .cross(&mobile)
            .iter()
            .fold(norm, |norm, value| norm.hypot(*value))
    })
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
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

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct ResidueSpan {
    start: ResidueBound,
    end: ResidueBound,
}

impl ResidueSpan {
    fn parse(value: &str) -> ArpeggiaResult<Self> {
        if let Some(bound) = ResidueBound::parse(value) {
            return Ok(Self {
                start: bound.clone(),
                end: bound,
            });
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
            if start.number > end.number
                || (start.number == end.number
                    && end.insertion.as_ref().is_some_and(|end| {
                        start.insertion.as_ref().is_some_and(|start| start > end)
                    }))
            {
                return Err(ArpeggiaError::InvalidArgument(format!(
                    "residue range starts after it ends: {value}"
                )));
            }
            return Ok(Self { start, end });
        }
        Err(ArpeggiaError::InvalidArgument(format!(
            "invalid residue selection: {value}"
        )))
    }

    fn matches(&self, number: isize, insertion: Option<&str>) -> bool {
        lower_bound_matches(number, insertion, &self.start)
            && upper_bound_matches(number, insertion, &self.end)
    }
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
        for spans in chains.values_mut().flatten() {
            spans.sort_unstable();
            let mut merged = Vec::<ResidueSpan>::with_capacity(spans.len());
            for span in std::mem::take(spans) {
                if let Some(previous) = merged.last_mut()
                    && upper_bound_matches(
                        span.start.number,
                        span.start.insertion.as_deref(),
                        &previous.end,
                    )
                {
                    if upper_bound_precedes(&previous.end, &span.end) {
                        previous.end = span.end;
                    }
                } else {
                    merged.push(span);
                }
            }
            *spans = merged;
        }
        Ok(Self { chains })
    }

    fn matches(&self, chain: &str, number: isize, insertion: Option<&str>) -> bool {
        if self.chains.is_empty() {
            return true;
        }
        match self.chains.get(chain) {
            Some(None) => true,
            Some(Some(spans)) => {
                let index = spans
                    .partition_point(|span| lower_bound_matches(number, insertion, &span.start));
                index > 0 && spans[index - 1].matches(number, insertion)
            }
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

/// Validate the chain and author-residue selection grammar without structure I/O.
pub fn validate_residue_selection(value: &str) -> ArpeggiaResult<()> {
    ResidueSelector::parse(value).map(|_| ())
}

fn upper_bound_precedes(left: &ResidueBound, right: &ResidueBound) -> bool {
    left.number < right.number
        || left.number == right.number
            && match (&left.insertion, &right.insertion) {
                (None, _) => false,
                (Some(_), None) => true,
                (Some(left), Some(right)) => left < right,
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
    fn kabsch_is_stable_for_coplanar_points_and_large_offsets() {
        let reference = [
            [1_000_000.0, -2_000_000.0, 3_000_000.0],
            [1_000_001.0, -2_000_000.0, 3_000_000.0],
            [1_000_000.0, -1_999_998.0, 3_000_000.0],
            [1_000_002.0, -1_999_997.0, 3_000_000.0],
        ];
        let mobile = reference.map(|[x, y, z]| [-y + 7.0, x - 11.0, z + 5.0]);
        assert!(kabsch_rmsd(&reference, &mobile).unwrap() < 1e-9);
    }

    #[test]
    fn kabsch_is_translation_invariant_at_large_offsets() {
        let reference = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ];
        for translation in [1e6, 1e15] {
            let mobile =
                reference.map(|[x, y, z]| [-y + translation, x + translation, z + translation]);
            assert!(kabsch_rmsd(&reference, &mobile).unwrap() < 1e-12);
        }
    }

    #[test]
    fn kabsch_keeps_extreme_finite_coordinates_finite() {
        for scale in [1e-200, 1e200] {
            let reference = [[0.0, 0.0, 0.0], [scale, 0.0, 0.0], [0.0, scale, 0.0]];
            let mobile = reference.map(|[x, y, z]| [-y, x, z]);
            let rmsd = kabsch_rmsd(&reference, &mobile).unwrap();
            assert!(rmsd.is_finite());
            assert!(rmsd / scale < 1e-12);
        }
    }

    #[test]
    fn kabsch_centers_near_f64_max_without_overflow() {
        let magnitude = 1.3e308;
        let points = [
            [-magnitude, 0.0, 0.0],
            [magnitude, magnitude, 0.0],
            [magnitude, -magnitude, 0.0],
            [magnitude, 0.0, magnitude],
        ];
        assert_eq!(kabsch_rmsd(&points, &points).unwrap(), 0.0);
    }

    #[test]
    fn kabsch_noise_is_positive_and_symmetric() {
        let reference = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ];
        let mut noisy = reference;
        noisy[3][2] += 0.1;
        let forward = kabsch_rmsd(&reference, &noisy).unwrap();
        let reverse = kabsch_rmsd(&noisy, &reference).unwrap();
        assert!(forward > 0.0);
        assert!((forward - reverse).abs() < 1e-12);
    }

    #[test]
    fn kabsch_preserves_small_residuals_at_large_coordinate_scales() {
        let magnitude = 1e200;
        let reference = [
            [0.0, 0.0, 0.0],
            [magnitude, 0.0, 0.0],
            [0.0, magnitude, 0.0],
            [0.0, 0.0, magnitude],
        ];
        let mut mobile = reference;
        mobile[1][1] = 1.0;
        mobile[1][2] = 1.0;
        mobile[2][0] = 1.0;
        mobile[2][2] = 1.0;
        mobile[3][0] = 1.0;
        mobile[3][1] = 1.0;
        let expected = (3.0_f64 / 4.0).sqrt();
        let forward = kabsch_rmsd(&reference, &mobile).unwrap();
        let reverse = kabsch_rmsd(&mobile, &reference).unwrap();
        assert!(
            (forward - expected).abs() < 1e-12,
            "{forward} != {expected}"
        );
        assert!(
            (reverse - expected).abs() < 1e-12,
            "{reverse} != {expected}"
        );
    }

    #[test]
    fn kabsch_rejects_unreliable_extreme_scale_rotation() {
        let magnitude = 1e200;
        let reference = [
            [0.0, 0.0, 0.0],
            [magnitude, 0.0, 0.0],
            [0.0, magnitude, 0.0],
            [0.0, 0.0, magnitude],
        ];
        let rotated = reference.map(|[x, y, z]| [-y, x, z]);
        let infinitesimally_rotated = [
            [0.0, 0.0, 0.0],
            [magnitude, 1.0, 0.0],
            [-1.0, magnitude, 0.0],
            [0.0, 0.0, magnitude],
        ];
        let mixed_strain_and_rotation = [
            [0.0, 0.0, 0.0],
            [magnitude, 1.0, 0.0],
            [1.0, magnitude, 0.0],
            [0.0, 0.0, magnitude],
        ];
        for mobile in [rotated, infinitesimally_rotated, mixed_strain_and_rotation] {
            assert!(matches!(
                kabsch_rmsd(&reference, &mobile),
                Err(ArpeggiaError::Calculation(message)) if message.contains("reliable Kabsch residual")
            ));
        }
    }

    #[test]
    fn kabsch_does_not_fit_scale() {
        let reference = [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let scaled = reference.map(|point| point.map(|value| value * 2.0));
        assert!((kabsch_rmsd(&reference, &scaled).unwrap() - 1.0).abs() < 1e-12);
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
    fn residue_selector_merges_overlapping_ranges() {
        let selector = ResidueSelector::parse("A:1-100,A:20-40,A:50-120,A:1-100").unwrap();
        assert_eq!(selector.chains["A"].as_ref().unwrap().len(), 1);
        assert!(selector.matches("A", 110, None));
        assert!(!selector.matches("A", 121, None));
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
    fn residue_selector_treats_bare_upper_bound_as_last_insertion() {
        let selector = ResidueSelector::parse("A:10A-10").unwrap();
        assert!(!selector.matches("A", 10, None));
        assert!(selector.matches("A", 10, Some("A")));
        assert!(selector.matches("A", 10, Some("B")));
        assert!(!selector.matches("A", 11, None));
    }

    #[test]
    fn residue_selector_rejects_malformed_and_unknown_chains() {
        for invalid in ["A:", "A:5-1", "A,,B", ":1"] {
            assert!(ResidueSelector::parse(invalid).is_err());
        }
        let selector = ResidueSelector::parse("missing").unwrap();
        assert!(selector.validate_chains(["A", "B"].into_iter()).is_err());
    }

    #[test]
    fn atom_subsets_include_caps_but_exclude_nonprotein_groups() {
        let input = std::env::temp_dir().join(format!(
            "arpeggia-rmsd-atom-subsets-{}.pdb",
            std::process::id()
        ));
        std::fs::write(
            &input,
            "HETATM    1  CH3 ACE A   0       0.000   0.000   0.000  1.00 20.00           C  \n\
             HETATM    2  H1  ACE A   0       0.000   1.000   0.000  1.00 20.00           H  \n\
             ATOM      3  N   ALA A   1       1.000   0.000   0.000  1.00 20.00           N  \n\
             ATOM      4  CA  ALA A   1       0.000   0.000   1.000  1.00 20.00           C  \n\
             HETATM    5  O   HOH A   2       3.000   3.000   3.000  1.00 20.00           O  \n\
             END\n",
        )
        .unwrap();
        let pdb = crate::load_model(input.to_str().unwrap()).unwrap().value;
        let selector = ResidueSelector::default();
        assert_eq!(
            select_coordinates(&pdb, 0, &selector, AtomSubset::Heavy)
                .unwrap()
                .keys
                .len(),
            3
        );
        assert_eq!(
            select_coordinates(&pdb, 0, &selector, AtomSubset::All)
                .unwrap()
                .keys
                .len(),
            4
        );
        assert_eq!(
            select_coordinates(&pdb, 0, &selector, AtomSubset::Backbone)
                .unwrap()
                .keys
                .len(),
            2
        );
    }

    #[test]
    fn heavy_subset_infers_digit_leading_hydrogen_names() {
        let hydrogen = pdbtbx::Atom::new(false, 1, "1HB", 0.0, 0.0, 0.0, 1.0, 20.0, "", 0).unwrap();
        assert!(hydrogen.element().is_none());
        assert!(!atom_in_subset(
            hydrogen.name(),
            hydrogen.element(),
            AtomSubset::Heavy
        ));
        assert!(!atom_in_subset("D1", None, AtomSubset::Heavy));
        assert!(atom_in_subset("CB", None, AtomSubset::Heavy));
    }

    #[test]
    fn exact_correspondence_rejects_missing_or_different_atoms() {
        let atom = AtomIdentity {
            chain: "A".into(),
            residue_number: 1,
            insertion: String::new(),
            residue_name: "ALA".into(),
            atom_name: "CA".into(),
        };
        let mut different = atom.clone();
        different.atom_name = "N".into();
        assert!(validate_correspondence(std::slice::from_ref(&atom), &[]).is_err());
        assert!(validate_correspondence(&[atom], &[different]).is_err());
    }

    #[test]
    fn identical_structure_has_zero_rmsd() {
        let input = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
        let pdb = crate::load_model(&input).unwrap().value;
        let analysis = get_rmsd(pdb.clone(), pdb, 0, "A:1-20", AtomSubset::Ca).unwrap();
        assert!(analysis.value < 1e-12);
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
        let mobile = reference.clone();
        let analysis = get_rmsd(reference, mobile, 0, "", AtomSubset::Ca).unwrap();
        assert_eq!(analysis.value, 0.0);
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
