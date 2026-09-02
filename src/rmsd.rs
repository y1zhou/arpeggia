//! Rigid protein-structure superposition and atom correspondence.

use crate::contacts::one_letter_code;
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

/// Superpose and calculate RMSD between two structures under exact correspondence.
///
/// The structures are consumed so alternate conformers can be selected without cloning.
pub fn get_rmsd(
    mut reference: PDB,
    mut mobile: PDB,
    model_num: usize,
    superpose_residues: &str,
    rmsd_residues: &str,
    atoms: AtomSubset,
) -> ArpeggiaResult<Analysis<f64>> {
    let mut warnings = select_conformers(&mut reference);
    warnings.extend(select_conformers(&mut mobile));
    let superpose_selector = ResidueSelector::parse(superpose_residues)?;
    let rmsd_selector = ResidueSelector::parse(rmsd_residues)?;
    let reference = select_coordinate_union(
        &reference,
        model_num,
        &superpose_selector,
        &rmsd_selector,
        atoms,
    )?;
    let mobile = select_coordinate_union(
        &mobile,
        model_num,
        &superpose_selector,
        &rmsd_selector,
        atoms,
    )?;

    // TODO: Future correspondence may use sequence and/or structural alignment
    // before creating these equal-shaped coordinate arrays, similar to PyMOL.
    validate_selection_correspondence(&reference, &mobile)?;
    let reference_coordinates = prepare_coordinate_union(
        &reference.coordinates,
        reference.superpose_end,
        reference.rmsd_start,
    )?;
    let mobile_coordinates =
        prepare_coordinate_union(&mobile.coordinates, mobile.superpose_end, mobile.rmsd_start)?;
    let rmsd = if superpose_selector == rmsd_selector {
        kabsch_prepared_rmsd(&reference_coordinates, &mobile_coordinates)?
    } else {
        kabsch_prepared_selected_rmsd(&reference_coordinates, &mobile_coordinates)?
    };
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
pub(crate) struct SelectedCoordinateUnion {
    pub(crate) keys: Vec<AtomIdentity>,
    pub(crate) coordinates: Vec<[f64; 3]>,
    pub(crate) superpose_end: usize,
    pub(crate) rmsd_start: usize,
    pub(crate) warnings: Vec<AnalysisWarning>,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum SelectionMembership {
    SuperposeOnly,
    Both,
    RmsdOnly,
}

pub(crate) fn select_coordinate_union(
    pdb: &PDB,
    model_num: usize,
    superpose_selector: &ResidueSelector,
    rmsd_selector: &ResidueSelector,
    atoms: AtomSubset,
) -> ArpeggiaResult<SelectedCoordinateUnion> {
    let model = selected_model(pdb, model_num)?;
    superpose_selector
        .validate_chains(model.chains().map(|chain| chain.id()))
        .map_err(|error| selection_error("Superposition Selection", error))?;
    rmsd_selector
        .validate_chains(model.chains().map(|chain| chain.id()))
        .map_err(|error| selection_error("RMSD Selection", error))?;

    let mut selected = Vec::new();
    for chain in model.chains() {
        for residue in chain.residues() {
            let superpose = superpose_selector.matches(
                chain.id(),
                residue.serial_number(),
                residue.insertion_code(),
            );
            let rmsd = rmsd_selector.matches(
                chain.id(),
                residue.serial_number(),
                residue.insertion_code(),
            );
            let membership = match (superpose, rmsd) {
                (true, true) => SelectionMembership::Both,
                (true, false) => SelectionMembership::SuperposeOnly,
                (false, true) => SelectionMembership::RmsdOnly,
                (false, false) => continue,
            };
            let residue_name = residue.name().unwrap_or("");
            // TODO: Add an exact-correspondence all-atom mode that retains arbitrary
            // polymers, ligands, and modified residues when both inputs match.
            let is_protein = residue.name().and_then(one_letter_code).is_some();
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
                    membership,
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
    selected.sort_unstable_by(|left, right| (&left.0, &left.1).cmp(&(&right.0, &right.1)));
    if let Some(duplicate) = selected
        .windows(2)
        .find(|pair| pair[0].1 == pair[1].1)
        .map(|pair| &pair[0].1)
    {
        return Err(ArpeggiaError::Calculation(format!(
            "selected atom identity occurs more than once: {duplicate:?}"
        )));
    }
    let superpose_end =
        selected.partition_point(|(membership, _, _)| *membership != SelectionMembership::RmsdOnly);
    let rmsd_start = selected
        .partition_point(|(membership, _, _)| *membership == SelectionMembership::SuperposeOnly);
    if superpose_end == 0 {
        return Err(ArpeggiaError::InvalidArgument(
            "Superposition Selection contains no atoms".into(),
        ));
    }
    if rmsd_start == selected.len() {
        return Err(ArpeggiaError::InvalidArgument(
            "RMSD Selection contains no atoms".into(),
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
    let (keys, coordinates) = selected
        .into_iter()
        .map(|(_, identity, coordinate)| (identity, coordinate))
        .unzip();
    Ok(SelectedCoordinateUnion {
        keys,
        coordinates,
        superpose_end,
        rmsd_start,
        warnings,
    })
}

fn selection_error(selection: &str, error: ArpeggiaError) -> ArpeggiaError {
    match error {
        ArpeggiaError::InvalidArgument(message) => {
            ArpeggiaError::InvalidArgument(format!("{selection}: {message}"))
        }
        error => error,
    }
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

pub(crate) fn validate_selection_correspondence(
    reference: &SelectedCoordinateUnion,
    mobile: &SelectedCoordinateUnion,
) -> ArpeggiaResult<()> {
    validate_selection_keys(
        &reference.keys,
        reference.superpose_end,
        reference.rmsd_start,
        &mobile.keys,
        mobile.superpose_end,
        mobile.rmsd_start,
    )
}

pub(crate) fn validate_selection_keys(
    reference: &[AtomIdentity],
    reference_superpose_end: usize,
    reference_rmsd_start: usize,
    mobile: &[AtomIdentity],
    mobile_superpose_end: usize,
    mobile_rmsd_start: usize,
) -> ArpeggiaResult<()> {
    validate_correspondence(
        &reference[..reference_superpose_end],
        &mobile[..mobile_superpose_end],
    )
    .map_err(|error| selection_error("Superposition Selection", error))?;
    validate_correspondence(
        &reference[reference_rmsd_start..],
        &mobile[mobile_rmsd_start..],
    )
    .map_err(|error| selection_error("RMSD Selection", error))
}

pub(crate) struct PreparedCoordinates {
    points: Vec<[f64; 3]>,
    centroid: [f64; 3],
    scale: f64,
    superpose_end: usize,
    rmsd_start: usize,
}

impl PreparedCoordinates {
    pub(crate) fn len(&self) -> usize {
        self.points.len()
    }
}

pub(crate) fn prepare_coordinates(coordinates: &[[f64; 3]]) -> ArpeggiaResult<PreparedCoordinates> {
    prepare_coordinate_union(coordinates, coordinates.len(), 0)
}

pub(crate) fn prepare_coordinate_union(
    coordinates: &[[f64; 3]],
    superpose_end: usize,
    rmsd_start: usize,
) -> ArpeggiaResult<PreparedCoordinates> {
    if superpose_end > coordinates.len() || rmsd_start > superpose_end {
        return Err(ArpeggiaError::Calculation(
            "invalid coordinate-selection boundaries".into(),
        ));
    }
    if superpose_end < 3 {
        return Err(ArpeggiaError::InvalidArgument(
            "Superposition Selection requires at least three coordinate pairs".into(),
        ));
    }
    if rmsd_start == coordinates.len() {
        return Err(ArpeggiaError::InvalidArgument(
            "RMSD Selection requires at least one coordinate pair".into(),
        ));
    }
    if coordinates.iter().flatten().any(|value| !value.is_finite()) {
        return Err(ArpeggiaError::InvalidArgument(
            "coordinates must be finite".into(),
        ));
    }
    let (points, scale) = normalized_displacements(coordinates, superpose_end)
        .ok_or_else(|| ArpeggiaError::Calculation("coordinates have no finite scale".into()))?;
    if points.iter().flatten().any(|value| !value.is_finite()) {
        return Err(ArpeggiaError::Calculation(
            "RMSD Selection cannot be represented in the Superposition Selection coordinate frame"
                .into(),
        ));
    }
    if !is_non_collinear(&points[..superpose_end]) {
        return Err(ArpeggiaError::InvalidArgument(
            "Superposition Selection requires at least three non-collinear points in each structure"
                .into(),
        ));
    }
    let n = superpose_end as f64;
    let centroid = points[..superpose_end]
        .iter()
        .fold([0.0; 3], |mut sum, point| {
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
        superpose_end,
        rmsd_start,
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

fn normalized_displacements(
    points: &[[f64; 3]],
    superpose_end: usize,
) -> Option<(Vec<[f64; 3]>, f64)> {
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
    let mut scale = displacements[..superpose_end]
        .iter()
        .flatten()
        .map(|value| value.abs())
        .max_by(f64::total_cmp)?;
    if displacements
        .iter()
        .flatten()
        .any(|value| !value.is_finite())
    {
        scale = points[..superpose_end]
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
    let transform = fit_prepared_transform(reference, mobile)?;
    finish_rmsd(
        transform.fit_residual_norm,
        transform.scale,
        reference.superpose_end,
    )
}

pub(crate) fn kabsch_prepared_selected_rmsd(
    reference: &PreparedCoordinates,
    mobile: &PreparedCoordinates,
) -> ArpeggiaResult<f64> {
    if reference.scale == mobile.scale && reference.points == mobile.points {
        return Ok(0.0);
    }
    let transform = fit_prepared_transform(reference, mobile)?;
    let reference_rmsd = &reference.points[reference.rmsd_start..];
    let mobile_rmsd = &mobile.points[mobile.rmsd_start..];
    if reference_rmsd.len() != mobile_rmsd.len() {
        return Err(ArpeggiaError::Calculation(format!(
            "RMSD Selection coordinate mismatch: reference has {} atoms but mobile has {}",
            reference_rmsd.len(),
            mobile_rmsd.len()
        )));
    }
    let residual_norm = fixed_transform_residual_norm(
        reference_rmsd,
        transform.reference_factor,
        mobile_rmsd,
        transform.mobile_factor,
        &transform.rotation,
        &transform.residual_centroid,
    );
    let scoring_radius = rms_radius(
        mobile_rmsd,
        transform.mobile_factor,
        na::Vector3::from(mobile.centroid) * transform.mobile_factor,
    );
    if scoring_radius > 0.0
        && (!transform.angular_error_bound.is_finite()
            || transform.angular_error_bound * scoring_radius * transform.scale > 1e-6)
    {
        return Err(ArpeggiaError::Calculation(
            "coordinate scale prevents a reliable Kabsch residual in Angstroms".into(),
        ));
    }
    finish_rmsd(residual_norm, transform.scale, reference_rmsd.len())
}

struct PreparedTransform {
    rotation: na::Matrix3<f64>,
    residual_centroid: na::Vector3<f64>,
    reference_factor: f64,
    mobile_factor: f64,
    scale: f64,
    fit_residual_norm: f64,
    angular_error_bound: f64,
}

fn fit_prepared_transform(
    reference: &PreparedCoordinates,
    mobile: &PreparedCoordinates,
) -> ArpeggiaResult<PreparedTransform> {
    let reference_points = &reference.points[..reference.superpose_end];
    let mobile_points = &mobile.points[..mobile.superpose_end];
    if reference_points.len() != mobile_points.len() {
        return Err(ArpeggiaError::Calculation(format!(
            "Superposition Selection coordinate mismatch: reference has {} atoms but mobile has {}",
            reference_points.len(),
            mobile_points.len()
        )));
    }
    let scale = reference.scale.max(mobile.scale);
    let reference_factor = reference.scale / scale;
    let mobile_factor = mobile.scale / scale;
    let reference_centroid = na::Vector3::from(reference.centroid) * reference_factor;
    let mobile_centroid = na::Vector3::from(mobile.centroid) * mobile_factor;
    let mut covariance = na::Matrix3::zeros();
    for (reference, mobile) in reference_points.iter().zip(mobile_points) {
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
    // [WARNING] Reflection-corrected fits need signed singular-value stiffness
    // before this bound can diagnose an ambiguous proper rotation.
    let minimum_rotational_inertia = svd.singular_values.sum() - svd.singular_values.max();
    let v = v_t.transpose();
    let mut correction = na::Matrix3::identity();
    let correction_sign = (v * u.transpose()).determinant().signum();
    correction[(2, 2)] = correction_sign;
    let rotation = v * correction * u.transpose();

    let (fitted_residual, fitted_centroid) = aligned_residual_norm(
        reference_points,
        reference_factor,
        mobile_points,
        mobile_factor,
        &rotation,
    );
    let identity = na::Matrix3::identity();
    let (identity_residual, identity_centroid) = aligned_residual_norm(
        reference_points,
        reference_factor,
        mobile_points,
        mobile_factor,
        &identity,
    );
    let (residual_norm, fitted, rotation, residual_centroid) =
        if fitted_residual < identity_residual {
            (fitted_residual, true, rotation, fitted_centroid)
        } else {
            (identity_residual, false, identity, identity_centroid)
        };
    let rmsd = scaled_rmsd(residual_norm, scale, reference_points.len());
    let solver_tolerance = f64::EPSILON * (reference_points.len() as f64).sqrt() * 64.0;
    let unreliable_fitted_residual = fitted && residual_norm <= solver_tolerance;
    let identity_rotation = if fitted {
        None
    } else {
        identity_rotational_residual(
            reference_points,
            reference_factor,
            mobile_points,
            mobile.centroid,
            mobile_factor,
        )
    };
    let unreliable_identity_rotation =
        !fitted && residual_norm <= solver_tolerance && rmsd > 1e-6 && {
            let rotation = identity_rotation.map_or(0.0, |(_, norm)| norm);
            rotation <= solver_tolerance
                && scaled_rmsd(rotation, scale, reference_points.len()) > 1e-6
        };
    if (unreliable_fitted_residual && rmsd > 1e-6) || unreliable_identity_rotation {
        return Err(ArpeggiaError::Calculation(
            "coordinate scale prevents a reliable Kabsch residual in Angstroms".into(),
        ));
    }
    let solver_angular_error = solver_tolerance / minimum_rotational_inertia.max(0.0).sqrt();
    let angular_error_bound = if correction_sign < 0.0 {
        f64::INFINITY
    } else if fitted {
        solver_angular_error
    } else {
        identity_rotation.map_or(f64::INFINITY, |(angle, _)| solver_angular_error + angle)
    };
    Ok(PreparedTransform {
        rotation,
        residual_centroid,
        reference_factor,
        mobile_factor,
        scale,
        fit_residual_norm: residual_norm,
        angular_error_bound,
    })
}

fn rms_radius(points: &[[f64; 3]], factor: f64, centroid: na::Vector3<f64>) -> f64 {
    points.iter().fold(0.0_f64, |norm, point| {
        (na::Vector3::from(*point) * factor - centroid)
            .iter()
            .fold(norm, |norm, value| norm.hypot(*value))
    }) / (points.len() as f64).sqrt()
}

fn finish_rmsd(residual_norm: f64, scale: f64, count: usize) -> ArpeggiaResult<f64> {
    let rmsd = scaled_rmsd(residual_norm, scale, count);
    if rmsd.is_finite() {
        Ok(rmsd)
    } else {
        Err(ArpeggiaError::Calculation(
            "Kabsch RMSD produced a non-finite result".into(),
        ))
    }
}

fn scaled_rmsd(residual_norm: f64, scale: f64, count: usize) -> f64 {
    residual_norm / (count as f64).sqrt() * scale
}

fn aligned_residual_norm(
    reference: &[[f64; 3]],
    reference_factor: f64,
    mobile: &[[f64; 3]],
    mobile_factor: f64,
    rotation: &na::Matrix3<f64>,
) -> (f64, na::Vector3<f64>) {
    let n = reference.len() as f64;
    let residual_centroid =
        reference
            .iter()
            .zip(mobile)
            .fold(na::Vector3::zeros(), |sum, (reference, mobile)| {
                let reference = na::Vector3::from(*reference) * reference_factor;
                let mobile = na::Vector3::from(*mobile) * mobile_factor;
                sum + (reference - rotation * mobile) / n
            });
    let norm = fixed_transform_residual_norm(
        reference,
        reference_factor,
        mobile,
        mobile_factor,
        rotation,
        &residual_centroid,
    );
    (norm, residual_centroid)
}

fn fixed_transform_residual_norm(
    reference: &[[f64; 3]],
    reference_factor: f64,
    mobile: &[[f64; 3]],
    mobile_factor: f64,
    rotation: &na::Matrix3<f64>,
    residual_centroid: &na::Vector3<f64>,
) -> f64 {
    reference
        .iter()
        .zip(mobile)
        .fold(0.0_f64, |norm, (reference, mobile)| {
            let reference = na::Vector3::from(*reference) * reference_factor;
            let mobile = na::Vector3::from(*mobile) * mobile_factor;
            let delta = reference - rotation * mobile - residual_centroid;
            delta.iter().fold(norm, |norm, value| norm.hypot(*value))
        })
}

fn identity_rotational_residual(
    reference: &[[f64; 3]],
    reference_factor: f64,
    mobile: &[[f64; 3]],
    mobile_centroid: [f64; 3],
    mobile_factor: f64,
) -> Option<(f64, f64)> {
    let n = reference.len() as f64;
    let residual_centroid =
        reference
            .iter()
            .zip(mobile)
            .fold(na::Vector3::zeros(), |sum, (reference, mobile)| {
                let reference = na::Vector3::from(*reference) * reference_factor;
                let mobile = na::Vector3::from(*mobile) * mobile_factor;
                sum + (reference - mobile) / n
            });
    let mobile_centroid = na::Vector3::from(mobile_centroid) * mobile_factor;
    let (inertia, torque) = reference.iter().zip(mobile).fold(
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
    let rotation = inertia.lu().solve(&torque)?;
    let residual_norm = mobile.iter().fold(0.0_f64, |norm, mobile| {
        let mobile = na::Vector3::from(*mobile) * mobile_factor - mobile_centroid;
        rotation
            .cross(&mobile)
            .iter()
            .fold(norm, |norm, value| norm.hypot(*value))
    });
    Some((rotation.norm(), residual_norm))
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
                    && (upper_bound_matches(
                        span.start.number,
                        span.start.insertion.as_deref(),
                        &previous.end,
                    ) || previous.end.insertion.is_none()
                        && span.start.insertion.is_none()
                        && previous.end.number.checked_add(1) == Some(span.start.number))
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

/// Validate both RMSD residue-selection grammars without structure I/O.
pub fn validate_rmsd_selections(
    superpose_residues: &str,
    rmsd_residues: &str,
) -> ArpeggiaResult<()> {
    ResidueSelector::parse(superpose_residues)
        .map_err(|error| selection_error("Superposition Selection", error))?;
    ResidueSelector::parse(rmsd_residues)
        .map_err(|error| selection_error("RMSD Selection", error))?;
    Ok(())
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
    fn generalized_equal_selection_matches_prepared_kabsch_bit_for_bit() {
        let reference = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ];
        let mut mobile = reference.map(|[x, y, z]| [-y + 10.0, x - 4.0, z + 2.0]);
        mobile[3][0] += 0.125;
        let reference = prepare_coordinates(&reference).unwrap();
        let mobile = prepare_coordinates(&mobile).unwrap();
        assert_eq!(
            kabsch_prepared_selected_rmsd(&reference, &mobile)
                .unwrap()
                .to_bits(),
            kabsch_prepared_rmsd(&reference, &mobile).unwrap().to_bits()
        );
    }

    #[test]
    fn fixed_superposition_transform_scores_a_disjoint_single_atom() {
        let reference = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mobile = [
            [10.0, -4.0, 2.0],
            [10.0, -3.0, 2.0],
            [9.0, -4.0, 2.0],
            [10.0, -4.0, 5.0],
        ];
        let reference = prepare_coordinate_union(&reference, 3, 3).unwrap();
        let mobile = prepare_coordinate_union(&mobile, 3, 3).unwrap();
        let rmsd = kabsch_prepared_selected_rmsd(&reference, &mobile).unwrap();
        let reverse = kabsch_prepared_selected_rmsd(&mobile, &reference).unwrap();
        assert!((rmsd - 2.0).abs() < 1e-12);
        assert!((rmsd - reverse).abs() < 1e-12);
    }

    #[test]
    fn distant_rmsd_atoms_must_not_amplify_solver_scale_rotation_error() {
        let distance = 1e12;
        let reference = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
            [distance, 0.3 * distance, -0.2 * distance],
        ];
        let mobile = reference.map(|[x, y, z]| [-y + 10.0, x - 4.0, z + 2.0]);
        let reference = prepare_coordinate_union(&reference, 4, 4).unwrap();
        let mobile = prepare_coordinate_union(&mobile, 4, 4).unwrap();
        assert!(matches!(
            kabsch_prepared_selected_rmsd(&reference, &mobile),
            Err(ArpeggiaError::Calculation(message))
                if message.contains("reliable Kabsch residual")
        ));
    }

    #[test]
    fn distant_rmsd_atoms_reject_identity_fallback_rotation_error() {
        let distance = 1e12;
        let angle = 1e-17;
        let reference = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
            [distance, 0.0, 0.0],
        ];
        let mobile =
            reference.map(|[x, y, z]| [x - angle * y + 10.0, angle * x + y - 4.0, z + 2.0]);
        let reference = prepare_coordinate_union(&reference, 4, 4).unwrap();
        let mobile = prepare_coordinate_union(&mobile, 4, 4).unwrap();
        let transform = fit_prepared_transform(&reference, &mobile).unwrap();
        assert_eq!(transform.rotation, na::Matrix3::identity());
        assert!(matches!(
            kabsch_prepared_selected_rmsd(&reference, &mobile),
            Err(ArpeggiaError::Calculation(message))
                if message.contains("reliable Kabsch residual")
        ));
    }

    #[test]
    fn distant_rmsd_atoms_reject_solver_error_with_fit_residual() {
        let distance = 1e12;
        let reference = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
            [distance, 0.3 * distance, -0.2 * distance],
        ];
        let mut mobile = reference.map(|[x, y, z]| [-y + 10.0, x - 4.0, z + 2.0]);
        mobile[3][0] += 0.01;
        let reference = prepare_coordinate_union(&reference, 4, 4).unwrap();
        let mobile = prepare_coordinate_union(&mobile, 4, 4).unwrap();
        let transform = fit_prepared_transform(&reference, &mobile).unwrap();
        let solver_tolerance = f64::EPSILON * 4.0_f64.sqrt() * 64.0;
        assert!(transform.fit_residual_norm > solver_tolerance);
        assert!(matches!(
            kabsch_prepared_selected_rmsd(&reference, &mobile),
            Err(ArpeggiaError::Calculation(message))
                if message.contains("reliable Kabsch residual")
        ));
    }

    #[test]
    fn rmsd_only_coordinates_must_fit_the_superposition_frame() {
        let coordinates = [
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [f64::MAX, 0.0, 0.0],
        ];
        assert!(matches!(
            prepare_coordinate_union(&coordinates, 3, 3),
            Err(ArpeggiaError::Calculation(message))
                if message.contains("RMSD Selection")
        ));
    }

    #[test]
    fn rmsd_only_subtraction_overflow_uses_the_superposition_scale() {
        let magnitude = 1.3e308;
        let coordinates = [
            [magnitude, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [magnitude, magnitude, 0.0],
            [magnitude, 0.0, magnitude],
            [-magnitude, 0.0, 0.0],
        ];
        let prepared = prepare_coordinate_union(&coordinates, 4, 4).unwrap();
        assert_eq!(
            kabsch_prepared_selected_rmsd(&prepared, &prepared).unwrap(),
            0.0
        );
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
    fn rmsd_rescaling_preserves_subnormal_results() {
        let minimum_subnormal = f64::from_bits(1);
        assert_eq!(scaled_rmsd(2.0, minimum_subnormal, 4), minimum_subnormal);
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
    fn residue_selector_merges_only_truly_adjacent_ranges() {
        let adjacent = ResidueSelector::parse("A:1-5,A:6-10").unwrap();
        assert_eq!(adjacent, ResidueSelector::parse("A:1-10").unwrap());

        let insertion_gap = ResidueSelector::parse("A:1-5A,A:6-10").unwrap();
        assert_eq!(insertion_gap.chains["A"].as_ref().unwrap().len(), 2);
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
            select_coordinate_union(&pdb, 0, &selector, &selector, AtomSubset::Heavy)
                .unwrap()
                .keys
                .len(),
            3
        );
        assert_eq!(
            select_coordinate_union(&pdb, 0, &selector, &selector, AtomSubset::All)
                .unwrap()
                .keys
                .len(),
            4
        );
        assert_eq!(
            select_coordinate_union(&pdb, 0, &selector, &selector, AtomSubset::Backbone)
                .unwrap()
                .keys
                .len(),
            2
        );
    }

    #[test]
    fn coordinate_union_uses_three_contiguous_membership_blocks() {
        let input = std::env::temp_dir().join(format!(
            "arpeggia-rmsd-selection-union-{}.pdb",
            std::process::id()
        ));
        std::fs::write(
            &input,
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
             ATOM      2  CA  ALA B   1       1.000   0.000   0.000  1.00 20.00           C  \n\
             ATOM      3  CA  ALA C   1       0.000   1.000   0.000  1.00 20.00           C  \n\
             END\n",
        )
        .unwrap();
        let pdb = crate::load_model(input.to_str().unwrap()).unwrap().value;
        let superpose = ResidueSelector::parse("A,B").unwrap();
        let rmsd = ResidueSelector::parse("B,C").unwrap();
        let selected = select_coordinate_union(&pdb, 0, &superpose, &rmsd, AtomSubset::Ca).unwrap();
        assert_eq!(selected.superpose_end, 2);
        assert_eq!(selected.rmsd_start, 1);
        assert_eq!(
            selected.coordinates.len() * std::mem::size_of::<[f64; 3]>(),
            72
        );
        assert_eq!(
            selected
                .keys
                .iter()
                .map(|identity| identity.chain.as_str())
                .collect::<Vec<_>>(),
            ["A", "B", "C"]
        );

        let missing = ResidueSelector::parse("missing").unwrap();
        assert!(matches!(
            select_coordinate_union(&pdb, 0, &missing, &rmsd, AtomSubset::Ca),
            Err(ArpeggiaError::InvalidArgument(message))
                if message.starts_with("Superposition Selection:")
        ));
        assert!(matches!(
            select_coordinate_union(&pdb, 0, &superpose, &missing, AtomSubset::Ca),
            Err(ArpeggiaError::InvalidArgument(message))
                if message.starts_with("RMSD Selection:")
        ));
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
        let analysis = get_rmsd(pdb.clone(), pdb, 0, "A:1-20", "A:1-20", AtomSubset::Ca).unwrap();
        assert!(analysis.value < 1e-12);
    }

    #[test]
    fn public_rmsd_superposes_one_chain_and_scores_others() {
        let directory = std::env::temp_dir().join(format!(
            "arpeggia-relative-domain-rmsd-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&directory).unwrap();
        let reference_path = directory.join("reference.pdb");
        let mobile_path = directory.join("mobile.pdb");
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
            &mobile_path,
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
        let mobile = crate::load_model(mobile_path.to_str().unwrap())
            .unwrap()
            .value;
        let selected = get_rmsd(
            reference.clone(),
            mobile.clone(),
            0,
            "A",
            "B,C",
            AtomSubset::Ca,
        )
        .unwrap()
        .value;
        assert!((selected - 2.0).abs() < 1e-12);

        let default_all = get_rmsd(reference, mobile, 0, "A", "", AtomSubset::Ca)
            .unwrap()
            .value;
        assert!((default_all - (8.0_f64 / 5.0).sqrt()).abs() < 1e-12);
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
        let analysis = get_rmsd(reference, mobile, 0, "", "", AtomSubset::Ca).unwrap();
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
