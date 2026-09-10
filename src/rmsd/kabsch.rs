//! Coordinate conditioning, rigid fitting, and fixed-transform residuals.

use crate::{ArpeggiaError, ArpeggiaResult};
use nalgebra as na;

/// Calculate optimally superposed RMSD between equal-sized coordinate arrays.
///
/// Coordinates are uniformly weighted. Reflection and scaling are forbidden.
pub fn kabsch_rmsd(reference: &[[f64; 3]], query: &[[f64; 3]]) -> ArpeggiaResult<f64> {
    if reference.len() != query.len() {
        return Err(ArpeggiaError::InvalidArgument(format!(
            "coordinate arrays must have equal lengths, got {} and {}",
            reference.len(),
            query.len()
        )));
    }
    let reference = prepare_coordinates(reference)?;
    let query = prepare_coordinates(query)?;
    kabsch_prepared_rmsd(&reference, &query)
}

pub(super) struct PreparedCoordinates {
    pub(super) points: Vec<[f64; 3]>,
    centroid: [f64; 3],
    pub(super) scale: f64,
    rmsd_radius: f64,
    superpose_end: usize,
    rmsd_start: usize,
}

impl PreparedCoordinates {
    pub(super) fn len(&self) -> usize {
        self.points.len()
    }
}

pub(super) fn prepare_coordinates(coordinates: &[[f64; 3]]) -> ArpeggiaResult<PreparedCoordinates> {
    prepare_coordinate_union(coordinates, coordinates.len(), 0)
}

pub(super) fn prepare_coordinate_union(
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
    let rmsd_radius = scaled_rmsd(
        rms_radius(&points[rmsd_start..], na::Vector3::from(centroid)),
        scale,
        1,
    );
    Ok(PreparedCoordinates {
        points,
        centroid,
        scale,
        rmsd_radius,
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

pub(super) fn kabsch_prepared_rmsd(
    reference: &PreparedCoordinates,
    query: &PreparedCoordinates,
) -> ArpeggiaResult<f64> {
    if reference.scale == query.scale && reference.points == query.points {
        return Ok(0.0);
    }
    let transform = fit_prepared_transform(reference, query)?;
    finish_rmsd(
        transform.fit_residual_norm,
        transform.scale,
        reference.superpose_end,
    )
}

pub(super) fn kabsch_prepared_selected_rmsd(
    reference: &PreparedCoordinates,
    query: &PreparedCoordinates,
) -> ArpeggiaResult<f64> {
    if reference.scale == query.scale && reference.points == query.points {
        return Ok(0.0);
    }
    let transform = fit_prepared_transform(reference, query)?;
    let reference_rmsd = &reference.points[reference.rmsd_start..];
    let query_rmsd = &query.points[query.rmsd_start..];
    if reference_rmsd.len() != query_rmsd.len() {
        return Err(ArpeggiaError::Calculation(format!(
            "RMSD Selection coordinate mismatch: reference has {} atoms but query has {}",
            reference_rmsd.len(),
            query_rmsd.len()
        )));
    }
    let residual_norm = fixed_transform_residual_norm(
        reference_rmsd,
        transform.reference_factor,
        query_rmsd,
        transform.query_factor,
        &transform.rotation,
        &transform.residual_centroid,
    );
    let scoring_radius = query.rmsd_radius;
    if scoring_radius > 0.0
        && (!transform.angular_error_bound.is_finite()
            || transform.angular_error_bound * scoring_radius > 1e-6)
    {
        return Err(ArpeggiaError::Calculation(
            "coordinate scale prevents a reliable Kabsch residual in Angstroms".into(),
        ));
    }
    finish_rmsd(residual_norm, transform.scale, reference_rmsd.len())
}

pub(super) struct PreparedTransform {
    pub(super) rotation: na::Matrix3<f64>,
    pub(super) residual_centroid: na::Vector3<f64>,
    pub(super) reference_factor: f64,
    pub(super) query_factor: f64,
    scale: f64,
    pub(super) fit_residual_norm: f64,
    angular_error_bound: f64,
}

pub(super) fn fit_prepared_transform(
    reference: &PreparedCoordinates,
    query: &PreparedCoordinates,
) -> ArpeggiaResult<PreparedTransform> {
    let reference_points = &reference.points[..reference.superpose_end];
    let query_points = &query.points[..query.superpose_end];
    if reference_points.len() != query_points.len() {
        return Err(ArpeggiaError::Calculation(format!(
            "Superposition Selection coordinate mismatch: reference has {} atoms but query has {}",
            reference_points.len(),
            query_points.len()
        )));
    }
    let scale = reference.scale.max(query.scale);
    let reference_factor = reference.scale / scale;
    let query_factor = query.scale / scale;
    let reference_centroid = na::Vector3::from(reference.centroid) * reference_factor;
    let query_centroid = na::Vector3::from(query.centroid) * query_factor;
    let mut covariance = na::Matrix3::zeros();
    for (reference, query) in reference_points.iter().zip(query_points) {
        let reference = na::Vector3::from(*reference) * reference_factor - reference_centroid;
        let query = na::Vector3::from(*query) * query_factor - query_centroid;
        covariance += query * reference.transpose();
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
        query_points,
        query_factor,
        &rotation,
    );
    let identity = na::Matrix3::identity();
    let (identity_residual, identity_centroid) = aligned_residual_norm(
        reference_points,
        reference_factor,
        query_points,
        query_factor,
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
            query_points,
            query.centroid,
            query_factor,
            identity_centroid,
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
        query_factor,
        scale,
        fit_residual_norm: residual_norm,
        angular_error_bound,
    })
}

fn rms_radius(points: &[[f64; 3]], centroid: na::Vector3<f64>) -> f64 {
    points.iter().fold(0.0_f64, |norm, point| {
        (na::Vector3::from(*point) - centroid)
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
    let divisor = (count as f64).sqrt();
    let normalized = residual_norm / divisor;
    if residual_norm != 0.0 && (normalized == 0.0 || (scale > 1.0 && normalized.is_subnormal())) {
        residual_norm * scale / divisor
    } else {
        normalized * scale
    }
}

fn aligned_residual_norm(
    reference: &[[f64; 3]],
    reference_factor: f64,
    query: &[[f64; 3]],
    query_factor: f64,
    rotation: &na::Matrix3<f64>,
) -> (f64, na::Vector3<f64>) {
    let n = reference.len() as f64;
    let residual_centroid =
        reference
            .iter()
            .zip(query)
            .fold(na::Vector3::zeros(), |sum, (reference, query)| {
                let reference = na::Vector3::from(*reference) * reference_factor;
                let query = na::Vector3::from(*query) * query_factor;
                sum + (reference - rotation * query) / n
            });
    let norm = fixed_transform_residual_norm(
        reference,
        reference_factor,
        query,
        query_factor,
        rotation,
        &residual_centroid,
    );
    (norm, residual_centroid)
}

fn fixed_transform_residual_norm(
    reference: &[[f64; 3]],
    reference_factor: f64,
    query: &[[f64; 3]],
    query_factor: f64,
    rotation: &na::Matrix3<f64>,
    residual_centroid: &na::Vector3<f64>,
) -> f64 {
    reference
        .iter()
        .zip(query)
        .fold(0.0_f64, |norm, (reference, query)| {
            let reference = na::Vector3::from(*reference) * reference_factor;
            let query = na::Vector3::from(*query) * query_factor;
            let delta = reference - rotation * query - residual_centroid;
            delta.iter().fold(norm, |norm, value| norm.hypot(*value))
        })
}

fn identity_rotational_residual(
    reference: &[[f64; 3]],
    reference_factor: f64,
    query: &[[f64; 3]],
    query_centroid: [f64; 3],
    query_factor: f64,
    residual_centroid: na::Vector3<f64>,
) -> Option<(f64, f64)> {
    let query_centroid = na::Vector3::from(query_centroid) * query_factor;
    let (inertia, torque) = reference.iter().zip(query).fold(
        (na::Matrix3::zeros(), na::Vector3::zeros()),
        |(inertia, torque), (reference, query)| {
            let reference = na::Vector3::from(*reference) * reference_factor;
            let query = na::Vector3::from(*query) * query_factor;
            let centered_query = query - query_centroid;
            let residual = reference - query - residual_centroid;
            (
                inertia + na::Matrix3::identity() * centered_query.norm_squared()
                    - centered_query * centered_query.transpose(),
                torque + centered_query.cross(&residual),
            )
        },
    );
    let rotation = inertia.lu().solve(&torque)?;
    let residual_norm = query.iter().fold(0.0_f64, |norm, query| {
        let query = na::Vector3::from(*query) * query_factor - query_centroid;
        rotation
            .cross(&query)
            .iter()
            .fold(norm, |norm, value| norm.hypot(*value))
    });
    Some((rotation.norm(), residual_norm))
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
        let query = reference.map(|[x, y, z]| [-y + 10.0, x - 4.0, z + 2.0]);
        let forward = kabsch_rmsd(&reference, &query).unwrap();
        let reverse = kabsch_rmsd(&query, &reference).unwrap();
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
        let mut query = reference.map(|[x, y, z]| [-y + 10.0, x - 4.0, z + 2.0]);
        query[3][0] += 0.125;
        let reference = prepare_coordinates(&reference).unwrap();
        let query = prepare_coordinates(&query).unwrap();
        assert_eq!(
            kabsch_prepared_selected_rmsd(&reference, &query)
                .unwrap()
                .to_bits(),
            kabsch_prepared_rmsd(&reference, &query).unwrap().to_bits()
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
        let query = [
            [10.0, -4.0, 2.0],
            [10.0, -3.0, 2.0],
            [9.0, -4.0, 2.0],
            [10.0, -4.0, 5.0],
        ];
        let reference = prepare_coordinate_union(&reference, 3, 3).unwrap();
        let query = prepare_coordinate_union(&query, 3, 3).unwrap();
        let rmsd = kabsch_prepared_selected_rmsd(&reference, &query).unwrap();
        let reverse = kabsch_prepared_selected_rmsd(&query, &reference).unwrap();
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
        let query = reference.map(|[x, y, z]| [-y + 10.0, x - 4.0, z + 2.0]);
        let reference = prepare_coordinate_union(&reference, 4, 4).unwrap();
        let query = prepare_coordinate_union(&query, 4, 4).unwrap();
        assert!(matches!(
            kabsch_prepared_selected_rmsd(&reference, &query),
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
        let query = reference.map(|[x, y, z]| [x - angle * y + 10.0, angle * x + y - 4.0, z + 2.0]);
        let reference = prepare_coordinate_union(&reference, 4, 4).unwrap();
        let query = prepare_coordinate_union(&query, 4, 4).unwrap();
        let transform = fit_prepared_transform(&reference, &query).unwrap();
        assert_eq!(transform.rotation, na::Matrix3::identity());
        assert!(matches!(
            kabsch_prepared_selected_rmsd(&reference, &query),
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
        let mut query = reference.map(|[x, y, z]| [-y + 10.0, x - 4.0, z + 2.0]);
        query[3][0] += 0.01;
        let reference = prepare_coordinate_union(&reference, 4, 4).unwrap();
        let query = prepare_coordinate_union(&query, 4, 4).unwrap();
        let transform = fit_prepared_transform(&reference, &query).unwrap();
        let solver_tolerance = f64::EPSILON * 4.0_f64.sqrt() * 64.0;
        assert!(transform.fit_residual_norm > solver_tolerance);
        assert!(matches!(
            kabsch_prepared_selected_rmsd(&reference, &query),
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
        let query = reference.map(|[x, y, z]| [-y + 7.0, x - 11.0, z + 5.0]);
        assert!(kabsch_rmsd(&reference, &query).unwrap() < 1e-9);
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
            let query =
                reference.map(|[x, y, z]| [-y + translation, x + translation, z + translation]);
            assert!(kabsch_rmsd(&reference, &query).unwrap() < 1e-12);
        }
    }

    #[test]
    fn kabsch_keeps_extreme_finite_coordinates_finite() {
        for scale in [1e-200, 1e200] {
            let reference = [[0.0, 0.0, 0.0], [scale, 0.0, 0.0], [0.0, scale, 0.0]];
            let query = reference.map(|[x, y, z]| [-y, x, z]);
            let rmsd = kabsch_rmsd(&reference, &query).unwrap();
            assert!(rmsd.is_finite());
            assert!(rmsd / scale < 1e-12);
        }
    }

    #[test]
    fn rmsd_rescaling_preserves_subnormal_results() {
        let minimum_subnormal = f64::from_bits(1);
        assert_eq!(scaled_rmsd(2.0, minimum_subnormal, 4), minimum_subnormal);
        assert_eq!(
            scaled_rmsd(minimum_subnormal, 1e308, 4),
            minimum_subnormal * 1e308 / 2.0
        );
        assert_eq!(
            scaled_rmsd(2.0 * minimum_subnormal, 1e308, 9),
            2.0 * minimum_subnormal * 1e308 / 3.0
        );
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
        let mut query = reference;
        query[1][1] = 1.0;
        query[1][2] = 1.0;
        query[2][0] = 1.0;
        query[2][2] = 1.0;
        query[3][0] = 1.0;
        query[3][1] = 1.0;
        let expected = (3.0_f64 / 4.0).sqrt();
        let forward = kabsch_rmsd(&reference, &query).unwrap();
        let reverse = kabsch_rmsd(&query, &reference).unwrap();
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
        for query in [rotated, infinitesimally_rotated, mixed_strain_and_rotation] {
            assert!(matches!(
                kabsch_rmsd(&reference, &query),
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
}
