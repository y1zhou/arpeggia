//! K-medoids clustering over pairwise RMSD matrices.

use crate::pairwise_rmsd::PairwiseRmsdMatrix;
use crate::utils::polars_calculation_error;
use crate::{Analysis, AnalysisWarning, ArpeggiaError, ArpeggiaResult, WarningCode};
use clap::ValueEnum;
use kmedoids::ArrayAdapter;
use polars::prelude::*;

const IDENTICAL_RMSD_CUTOFF_ANGSTROM: f64 = 1e-12;

/// Clustering algorithms available for pairwise RMSD matrices.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, ValueEnum)]
pub enum ClusteringMethod {
    /// Partition structures around observed representative structures.
    #[default]
    KMedoids,
}

/// Options controlling clustering of a pairwise RMSD matrix.
#[derive(Clone, Debug)]
pub struct ClusterOptions {
    /// Clustering algorithm.
    pub method: ClusteringMethod,
    /// Exact number of clusters; takes priority over `max_clusters`.
    pub num_clusters: Option<usize>,
    /// Largest cluster count considered by automatic selection.
    pub max_clusters: Option<usize>,
    /// Maximum optimization iterations.
    pub max_iterations: usize,
}

impl Default for ClusterOptions {
    fn default() -> Self {
        Self {
            method: ClusteringMethod::KMedoids,
            num_clusters: None,
            max_clusters: None,
            max_iterations: 100,
        }
    }
}

impl ClusterOptions {
    /// Validate clustering options that do not depend on the input size.
    pub fn validate_without_structure_count(&self) -> ArpeggiaResult<()> {
        if self.max_iterations == 0 {
            return Err(ArpeggiaError::InvalidArgument(
                "max_iterations must be positive".into(),
            ));
        }
        if self.num_clusters == Some(0) {
            return Err(ArpeggiaError::InvalidArgument(
                "num_clusters must be positive".into(),
            ));
        }
        if self.num_clusters.is_none() {
            match self.max_clusters {
                None => {
                    return Err(ArpeggiaError::InvalidArgument(
                        "one of num_clusters or max_clusters is required".into(),
                    ));
                }
                Some(0 | 1) => {
                    return Err(ArpeggiaError::InvalidArgument(
                        "max_clusters must be at least 2".into(),
                    ));
                }
                Some(_) => {}
            }
        }
        Ok(())
    }

    /// Validate clustering bounds for a known structure count.
    pub fn validate(&self, structure_count: usize) -> ArpeggiaResult<()> {
        self.validate_without_structure_count()?;
        if structure_count < 3 {
            return Err(ArpeggiaError::InvalidArgument(
                "structure clustering requires at least three structures".into(),
            ));
        }
        if structure_count > u32::MAX as usize {
            return Err(ArpeggiaError::InvalidArgument(
                "structure count exceeds UInt32 cluster output capacity".into(),
            ));
        }
        if let Some(k) = self.num_clusters {
            if (1..=structure_count).contains(&k) {
                Ok(())
            } else {
                Err(ArpeggiaError::InvalidArgument(format!(
                    "num_clusters must be between 1 and {structure_count}"
                )))
            }
        } else if let Some(k) = self.max_clusters {
            if !(2..structure_count).contains(&k) {
                Err(ArpeggiaError::InvalidArgument(format!(
                    "max_clusters must be between 2 and {}",
                    structure_count - 1
                )))
            } else {
                Ok(())
            }
        } else {
            unreachable!("cluster bounds were validated before count-dependent checks")
        }
    }
}

struct ScaledPairwiseRmsdMatrix<'a> {
    matrix: &'a PairwiseRmsdMatrix,
    reciprocal: Option<f64>,
}

impl<'a> ScaledPairwiseRmsdMatrix<'a> {
    fn new(
        matrix: &'a PairwiseRmsdMatrix,
        maximum: f64,
        minimum_positive: Option<f64>,
    ) -> ArpeggiaResult<Self> {
        let safe_distance = f64::MAX / (4.0 * matrix.len() as f64);
        let scale = (maximum / safe_distance).max(1.0);
        if maximum > 0.0 && minimum_positive.is_some_and(|minimum| minimum / maximum == 0.0) {
            return Err(ArpeggiaError::Calculation(
                "pairwise RMSD dynamic range cannot be represented safely during clustering".into(),
            ));
        }
        Ok(Self {
            matrix,
            reciprocal: (scale > 1.0).then(|| scale.recip()),
        })
    }
}

impl ArrayAdapter<f64> for ScaledPairwiseRmsdMatrix<'_> {
    fn len(&self) -> usize {
        self.matrix.len()
    }

    // Required by kmedoids::ArrayAdapter; construction already guarantees a
    // complete packed symmetric matrix.
    fn is_square(&self) -> bool {
        true
    }

    fn get(&self, left: usize, right: usize) -> f64 {
        let distance = self.matrix.distance(left, right);
        self.reciprocal.map_or(distance, |scale| distance * scale)
    }
}

/// Cluster a validated packed RMSD matrix.
pub fn cluster_pairwise_rmsd(
    matrix: &PairwiseRmsdMatrix,
    options: &ClusterOptions,
) -> ArpeggiaResult<Analysis<DataFrame>> {
    let n = matrix.len();
    options.validate(n)?;

    let mut warnings = Vec::new();
    let (medoids, assignments) = if let Some(k) = options.num_clusters {
        if options.max_clusters.is_some() {
            warnings.push(AnalysisWarning::new(
                WarningCode::ArgumentIgnored,
                "num_clusters takes priority; max_clusters was ignored",
            ));
        }
        if k == n {
            let identity = (0..n).collect::<Vec<_>>();
            (identity.clone(), identity)
        } else {
            let (maximum, minimum_positive) = matrix.distance_range();
            run_fasterpam(
                &ScaledPairwiseRmsdMatrix::new(matrix, maximum, minimum_positive)?,
                k,
                options.max_iterations,
            )?
        }
    } else {
        let max_k = options
            .max_clusters
            .expect("cluster options were validated before use");
        let (maximum, minimum_positive) = matrix.distance_range();
        if maximum <= IDENTICAL_RMSD_CUTOFF_ANGSTROM {
            (vec![0], vec![0; n])
        } else {
            run_dynmsc(
                &ScaledPairwiseRmsdMatrix::new(matrix, maximum, minimum_positive)?,
                max_k,
                options.max_iterations,
            )?
        }
    };
    let value = cluster_dataframe(matrix, &medoids, &assignments)?;
    Ok(Analysis::new(value, warnings))
}

fn run_fasterpam(
    matrix: &ScaledPairwiseRmsdMatrix<'_>,
    k: usize,
    max_iterations: usize,
) -> ArpeggiaResult<(Vec<usize>, Vec<usize>)> {
    let (_, assignments, mut medoids) = kmedoids::pam_build::<_, f64, f64>(matrix, k);
    if k == 1 {
        return Ok((medoids, assignments));
    }
    if medoids.len() < k {
        let mut selected = vec![false; matrix.len()];
        for &medoid in &medoids {
            selected[medoid] = true;
        }
        for (index, selected) in selected.into_iter().enumerate() {
            if !selected {
                medoids.push(index);
                if medoids.len() == k {
                    break;
                }
            }
        }
    }
    let mut remaining_iterations = max_iterations;
    loop {
        let requested_iterations = remaining_iterations.max(1);
        let (_, mut assignments, iterations, swaps) =
            kmedoids::fasterpam::<_, f64, f64>(matrix, &mut medoids, requested_iterations);
        remaining_iterations = remaining_iterations.saturating_sub(iterations);
        if iterations == requested_iterations && swaps > 0 && remaining_iterations == 0 {
            let (_, diagnostic_assignments, _, diagnostic_swaps) =
                kmedoids::fasterpam::<_, f64, f64>(matrix, &mut medoids, 1);
            if diagnostic_swaps > 0 {
                return Err(ArpeggiaError::Calculation(format!(
                    "k-medoids reached the {max_iterations}-iteration limit"
                )));
            }
            assignments = diagnostic_assignments;
        }

        let mut previous_medoids = medoids.clone();
        previous_medoids.sort_unstable();
        let (canonical_medoids, canonical_assignments) =
            canonicalize_medoids(matrix, medoids, assignments);
        medoids = canonical_medoids;
        if medoids == previous_medoids {
            return Ok((medoids, canonical_assignments));
        }
        if remaining_iterations == 0 {
            let mut diagnostic_medoids = medoids.clone();
            let (_, _, _, diagnostic_swaps) =
                kmedoids::fasterpam::<_, f64, f64>(matrix, &mut diagnostic_medoids, 1);
            if diagnostic_swaps > 0 {
                return Err(ArpeggiaError::Calculation(format!(
                    "k-medoids reached the {max_iterations}-iteration limit"
                )));
            }
            return Ok((medoids, canonical_assignments));
        }
    }
}

fn canonicalize_medoids(
    matrix: &ScaledPairwiseRmsdMatrix<'_>,
    mut medoids: Vec<usize>,
    mut assignments: Vec<usize>,
) -> (Vec<usize>, Vec<usize>) {
    loop {
        let previous_memberships = assigned_medoid_ids(&medoids, &assignments);
        let (next_medoids, next_assignments) =
            canonicalize_medoids_once(matrix, medoids, &assignments);
        let next_memberships = assigned_medoid_ids(&next_medoids, &next_assignments);
        medoids = next_medoids;
        assignments = next_assignments;
        if previous_memberships == next_memberships {
            return (medoids, assignments);
        }
    }
}

fn canonicalize_medoids_once(
    matrix: &ScaledPairwiseRmsdMatrix<'_>,
    mut medoids: Vec<usize>,
    assignments: &[usize],
) -> (Vec<usize>, Vec<usize>) {
    for slot in 0..medoids.len() {
        let members = assignments
            .iter()
            .enumerate()
            .filter_map(|(index, assigned)| (*assigned == slot).then_some(index))
            .collect::<Vec<_>>();
        if members.is_empty() {
            continue;
        }
        let mut best = medoids[slot];
        let best_loss = members
            .iter()
            .map(|member| matrix.get(*member, best))
            .sum::<f64>();
        for &candidate in &members {
            if medoids
                .iter()
                .enumerate()
                .any(|(other_slot, medoid)| other_slot != slot && *medoid == candidate)
            {
                continue;
            }
            let loss = members
                .iter()
                .map(|member| matrix.get(*member, candidate))
                .sum::<f64>();
            if loss == best_loss && candidate < best {
                best = candidate;
            }
        }
        medoids[slot] = best;
    }

    medoids.sort_unstable();
    let assignments = (0..matrix.len())
        .map(|index| {
            if let Ok(slot) = medoids.binary_search(&index) {
                return slot;
            }
            medoids
                .iter()
                .enumerate()
                .min_by(|(left_slot, left), (right_slot, right)| {
                    matrix
                        .get(index, **left)
                        .total_cmp(&matrix.get(index, **right))
                        .then_with(|| left_slot.cmp(right_slot))
                })
                .map_or(0, |(slot, _)| slot)
        })
        .collect();
    (medoids, assignments)
}

fn assigned_medoid_ids(medoids: &[usize], assignments: &[usize]) -> Vec<usize> {
    assignments
        .iter()
        .map(|assignment| medoids[*assignment])
        .collect()
}

fn run_dynmsc(
    matrix: &ScaledPairwiseRmsdMatrix<'_>,
    max_k: usize,
    max_iterations: usize,
) -> ArpeggiaResult<(Vec<usize>, Vec<usize>)> {
    let (_, _, initial_medoids) = kmedoids::pam_build::<_, f64, f64>(matrix, max_k);
    let (_, assignments, iterations, _, medoids, _) =
        kmedoids::dynmsc::<_, f64, f64>(matrix, &initial_medoids, 2, max_iterations);
    let stage_limit = max_iterations.saturating_mul(initial_medoids.len() - 1);
    if iterations >= stage_limit {
        return Err(ArpeggiaError::Calculation(format!(
            "every automatic k-medoids stage reached the {max_iterations}-iteration limit"
        )));
    }
    Ok((medoids, assignments))
}

fn cluster_dataframe(
    matrix: &PairwiseRmsdMatrix,
    medoids: &[usize],
    assignments: &[usize],
) -> ArpeggiaResult<DataFrame> {
    if assignments.len() != matrix.len() || medoids.is_empty() {
        return Err(ArpeggiaError::Calculation(
            "k-medoids returned an incomplete assignment".into(),
        ));
    }
    let mut medoid_order = (0..medoids.len()).collect::<Vec<_>>();
    medoid_order.sort_unstable_by_key(|slot| medoids[*slot]);
    let mut cluster_ids = vec![0_u32; medoids.len()];
    for (cluster_id, slot) in medoid_order.into_iter().enumerate() {
        cluster_ids[slot] = u32::try_from(cluster_id)
            .map_err(|_| ArpeggiaError::Calculation("cluster identifier exceeds UInt32".into()))?;
    }

    let mut output_cluster_ids = Vec::with_capacity(matrix.len());
    let mut medoid_ids = Vec::with_capacity(matrix.len());
    let mut rmsd_to_medoid = Vec::with_capacity(matrix.len());
    for (index, &slot) in assignments.iter().enumerate() {
        let &medoid = medoids.get(slot).ok_or_else(|| {
            ArpeggiaError::Calculation("k-medoids returned an invalid medoid index".into())
        })?;
        output_cluster_ids.push(cluster_ids[slot]);
        medoid_ids.push(matrix.ids()[medoid].as_str());
        rmsd_to_medoid.push(matrix.distance(index, medoid));
    }
    df!(
        "id" => matrix.ids().iter().map(String::as_str).collect::<Vec<_>>(),
        "cluster_id" => output_cluster_ids,
        "medoid_id" => medoid_ids,
        "rmsd_to_medoid" => rmsd_to_medoid,
    )
    .map_err(polars_calculation_error)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    macro_rules! matrix {
        (ids: $ids:expr, $data:ident $(,)?) => {
            PairwiseRmsdMatrix::from_test_data($ids, $data)
        };
        (ids: $ids:expr, data: $data:expr $(,)?) => {
            PairwiseRmsdMatrix::from_test_data($ids, $data)
        };
    }

    fn one_dimensional_matrix(positions: &[f64]) -> PairwiseRmsdMatrix {
        let mut data = Vec::new();
        for row in 1..positions.len() {
            for column in 0..row {
                data.push((positions[row] - positions[column]).abs());
            }
        }
        matrix! {
            ids: (b'a'..)
                .take(positions.len())
                .map(char::from)
                .map(String::from)
                .collect(),
            data,
        }
    }

    #[test]
    fn fixed_k_medoids_separates_two_groups() {
        let matrix = matrix! {
            ids: ["a", "b", "c", "d", "e", "f"]
                .into_iter()
                .map(str::to_string)
                .collect(),
            data: vec![
                0.1, // b-a
                0.2, 0.1, // c-a, c-b
                5.0, 5.1, 5.0, // d-a..c
                5.1, 5.0, 5.1, 0.1, // e-a..d
                5.0, 5.1, 5.0, 0.2, 0.1, // f-a..e
            ],
        };
        let result = cluster_pairwise_rmsd(
            &matrix,
            &ClusterOptions {
                num_clusters: Some(2),
                ..ClusterOptions::default()
            },
        )
        .unwrap()
        .value;
        let clusters = result
            .column("cluster_id")
            .unwrap()
            .u32()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>();
        assert_eq!(&clusters[..3], &[0, 0, 0]);
        assert_eq!(&clusters[3..], &[1, 1, 1]);
        assert_eq!(result.column("id").unwrap().dtype(), &DataType::String);
        assert_eq!(
            result.column("rmsd_to_medoid").unwrap().dtype(),
            &DataType::Float64
        );
    }

    #[test]
    fn clustering_is_invariant_to_large_finite_distance_scales() {
        let matrix = matrix! {
            ids: ["a", "b", "c", "d", "e", "f"]
                .into_iter()
                .map(str::to_string)
                .collect(),
            data: vec![
                0.1, 0.2, 0.1, 5.0, 5.1, 5.0, 5.1, 5.0, 5.1, 0.1, 5.0, 5.1, 5.0, 0.2, 0.1,
            ],
        };
        let scaled = matrix! {
            ids: matrix.ids().to_vec(),
            data: matrix.packed_data().iter().map(|value| value * 1e307).collect(),
        };
        let options = ClusterOptions {
            num_clusters: Some(2),
            ..ClusterOptions::default()
        };
        let ordinary = cluster_pairwise_rmsd(&matrix, &options).unwrap().value;
        let extreme = cluster_pairwise_rmsd(&scaled, &options).unwrap().value;
        assert_eq!(
            ordinary.column("cluster_id").unwrap(),
            extreme.column("cluster_id").unwrap()
        );
        assert_eq!(
            ordinary.column("medoid_id").unwrap(),
            extreme.column("medoid_id").unwrap()
        );
    }

    #[test]
    fn clustering_rejects_unrepresentable_distance_dynamic_range() {
        let matrix = matrix! {
            ids: vec!["a".into(), "b".into(), "c".into(), "d".into()],
            data: vec![1e-16, 2e-16, 1e-16, 1e308, 1e308, 1e308],
        };
        assert!(matches!(
            cluster_pairwise_rmsd(
                &matrix,
                &ClusterOptions {
                    num_clusters: Some(2),
                    ..ClusterOptions::default()
                }
            ),
            Err(ArpeggiaError::Calculation(message))
                if message.contains("dynamic range")
        ));
    }

    #[test]
    fn automatic_k_medoids_short_circuits_identical_structures() {
        let matrix = matrix! {
            ids: vec!["a".into(), "b".into(), "c".into()],
            data: vec![IDENTICAL_RMSD_CUTOFF_ANGSTROM; 3],
        };
        let result = cluster_pairwise_rmsd(
            &matrix,
            &ClusterOptions {
                max_clusters: Some(2),
                ..ClusterOptions::default()
            },
        )
        .unwrap()
        .value;
        assert_eq!(
            result
                .column("cluster_id")
                .unwrap()
                .u32()
                .unwrap()
                .into_no_null_iter()
                .collect::<Vec<_>>(),
            [0, 0, 0]
        );
        let medoids = result.column("medoid_id").unwrap().str().unwrap();
        assert!((0..3).all(|index| medoids.get(index) == Some("a")));
    }

    #[test]
    fn automatic_k_medoids_selects_within_the_requested_bound() {
        let matrix = matrix! {
            ids: ["a", "b", "c", "d", "e", "f"]
                .into_iter()
                .map(str::to_string)
                .collect(),
            data: vec![
                0.1, 0.2, 0.1, 5.0, 5.1, 5.0, 5.1, 5.0, 5.1, 0.1, 5.0, 5.1, 5.0, 0.2, 0.1,
            ],
        };
        let result = cluster_pairwise_rmsd(
            &matrix,
            &ClusterOptions {
                max_clusters: Some(3),
                ..ClusterOptions::default()
            },
        )
        .unwrap()
        .value;
        let clusters = result
            .column("cluster_id")
            .unwrap()
            .u32()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>();
        assert_eq!(clusters, [0, 0, 0, 1, 1, 1]);
        let medoids = result
            .column("medoid_id")
            .unwrap()
            .str()
            .unwrap()
            .iter()
            .flatten()
            .collect::<Vec<_>>();
        assert_eq!(medoids, ["b", "b", "b", "e", "e", "e"]);
    }

    #[test]
    fn fixed_cluster_count_takes_priority_with_warning() {
        let matrix = matrix! {
            ids: vec!["a".into(), "b".into(), "c".into()],
            data: vec![1.0, 2.0, 1.0],
        };
        let result = cluster_pairwise_rmsd(
            &matrix,
            &ClusterOptions {
                num_clusters: Some(1),
                max_clusters: Some(3),
                ..ClusterOptions::default()
            },
        )
        .unwrap();
        assert_eq!(result.warnings[0].code, WarningCode::ArgumentIgnored);
    }

    #[test]
    fn fixed_k_medoids_supports_one_and_every_structure() {
        let matrix = matrix! {
            ids: vec!["a".into(), "b".into(), "c".into()],
            data: vec![1.0, 2.0, 1.0],
        };
        for k in [1, 3] {
            let result = cluster_pairwise_rmsd(
                &matrix,
                &ClusterOptions {
                    num_clusters: Some(k),
                    max_iterations: 1,
                    ..ClusterOptions::default()
                },
            )
            .unwrap()
            .value;
            let cluster_count = result
                .column("cluster_id")
                .unwrap()
                .u32()
                .unwrap()
                .max()
                .unwrap()
                + 1;
            assert_eq!(cluster_count as usize, k);
        }
    }

    #[test]
    fn fixed_k_medoids_keeps_requested_clusters_for_duplicate_observations() {
        let matrix = matrix! {
            ids: vec!["a".into(), "b".into(), "c".into()],
            data: vec![0.0; 3],
        };
        let result = cluster_pairwise_rmsd(
            &matrix,
            &ClusterOptions {
                num_clusters: Some(2),
                ..ClusterOptions::default()
            },
        )
        .unwrap()
        .value;
        assert_eq!(
            result
                .column("cluster_id")
                .unwrap()
                .u32()
                .unwrap()
                .into_no_null_iter()
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([0, 1])
        );
    }

    #[test]
    fn fixed_k_medoids_is_deterministic_for_ties() {
        let matrix = matrix! {
            ids: vec!["a".into(), "b".into(), "c".into(), "d".into()],
            data: vec![1.0; 6],
        };
        let options = ClusterOptions {
            num_clusters: Some(2),
            ..ClusterOptions::default()
        };
        let first = cluster_pairwise_rmsd(&matrix, &options).unwrap().value;
        let second = cluster_pairwise_rmsd(&matrix, &options).unwrap().value;
        assert_eq!(first, second);
    }

    #[test]
    fn fixed_k_medoids_uses_lowest_index_for_equal_loss() {
        let matrix = matrix! {
            ids: vec!["a".into(), "b".into(), "c".into(), "d".into()],
            data: vec![1.0, 8.0, 7.0, 8.0, 7.0, 0.0],
        };
        let result = cluster_pairwise_rmsd(
            &matrix,
            &ClusterOptions {
                num_clusters: Some(2),
                ..ClusterOptions::default()
            },
        )
        .unwrap()
        .value;
        let medoids = result
            .column("medoid_id")
            .unwrap()
            .str()
            .unwrap()
            .iter()
            .flatten()
            .collect::<Vec<_>>();
        assert_eq!(medoids, ["a", "a", "c", "c"]);
    }

    #[test]
    fn fixed_k_medoids_reoptimizes_after_tie_canonicalization() {
        let matrix = one_dimensional_matrix(&[1.0, 5.0, 7.0, 18.0, 23.0]);
        let result = cluster_pairwise_rmsd(
            &matrix,
            &ClusterOptions {
                num_clusters: Some(3),
                ..ClusterOptions::default()
            },
        )
        .unwrap()
        .value;
        let loss = result
            .column("rmsd_to_medoid")
            .unwrap()
            .f64()
            .unwrap()
            .sum()
            .unwrap();
        let medoids = result
            .column("medoid_id")
            .unwrap()
            .str()
            .unwrap()
            .iter()
            .flatten()
            .collect::<BTreeSet<_>>();
        assert_eq!(loss, 6.0);
        assert_eq!(medoids, BTreeSet::from(["b", "d", "e"]));
    }

    #[test]
    fn fixed_k_medoids_canonicalizes_after_assignment_ties() {
        let matrix = one_dimensional_matrix(&[0.0, 3.0, 1.0, 2.0]);
        let result = cluster_pairwise_rmsd(
            &matrix,
            &ClusterOptions {
                num_clusters: Some(2),
                ..ClusterOptions::default()
            },
        )
        .unwrap()
        .value;
        let medoids = result
            .column("medoid_id")
            .unwrap()
            .str()
            .unwrap()
            .iter()
            .flatten()
            .collect::<BTreeSet<_>>();
        assert_eq!(medoids, BTreeSet::from(["a", "b"]));
    }

    #[test]
    fn automatic_k_medoids_preserves_dynmsc_objective() {
        let matrix = one_dimensional_matrix(&[0.0, 1.0, 3.0, 10.0, 14.0]);
        let result = cluster_pairwise_rmsd(
            &matrix,
            &ClusterOptions {
                max_clusters: Some(4),
                ..ClusterOptions::default()
            },
        )
        .unwrap()
        .value;
        let medoids = result
            .column("medoid_id")
            .unwrap()
            .str()
            .unwrap()
            .iter()
            .flatten()
            .collect::<BTreeSet<_>>();
        assert_eq!(medoids, BTreeSet::from(["b", "c", "d", "e"]));
    }

    #[test]
    fn fixed_k_medoids_accepts_convergence_after_the_final_swap_pass() {
        let points = [
            (23.0_f64, 39.0_f64),
            (90.0, 92.0),
            (34.0, 25.0),
            (77.0, 30.0),
            (27.0, 74.0),
            (43.0, 82.0),
            (4.0, 8.0),
            (97.0, 68.0),
        ];
        let mut data = Vec::new();
        for row in 1..points.len() {
            for column in 0..row {
                data.push(
                    ((points[row].0 - points[column].0).powi(2)
                        + (points[row].1 - points[column].1).powi(2))
                    .sqrt(),
                );
            }
        }
        let matrix = matrix! {
            ids: (0..points.len()).map(|index| index.to_string()).collect(),
            data,
        };
        assert!(
            cluster_pairwise_rmsd(
                &matrix,
                &ClusterOptions {
                    num_clusters: Some(3),
                    max_iterations: 1,
                    ..ClusterOptions::default()
                }
            )
            .is_ok()
        );
    }

    #[test]
    fn fixed_k_medoids_reports_genuine_iteration_exhaustion() {
        let points = [
            (1782426941.0_f64, 854532251.0_f64),
            (231763727.0, 2865926783.0),
            (2835306018.0, 3293277652.0),
            (1944078868.0, 2447144979.0),
            (3486048204.0, 4007635998.0),
            (898367503.0, 3128196544.0),
            (3651832499.0, 4201332114.0),
            (4181355587.0, 3675685636.0),
            (3958615361.0, 2731153469.0),
            (3948766392.0, 45079807.0),
            (2277873536.0, 212778365.0),
            (1819579428.0, 1598085891.0),
        ];
        let mut data = Vec::new();
        for row in 1..points.len() {
            for column in 0..row {
                data.push(f64::hypot(
                    points[row].0 - points[column].0,
                    points[row].1 - points[column].1,
                ));
            }
        }
        let matrix = matrix! {
            ids: (0..points.len()).map(|index| index.to_string()).collect(),
            data,
        };
        assert!(matches!(
            cluster_pairwise_rmsd(
                &matrix,
                &ClusterOptions {
                    num_clusters: Some(4),
                    max_iterations: 1,
                    ..ClusterOptions::default()
                },
            ),
            Err(ArpeggiaError::Calculation(_))
        ));
    }

    #[test]
    fn automatic_k_medoids_counts_only_actual_duplicate_stages() {
        let positions = [0.0_f64, 0.0, 10.0, 10.0, 20.0, 20.0];
        let mut data = Vec::new();
        for row in 1..positions.len() {
            for column in 0..row {
                data.push((positions[row] - positions[column]).abs());
            }
        }
        let matrix = matrix! {
            ids: (0..positions.len())
                .map(|index| index.to_string())
                .collect(),
            data,
        };
        assert!(matches!(
            cluster_pairwise_rmsd(
                &matrix,
                &ClusterOptions {
                    max_clusters: Some(5),
                    max_iterations: 1,
                    ..ClusterOptions::default()
                }
            ),
            Err(ArpeggiaError::Calculation(_))
        ));
    }

    #[test]
    fn cluster_options_require_a_valid_bound() {
        assert!(ClusterOptions::default().validate(3).is_err());
        assert!(
            ClusterOptions {
                max_clusters: Some(3),
                ..ClusterOptions::default()
            }
            .validate(3)
            .is_err()
        );
        assert!(
            ClusterOptions {
                num_clusters: Some(2),
                ..ClusterOptions::default()
            }
            .validate(3)
            .is_ok()
        );
    }
}
