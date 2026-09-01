//! Ensemble structure input, pairwise RMSD, and clustering.

use crate::rmsd::{
    AtomIdentity, PreparedCoordinates, ResidueSelector, kabsch_prepared_rmsd, prepare_coordinates,
    select_coordinates, validate_correspondence,
};
use crate::utils::polars_calculation_error;
use crate::{
    Analysis, AnalysisWarning, ArpeggiaError, ArpeggiaResult, AtomSubset, WarningCode, load_model,
};
use clap::ValueEnum;
use kmedoids::ArrayAdapter;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{Cursor, Read, Seek};
use std::path::{Path, PathBuf};

const FALLBACK_WARNING_BYTES: u64 = 8 * 1024 * 1024 * 1024;
const IDENTICAL_RMSD_CUTOFF_ANGSTROM: f64 = 1e-12;

/// One named structure participating in an ensemble calculation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StructureObservation {
    /// Case-sensitive identifier used in result tables.
    id: String,
    /// Canonical path to the structure file.
    path: PathBuf,
}

impl StructureObservation {
    /// Case-sensitive identifier used in result tables.
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Canonical structure path.
    pub fn path(&self) -> &Path {
        &self.path
    }
}

/// Options controlling pairwise structural superposition.
#[derive(Clone, Debug, Default)]
pub struct PairwiseRmsdOptions {
    /// Coordinate model number, with zero selecting the first model.
    pub model_num: usize,
    /// Chain and author-residue selection grammar.
    pub residues: String,
    /// Atom population used for the fit and RMSD.
    pub atoms: AtomSubset,
    /// Worker count, with zero selecting all available processors.
    pub num_threads: usize,
    /// Skip heuristic RAM warnings and failures.
    pub bypass_mem_check: bool,
}

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

/// Read structure observations from a directory or tabular manifest.
pub fn read_structure_observations(
    input: &Path,
    id_column: &str,
    path_column: &str,
) -> ArpeggiaResult<Vec<StructureObservation>> {
    let input = input.canonicalize().map_err(|error| {
        ArpeggiaError::Io(std::io::Error::new(
            error.kind(),
            format!("cannot resolve input {}: {error}", input.display()),
        ))
    })?;
    let observations = if input.is_dir() {
        read_structure_directory(&input)?
    } else {
        read_structure_manifest(&input, id_column, path_column)?
    };
    validate_observations(observations)
}

/// Calculate every unordered RMSD pair and return a long Polars table.
pub fn get_pairwise_rmsd(
    input: &Path,
    id_column: &str,
    path_column: &str,
    options: &PairwiseRmsdOptions,
) -> ArpeggiaResult<Analysis<DataFrame>> {
    let observations = read_structure_observations(input, id_column, path_column)?;
    get_pairwise_rmsd_matrix(&observations, options).and_then(|analysis| {
        analysis
            .value
            .into_dataframe()
            .map(|value| Analysis::new(value, analysis.warnings))
    })
}

#[derive(Clone, Debug)]
/// Symmetric pairwise RMSD matrix stored as a packed lower triangle.
pub struct PairwiseRmsdMatrix {
    ids: Vec<String>,
    data: Vec<f64>,
}

impl PairwiseRmsdMatrix {
    /// Case-sensitive structure IDs in matrix order.
    pub fn ids(&self) -> &[String] {
        &self.ids
    }

    /// Number of structures represented by the matrix.
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// Whether the matrix contains no structures.
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    /// RMSD between two matrix indices.
    pub fn get(&self, left: usize, right: usize) -> Option<f64> {
        (left < self.len() && right < self.len()).then(|| self.distance(left, right))
    }

    fn distance(&self, left: usize, right: usize) -> f64 {
        if left == right {
            0.0
        } else {
            let (row, column) = (left.max(right), left.min(right));
            self.data[row * (row - 1) / 2 + column]
        }
    }

    fn distance_range(&self) -> (f64, Option<f64>) {
        self.data
            .iter()
            .copied()
            .fold((0.0, None), |(maximum, minimum_positive), distance| {
                (
                    maximum.max(distance),
                    if distance > 0.0 {
                        Some(minimum_positive.map_or(distance, |minimum| minimum.min(distance)))
                    } else {
                        minimum_positive
                    },
                )
            })
    }

    /// Convert the packed matrix to one row per unordered pair.
    pub fn to_dataframe(&self) -> ArpeggiaResult<DataFrame> {
        pairwise_dataframe(&self.ids, self.data.clone())
    }

    fn into_dataframe(self) -> ArpeggiaResult<DataFrame> {
        pairwise_dataframe(&self.ids, self.data)
    }

    /// Validate a complete long pairwise table and pack it for clustering.
    pub fn from_dataframe(
        dataframe: &DataFrame,
        expected_ids: Option<&[String]>,
        minimum_ids: usize,
    ) -> ArpeggiaResult<Self> {
        let left = dataframe
            .column("id_1")
            .map_err(|_| missing_table_column("pairwise", "id_1"))?
            .str()
            .map_err(|_| invalid_table_type("pairwise", "id_1", "String"))?;
        let right = dataframe
            .column("id_2")
            .map_err(|_| missing_table_column("pairwise", "id_2"))?
            .str()
            .map_err(|_| invalid_table_type("pairwise", "id_2", "String"))?;
        let rmsd = dataframe
            .column("rmsd")
            .map_err(|_| missing_table_column("pairwise", "rmsd"))?;
        if !rmsd.dtype().is_numeric() {
            return Err(ArpeggiaError::InvalidArgument(
                "pairwise rmsd must be numeric".into(),
            ));
        }
        if let Some(expected) = expected_ids {
            let pair_count = checked_pair_count(expected.len())?;
            if dataframe.height() != pair_count {
                return Err(ArpeggiaError::InvalidArgument(format!(
                    "pairwise table has {} rows; {} expected IDs require exactly {pair_count}",
                    dataframe.height(),
                    expected.len()
                )));
            }
        }
        if left.null_count() > 0 {
            return Err(ArpeggiaError::InvalidArgument(
                "pairwise id_1 contains null values".into(),
            ));
        }
        if right.null_count() > 0 {
            return Err(ArpeggiaError::InvalidArgument(
                "pairwise id_2 contains null values".into(),
            ));
        }
        if left.iter().chain(right.iter()).flatten().any(str::is_empty) {
            return Err(ArpeggiaError::InvalidArgument(
                "pairwise IDs cannot be empty".into(),
            ));
        }
        let mut id_values = left.clone();
        id_values.append(right).map_err(polars_input_error)?;
        let mut ids = id_values
            .unique()
            .map_err(polars_input_error)?
            .iter()
            .flatten()
            .map(str::to_string)
            .collect::<Vec<_>>();
        ids.sort_unstable();
        if ids.len() < minimum_ids {
            return Err(ArpeggiaError::InvalidArgument(format!(
                "pairwise table contains {} distinct IDs but at least {minimum_ids} are required",
                ids.len()
            )));
        }
        if let Some(expected) = expected_ids
            && ids != expected
        {
            return Err(ArpeggiaError::InvalidArgument(
                "pairwise cache IDs do not match current structure IDs; remove the cache file to recalculate"
                    .into(),
            ));
        }

        let id_index = ids
            .iter()
            .enumerate()
            .map(|(index, id)| (id.as_str(), index))
            .collect::<HashMap<_, _>>();
        let pair_count = checked_pair_count(ids.len())?;
        if dataframe.height() != pair_count {
            return Err(ArpeggiaError::InvalidArgument(format!(
                "pairwise table has {} rows; {} IDs require exactly {pair_count}",
                dataframe.height(),
                ids.len()
            )));
        }
        if rmsd.null_count() > 0 {
            return Err(ArpeggiaError::InvalidArgument(
                "pairwise rmsd contains null values".into(),
            ));
        }
        let rmsd = rmsd.cast(&DataType::Float64).map_err(polars_input_error)?;
        let rmsd = rmsd.f64().map_err(polars_input_error)?;
        let mut data = vec![f64::NAN; pair_count];
        for (input_row, ((left_id, right_id), value)) in left
            .iter()
            .flatten()
            .zip(right.iter().flatten())
            .zip(rmsd.into_no_null_iter())
            .enumerate()
        {
            if !value.is_finite() || value < 0.0 {
                return Err(ArpeggiaError::InvalidArgument(format!(
                    "pairwise rmsd at row {input_row} must be finite and non-negative"
                )));
            }
            let left_index = id_index[left_id];
            let right_index = id_index[right_id];
            if left_index == right_index {
                return Err(ArpeggiaError::InvalidArgument(format!(
                    "pairwise table contains diagonal row for {left_id}"
                )));
            }
            let (row, column) = (left_index.max(right_index), left_index.min(right_index));
            let output_index = row * (row - 1) / 2 + column;
            if !data[output_index].is_nan() {
                return Err(ArpeggiaError::InvalidArgument(format!(
                    "pairwise table contains duplicate pair {left_id}/{right_id}"
                )));
            }
            data[output_index] = value;
        }
        Ok(Self { ids, data })
    }
}

fn pairwise_dataframe(ids: &[String], rmsd: Vec<f64>) -> ArpeggiaResult<DataFrame> {
    let pair_count = rmsd.len();
    let mut id_1 = StringChunkedBuilder::new("id_1".into(), pair_count);
    let mut id_2 = StringChunkedBuilder::new("id_2".into(), pair_count);
    for row in 1..ids.len() {
        for column in 0..row {
            id_1.append_value(&ids[column]);
            id_2.append_value(&ids[row]);
        }
    }
    DataFrame::new(
        pair_count,
        vec![
            id_1.finish().into_column(),
            id_2.finish().into_column(),
            Column::new("rmsd".into(), rmsd),
        ],
    )
    .map_err(polars_calculation_error)
}

impl ArrayAdapter<f64> for PairwiseRmsdMatrix {
    fn len(&self) -> usize {
        self.ids.len()
    }

    fn is_square(&self) -> bool {
        checked_pair_count(self.ids.len()).is_ok_and(|count| count == self.data.len())
    }

    fn get(&self, left: usize, right: usize) -> f64 {
        self.distance(left, right)
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

    fn is_square(&self) -> bool {
        self.matrix.is_square()
    }

    fn get(&self, left: usize, right: usize) -> f64 {
        let distance = self.matrix.distance(left, right);
        self.reciprocal.map_or(distance, |scale| distance * scale)
    }
}

/// Read and validate a cached pairwise RMSD table for the expected IDs.
///
/// The packed-matrix memory guard runs before the table is decoded unless
/// `bypass_mem_check` is true.
pub fn read_pairwise_matrix(
    path: &Path,
    expected_ids: &[String],
    bypass_mem_check: bool,
) -> ArpeggiaResult<Analysis<PairwiseRmsdMatrix>> {
    let expected_rows = checked_pair_count(expected_ids.len())?;
    let warnings = check_packed_matrix_memory(expected_rows, bypass_mem_check)?;
    let dataframe = read_dataframe(path, &["id_1", "id_2", "rmsd"], Some(expected_rows))?;
    let matrix = PairwiseRmsdMatrix::from_dataframe(&dataframe, Some(expected_ids), 2)?;
    Ok(Analysis::new(matrix, warnings))
}

/// Calculate every unordered RMSD pair into a packed matrix.
pub fn get_pairwise_rmsd_matrix(
    observations: &[StructureObservation],
    options: &PairwiseRmsdOptions,
) -> ArpeggiaResult<Analysis<PairwiseRmsdMatrix>> {
    let observations = validate_observations(observations.to_vec())?;
    if observations.len() < 2 {
        return Err(ArpeggiaError::InvalidArgument(
            "pairwise RMSD requires at least two structures".into(),
        ));
    }
    let pair_count = checked_pair_count(observations.len())?;
    let mut warnings = check_packed_matrix_memory(pair_count, options.bypass_mem_check)?;
    let matrix_bytes = pair_count.checked_mul(size_of::<f64>()).ok_or_else(|| {
        ArpeggiaError::InvalidArgument("pairwise matrix size overflows usize".into())
    })? as u64;

    let selector = ResidueSelector::parse(&options.residues)?;
    let first_loaded = load_observation(&observations[0], options, &selector, None)?;
    let reference_keys = first_loaded
        .keys
        .expect("the first Structure Observation retains its atom keys");
    let atom_count = first_loaded.coordinates.len();
    let coordinate_bytes = observations
        .len()
        .checked_mul(atom_count)
        .and_then(|value| value.checked_mul(3 * size_of::<f64>()))
        .ok_or_else(|| {
            ArpeggiaError::InvalidArgument("selected-coordinate size overflows usize".into())
        })? as u64;
    let combined_bytes = matrix_bytes.checked_add(coordinate_bytes).ok_or_else(|| {
        ArpeggiaError::InvalidArgument("combined memory estimate overflows u64".into())
    })?;
    warnings.extend(check_memory(
        combined_bytes,
        options.bypass_mem_check,
        "packed RMSD matrix plus selected coordinates",
    )?);
    warnings.extend(first_loaded.warnings);

    let parse_threads = effective_thread_count(options.num_threads).min(8);
    let remaining = {
        let parse_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(parse_threads)
            .build()
            .map_err(|error| ArpeggiaError::Calculation(error.to_string()))?;
        parse_pool.install(|| {
            observations[1..]
                .par_iter()
                .map(|observation| {
                    load_observation(observation, options, &selector, Some(&reference_keys))
                })
                .collect::<ArpeggiaResult<Vec<_>>>()
        })
    }?;

    let mut coordinates = Vec::with_capacity(observations.len());
    coordinates.push(first_loaded.coordinates);
    for prepared in remaining {
        warnings.extend(prepared.warnings);
        coordinates.push(prepared.coordinates);
    }

    let pair_threads = effective_thread_count(options.num_threads).min(pair_count);
    let chunk_size = pair_count.div_ceil(pair_threads.saturating_mul(8));
    let mut data = vec![0.0; pair_count];
    let pair_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(pair_threads)
        .build()
        .map_err(|error| ArpeggiaError::Calculation(error.to_string()))?;
    pair_pool.install(|| {
        data.par_chunks_mut(chunk_size)
            .enumerate()
            .try_for_each(|(chunk_index, output)| {
                let start = chunk_index * chunk_size;
                let (mut row, mut column) = pair_at_index(start, observations.len());
                for value in output {
                    *value = kabsch_prepared_rmsd(&coordinates[column], &coordinates[row])?;
                    column += 1;
                    if column == row {
                        row += 1;
                        column = 0;
                    }
                }
                ArpeggiaResult::Ok(())
            })
    })?;

    Ok(Analysis::new(
        PairwiseRmsdMatrix {
            ids: observations
                .iter()
                .map(|observation| observation.id.clone())
                .collect(),
            data,
        },
        warnings,
    ))
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
        medoid_ids.push(matrix.ids[medoid].as_str());
        rmsd_to_medoid.push(matrix.distance(index, medoid));
    }
    df!(
        "id" => matrix.ids.iter().map(String::as_str).collect::<Vec<_>>(),
        "cluster_id" => output_cluster_ids,
        "medoid_id" => medoid_ids,
        "rmsd_to_medoid" => rmsd_to_medoid,
    )
    .map_err(polars_calculation_error)
}

struct PreparedObservation {
    keys: Option<Vec<AtomIdentity>>,
    coordinates: PreparedCoordinates,
    warnings: Vec<AnalysisWarning>,
}

fn load_observation(
    observation: &StructureObservation,
    options: &PairwiseRmsdOptions,
    selector: &ResidueSelector,
    reference_keys: Option<&[AtomIdentity]>,
) -> ArpeggiaResult<PreparedObservation> {
    let path = observation.path.to_str().ok_or_else(|| {
        ArpeggiaError::InvalidArgument(format!(
            "structure path is not valid UTF-8: {}",
            observation.path.display()
        ))
    })?;
    let loaded = load_model(path)?;
    let selected = select_coordinates(&loaded.value, options.model_num, selector, options.atoms)?;
    let keys = if let Some(reference_keys) = reference_keys {
        validate_correspondence(reference_keys, &selected.keys)?;
        None
    } else {
        Some(selected.keys)
    };
    let coordinates = prepare_coordinates(&selected.coordinates)?;
    let mut warnings = loaded.warnings;
    warnings.extend(selected.warnings);
    Ok(PreparedObservation {
        keys,
        coordinates,
        warnings,
    })
}

fn read_structure_directory(directory: &Path) -> ArpeggiaResult<Vec<StructureObservation>> {
    let mut observations = Vec::new();
    for entry in std::fs::read_dir(directory)? {
        let entry = entry?;
        let path = entry.path();
        if !path.is_file() || !is_structure_path(&path) {
            continue;
        }
        let id = path
            .file_stem()
            .and_then(|stem| stem.to_str())
            .ok_or_else(|| {
                ArpeggiaError::InvalidArgument(format!(
                    "structure filename is not valid UTF-8: {}",
                    path.display()
                ))
            })?
            .to_string();
        observations.push(StructureObservation {
            id,
            path: path.canonicalize()?,
        });
    }
    Ok(observations)
}

fn read_structure_manifest(
    manifest: &Path,
    id_column: &str,
    path_column: &str,
) -> ArpeggiaResult<Vec<StructureObservation>> {
    let dataframe = read_dataframe(manifest, &[id_column, path_column], None)?;
    let ids = dataframe
        .column(id_column)
        .map_err(|_| missing_table_column("manifest", id_column))?
        .str()
        .map_err(|_| invalid_table_type("manifest", id_column, "String"))?;
    let paths = dataframe
        .column(path_column)
        .map_err(|_| missing_table_column("manifest", path_column))?
        .str()
        .map_err(|_| invalid_table_type("manifest", path_column, "String"))?;
    for (row, (id, path)) in ids.iter().zip(paths.iter()).enumerate() {
        let id = id.ok_or_else(|| {
            ArpeggiaError::InvalidArgument(format!(
                "manifest column {id_column} contains null at row {row}"
            ))
        })?;
        let path = path.ok_or_else(|| {
            ArpeggiaError::InvalidArgument(format!(
                "manifest column {path_column} contains null at row {row}"
            ))
        })?;
        if id.is_empty() || path.is_empty() {
            return Err(ArpeggiaError::InvalidArgument(format!(
                "manifest row {row} contains an empty ID or path"
            )));
        }
    }
    let base = manifest.parent().unwrap_or_else(|| Path::new("."));
    ids.iter()
        .flatten()
        .zip(paths.iter().flatten())
        .map(|(id, path)| {
            let path = base.join(path);
            Ok(StructureObservation {
                id: id.to_string(),
                path: path.canonicalize().map_err(|error| {
                    ArpeggiaError::Io(std::io::Error::new(
                        error.kind(),
                        format!("cannot resolve manifest path {}: {error}", path.display()),
                    ))
                })?,
            })
        })
        .collect()
}

fn read_dataframe(
    path: &Path,
    columns: &[&str],
    expected_rows: Option<usize>,
) -> ArpeggiaResult<DataFrame> {
    if !path.is_file() {
        return Err(ArpeggiaError::InvalidArgument(format!(
            "tabular input is not a regular file: {}",
            path.display()
        )));
    }
    let extension = path
        .extension()
        .and_then(|extension| extension.to_str())
        .map(str::to_ascii_lowercase)
        .ok_or_else(|| {
            ArpeggiaError::InvalidArgument(format!(
                "tabular input has no supported extension: {}",
                path.display()
            ))
        })?;
    let mut file = File::open(path)?;
    let projection = || {
        columns
            .iter()
            .map(|column| PlSmallStr::from_str(column))
            .collect::<Vec<_>>()
    };
    match extension.as_str() {
        "csv" => CsvReadOptions::default()
            .with_n_rows(expected_rows.map(|rows| rows.saturating_add(1)))
            .with_columns(Some(projection().into()))
            .into_reader_with_file_handle(file)
            .finish()
            .map_err(polars_input_error),
        "parquet" => {
            let mut reader = ParquetReader::new(file);
            if let Some(expected) = expected_rows {
                let actual = reader.num_rows().map_err(polars_input_error)?;
                if actual != expected {
                    return Err(ArpeggiaError::InvalidArgument(format!(
                        "pairwise table has {actual} rows; expected exactly {expected}"
                    )));
                }
            }
            reader
                .with_columns(Some(
                    columns.iter().map(|column| column.to_string()).collect(),
                ))
                .finish()
                .map_err(polars_input_error)
        }
        "json" => {
            preflight_json_rows(&mut file, false, expected_rows)?;
            JsonReader::new(file)
                .with_json_format(JsonFormat::Json)
                .with_projection(Some(projection()))
                .finish()
                .map_err(polars_input_error)
        }
        "ndjson" => {
            if let Some(bytes) = preflight_json_rows(&mut file, true, expected_rows)? {
                JsonReader::new(Cursor::new(bytes))
                    .with_json_format(JsonFormat::JsonLines)
                    .with_projection(Some(projection()))
                    .finish()
                    .map_err(polars_input_error)
            } else {
                LazyJsonLineReader::new(PlRefPath::try_from_path(path).map_err(polars_input_error)?)
                    .finish()
                    .map_err(polars_input_error)?
                    .select(
                        columns
                            .iter()
                            .map(|column| col(*column))
                            .collect::<Vec<_>>(),
                    )
                    .collect()
                    .map_err(polars_input_error)
            }
        }
        _ => Err(ArpeggiaError::InvalidArgument(format!(
            "unsupported table format .{extension}; expected csv, parquet, json, or ndjson"
        ))),
    }
}

fn polars_input_error(error: PolarsError) -> ArpeggiaError {
    ArpeggiaError::InvalidArgument(error.to_string())
}

fn preflight_json_rows(
    file: &mut File,
    lines: bool,
    expected_rows: Option<usize>,
) -> ArpeggiaResult<Option<Vec<u8>>> {
    if let Some(expected) = expected_rows {
        let mut bytes = lines.then(Vec::new);
        let actual =
            bounded_json_row_count(file, lines, expected.saturating_add(1), bytes.as_mut())?;
        if actual != expected {
            return Err(ArpeggiaError::InvalidArgument(format!(
                "pairwise table has {actual} rows; expected exactly {expected}"
            )));
        }
        if !lines {
            file.rewind()?;
        }
        return Ok(bytes);
    }
    Ok(None)
}

fn bounded_json_row_count(
    file: &mut File,
    lines: bool,
    limit: usize,
    mut snapshot: Option<&mut Vec<u8>>,
) -> ArpeggiaResult<usize> {
    let mut buffer = [0_u8; 8192];
    let mut count = 0;
    let mut line_has_content = false;
    let mut array_started = false;
    let mut depth = 0_usize;
    let mut element_started = false;
    let mut in_string = false;
    let mut escaped = false;
    loop {
        let bytes_read = file.read(&mut buffer)?;
        if bytes_read == 0 {
            break;
        }
        if let Some(snapshot) = snapshot.as_deref_mut() {
            snapshot.extend_from_slice(&buffer[..bytes_read]);
        }
        for &byte in &buffer[..bytes_read] {
            if lines {
                if byte == b'\n' {
                    if line_has_content {
                        count += 1;
                        if count >= limit {
                            return Ok(count);
                        }
                    }
                    line_has_content = false;
                } else if !byte.is_ascii_whitespace() {
                    line_has_content = true;
                }
                continue;
            }
            if !array_started {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                if byte != b'[' {
                    return Ok(0);
                }
                array_started = true;
                depth = 1;
                continue;
            }
            if in_string {
                if escaped {
                    escaped = false;
                } else if byte == b'\\' {
                    escaped = true;
                } else if byte == b'"' {
                    in_string = false;
                }
                continue;
            }
            match byte {
                b'"' => {
                    element_started |= depth == 1;
                    in_string = true;
                }
                b'{' | b'[' => {
                    element_started |= depth == 1;
                    depth += 1;
                }
                b'}' if depth > 1 => depth -= 1,
                b']' if depth == 1 => {
                    if element_started {
                        count += 1;
                    }
                    return Ok(count);
                }
                b']' if depth > 1 => depth -= 1,
                b',' if depth == 1 => {
                    if element_started {
                        count += 1;
                        if count >= limit {
                            return Ok(count);
                        }
                    }
                    element_started = false;
                }
                byte if depth == 1 && !byte.is_ascii_whitespace() => element_started = true,
                _ => {}
            }
        }
    }
    if lines && line_has_content {
        count += 1;
    }
    Ok(count)
}

fn validate_observations(
    mut observations: Vec<StructureObservation>,
) -> ArpeggiaResult<Vec<StructureObservation>> {
    if observations.is_empty() {
        return Err(ArpeggiaError::InvalidArgument(
            "input contains no protein structures".into(),
        ));
    }
    observations.sort_unstable_by(|left, right| left.id.cmp(&right.id));
    let mut ids = BTreeSet::new();
    let mut paths = BTreeSet::new();
    for observation in &observations {
        if observation.id.is_empty() {
            return Err(ArpeggiaError::InvalidArgument(
                "structure IDs cannot be empty".into(),
            ));
        }
        if !ids.insert(&observation.id) {
            return Err(ArpeggiaError::InvalidArgument(format!(
                "duplicate structure ID: {}",
                observation.id
            )));
        }
        if !paths.insert(&observation.path) {
            return Err(ArpeggiaError::InvalidArgument(format!(
                "duplicate structure path: {}",
                observation.path.display()
            )));
        }
        if !observation.path.is_file() {
            return Err(ArpeggiaError::InvalidArgument(format!(
                "structure path is not a regular file: {}",
                observation.path.display()
            )));
        }
        if !is_structure_path(&observation.path) {
            return Err(ArpeggiaError::InvalidArgument(format!(
                "unsupported structure format: {}",
                observation.path.display()
            )));
        }
    }
    Ok(observations)
}

fn is_structure_path(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| {
            ["pdb", "cif", "mmcif"]
                .iter()
                .any(|expected| extension.eq_ignore_ascii_case(expected))
        })
}

fn checked_pair_count(n: usize) -> ArpeggiaResult<usize> {
    n.checked_mul(n.saturating_sub(1))
        .map(|value| value / 2)
        .ok_or_else(|| ArpeggiaError::InvalidArgument("pair count overflows usize".into()))
}

pub(crate) fn check_packed_matrix_memory(
    pair_count: usize,
    bypass_mem_check: bool,
) -> ArpeggiaResult<Vec<AnalysisWarning>> {
    let matrix_bytes = pair_count.checked_mul(size_of::<f64>()).ok_or_else(|| {
        ArpeggiaError::InvalidArgument("pairwise matrix size overflows usize".into())
    })? as u64;
    check_memory(matrix_bytes, bypass_mem_check, "packed RMSD matrix")
}

fn pair_at_index(index: usize, n: usize) -> (usize, usize) {
    let mut low = 1;
    let mut high = n;
    while low + 1 < high {
        let middle = (low + high) / 2;
        if middle * (middle - 1) / 2 <= index {
            low = middle;
        } else {
            high = middle;
        }
    }
    (low, index - low * (low - 1) / 2)
}

fn effective_thread_count(requested: usize) -> usize {
    let available = std::thread::available_parallelism().map_or(1, usize::from);
    if requested == 0 {
        available
    } else {
        requested.min(available)
    }
}

fn effective_available_memory() -> Option<u64> {
    let mut system = sysinfo::System::new();
    system.refresh_memory_specifics(sysinfo::MemoryRefreshKind::nothing().with_ram());
    let host = system.available_memory();
    if host == 0 {
        return None;
    }
    Some(
        system
            .cgroup_limits()
            .map(|limits| limits.free_memory)
            .filter(|available| *available > 0)
            .map_or(host, |available| host.min(available)),
    )
}

fn check_memory(
    estimated_bytes: u64,
    bypass: bool,
    label: &str,
) -> ArpeggiaResult<Vec<AnalysisWarning>> {
    if bypass {
        Ok(Vec::new())
    } else {
        memory_warnings_or_error(estimated_bytes, effective_available_memory(), label)
    }
}

fn memory_warnings_or_error(
    estimated_bytes: u64,
    available_bytes: Option<u64>,
    label: &str,
) -> ArpeggiaResult<Vec<AnalysisWarning>> {
    if let Some(available) = available_bytes {
        let ceiling = available.saturating_mul(4) / 5;
        if estimated_bytes > ceiling {
            return Err(ArpeggiaError::Calculation(format!(
                "estimated {label} memory {} exceeds the 80% available-RAM ceiling {}; use bypass_mem_check to proceed",
                format_bytes(estimated_bytes),
                format_bytes(ceiling)
            )));
        }
        Ok(Vec::new())
    } else if estimated_bytes > FALLBACK_WARNING_BYTES {
        Ok(vec![AnalysisWarning::new(
            WarningCode::MemoryEstimate,
            format!(
                "cannot query available RAM; estimated {label} memory is {} and excludes transient and clustering allocations",
                format_bytes(estimated_bytes)
            ),
        )])
    } else {
        Ok(Vec::new())
    }
}

fn format_bytes(bytes: u64) -> String {
    format!("{:.2} GiB", bytes as f64 / 1024_f64.powi(3))
}

fn missing_table_column(table: &str, column: &str) -> ArpeggiaError {
    ArpeggiaError::InvalidArgument(format!("{table} table requires column {column}"))
}

fn invalid_table_type(table: &str, column: &str, expected: &str) -> ArpeggiaError {
    ArpeggiaError::InvalidArgument(format!("{table} column {column} must have type {expected}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Copy)]
    enum TestTableFormat {
        Csv,
        Parquet,
        Json,
        NDJson,
    }

    impl TestTableFormat {
        fn extension(self) -> &'static str {
            match self {
                Self::Csv => "csv",
                Self::Parquet => "parquet",
                Self::Json => "json",
                Self::NDJson => "ndjson",
            }
        }
    }

    fn write_test_table(dataframe: &mut DataFrame, base: &Path, format: TestTableFormat) {
        let mut file = File::create(base.with_extension(format.extension())).unwrap();
        match format {
            TestTableFormat::Csv => CsvWriter::new(&mut file).finish(dataframe).unwrap(),
            TestTableFormat::Parquet => {
                ParquetWriter::new(&mut file).finish(dataframe).unwrap();
            }
            TestTableFormat::Json => JsonWriter::new(&mut file)
                .with_json_format(JsonFormat::Json)
                .finish(dataframe)
                .unwrap(),
            TestTableFormat::NDJson => JsonWriter::new(&mut file)
                .with_json_format(JsonFormat::JsonLines)
                .finish(dataframe)
                .unwrap(),
        };
    }

    fn one_dimensional_matrix(positions: &[f64]) -> PairwiseRmsdMatrix {
        let mut data = Vec::new();
        for row in 1..positions.len() {
            for column in 0..row {
                data.push((positions[row] - positions[column]).abs());
            }
        }
        PairwiseRmsdMatrix {
            ids: (b'a'..)
                .take(positions.len())
                .map(char::from)
                .map(String::from)
                .collect(),
            data,
        }
    }

    #[test]
    fn pair_index_matches_lower_triangle_layout() {
        let expected = [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)];
        for (index, pair) in expected.into_iter().enumerate() {
            assert_eq!(pair_at_index(index, 4), pair);
        }
    }

    #[test]
    fn pairwise_dataframe_is_long_and_roundtrips() {
        let matrix = PairwiseRmsdMatrix {
            ids: vec!["a".into(), "b".into(), "c".into()],
            data: vec![1.0, 2.0, 3.0],
        };
        let dataframe = matrix.to_dataframe().unwrap();
        assert_eq!(dataframe.get_column_names(), ["id_1", "id_2", "rmsd"]);
        let roundtrip = PairwiseRmsdMatrix::from_dataframe(&dataframe, None, 3).unwrap();
        assert_eq!(roundtrip.ids, matrix.ids);
        assert_eq!(roundtrip.data, matrix.data);
    }

    #[test]
    fn pairwise_dataframe_normalizes_row_and_id_order() {
        let dataframe = df!(
            "id_1" => ["c", "a", "c"],
            "id_2" => ["a", "b", "b"],
            "rmsd" => [2.0, 1.0, 3.0]
        )
        .unwrap();
        let matrix = PairwiseRmsdMatrix::from_dataframe(&dataframe, None, 3).unwrap();
        assert_eq!(matrix.ids, ["a", "b", "c"]);
        assert_eq!(matrix.data, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn pairwise_dataframe_accepts_fragmented_columns() {
        let mut dataframe = df!("id_1" => ["a"], "id_2" => ["b"], "rmsd" => [1.0]).unwrap();
        dataframe
            .vstack_mut(&df!("id_1" => ["a"], "id_2" => ["c"], "rmsd" => [2.0]).unwrap())
            .unwrap();
        dataframe
            .vstack_mut(&df!("id_1" => ["b"], "id_2" => ["c"], "rmsd" => [3.0]).unwrap())
            .unwrap();
        assert!(dataframe.column("id_1").unwrap().n_chunks() > 1);
        let matrix = PairwiseRmsdMatrix::from_dataframe(&dataframe, None, 3).unwrap();
        assert_eq!(matrix.data, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn pairwise_dataframe_rejects_duplicate_pairs() {
        let dataframe = df!(
            "id_1" => ["a", "b", "a"],
            "id_2" => ["b", "a", "c"],
            "rmsd" => [1.0, 1.0, 2.0]
        )
        .unwrap();
        assert!(PairwiseRmsdMatrix::from_dataframe(&dataframe, None, 2).is_err());
    }

    #[test]
    fn pairwise_dataframe_requires_numeric_rmsd() {
        let dataframe = df!(
            "id_1" => ["a", "a", "b"],
            "id_2" => ["b", "c", "c"],
            "rmsd" => [true, false, true]
        )
        .unwrap();
        assert!(PairwiseRmsdMatrix::from_dataframe(&dataframe, None, 3).is_err());
    }

    #[test]
    fn pairwise_dataframe_rejects_invalid_distances() {
        for invalid in [None, Some(-1.0), Some(f64::NAN), Some(f64::INFINITY)] {
            let dataframe = df!(
                "id_1" => ["a", "a", "b"],
                "id_2" => ["b", "c", "c"],
                "rmsd" => [invalid, Some(2.0), Some(3.0)]
            )
            .unwrap();
            assert!(PairwiseRmsdMatrix::from_dataframe(&dataframe, None, 3).is_err());
        }
    }

    #[test]
    fn pairwise_dataframe_schema_errors_are_invalid_arguments() {
        let missing = df!("wrong" => [1_u32]).unwrap();
        assert!(matches!(
            PairwiseRmsdMatrix::from_dataframe(&missing, None, 3),
            Err(ArpeggiaError::InvalidArgument(_))
        ));
        let wrong_type = df!(
            "id_1" => [1_u32, 1, 2],
            "id_2" => [2_u32, 3, 3],
            "rmsd" => [1.0, 2.0, 3.0]
        )
        .unwrap();
        assert!(matches!(
            PairwiseRmsdMatrix::from_dataframe(&wrong_type, None, 3),
            Err(ArpeggiaError::InvalidArgument(_))
        ));
    }

    #[test]
    fn pairwise_cache_rejects_impossible_height_before_scanning_ids() {
        let dataframe = df!(
            "id_1" => ["bad"],
            "id_2" => ["bad"],
            "rmsd" => [0.0]
        )
        .unwrap();
        let expected = vec!["a".into(), "b".into(), "c".into()];
        assert!(matches!(
            PairwiseRmsdMatrix::from_dataframe(&dataframe, Some(&expected), 3),
            Err(ArpeggiaError::InvalidArgument(message)) if message.contains("3 expected IDs")
        ));
    }

    #[test]
    fn cached_table_readers_preflight_wrong_height() {
        let directory =
            std::env::temp_dir().join(format!("arpeggia-cache-height-{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        let table = df!(
            "id_1" => ["a", "a", "b", "a,quoted"],
            "id_2" => ["b", "c", "c", "c"],
            "rmsd" => [1.0, 2.0, 3.0, 4.0]
        )
        .unwrap();
        let expected = vec!["a".into(), "b".into(), "c".into()];
        for format in [
            TestTableFormat::Csv,
            TestTableFormat::Parquet,
            TestTableFormat::Json,
            TestTableFormat::NDJson,
        ] {
            let base = directory.join(format!("cache-{}", format.extension()));
            write_test_table(&mut table.clone(), &base, format);
            assert!(
                read_pairwise_matrix(&base.with_extension(format.extension()), &expected, false,)
                    .is_err()
            );
        }
    }

    #[test]
    fn ndjson_cache_preflight_returns_an_immutable_snapshot() {
        let path = std::env::temp_dir().join(format!(
            "arpeggia-ndjson-snapshot-{}.ndjson",
            std::process::id()
        ));
        std::fs::write(&path, "{\"id\":\"a\"}\n{\"id\":\"b\"}\n").unwrap();
        let mut file = File::open(&path).unwrap();
        let bytes = preflight_json_rows(&mut file, true, Some(2))
            .unwrap()
            .unwrap();
        std::fs::write(&path, "{\"id\":\"replacement\"}\n").unwrap();
        let dataframe = JsonReader::new(Cursor::new(bytes))
            .with_json_format(JsonFormat::JsonLines)
            .finish()
            .unwrap();
        assert_eq!(dataframe.height(), 2);
    }

    #[test]
    fn memory_guard_has_a_bypass_and_fallback_warning() {
        assert!(memory_warnings_or_error(90, Some(100), "test").is_err());
        assert!(check_memory(u64::MAX, true, "test").is_ok());
        assert_eq!(
            memory_warnings_or_error(FALLBACK_WARNING_BYTES + 1, None, "test").unwrap()[0].code,
            WarningCode::MemoryEstimate
        );
    }

    #[test]
    fn directory_input_is_sorted_and_case_insensitive() {
        let directory = std::env::temp_dir().join(format!(
            "arpeggia-ensemble-directory-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&directory).unwrap();
        std::fs::write(directory.join("b.PDB"), "END\n").unwrap();
        std::fs::write(directory.join("A.cIf"), "data_A\n#\n").unwrap();
        std::fs::write(directory.join("ignored.txt"), "ignored\n").unwrap();
        let observations = read_structure_observations(&directory, "id", "path").unwrap();
        assert_eq!(
            observations
                .iter()
                .map(|observation| observation.id.as_str())
                .collect::<Vec<_>>(),
            ["A", "b"]
        );
    }

    #[test]
    fn manifests_support_all_table_formats_and_relative_paths() {
        let directory = std::env::temp_dir().join(format!(
            "arpeggia-ensemble-manifests-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&directory).unwrap();
        std::fs::write(directory.join("a.pdb"), "END\n").unwrap();
        std::fs::write(directory.join("b.cif"), "data_b\n#\n").unwrap();
        let manifest = df!("name" => ["b", "a"], "file" => ["b.cif", "a.pdb"]).unwrap();
        for format in [
            TestTableFormat::Csv,
            TestTableFormat::Parquet,
            TestTableFormat::Json,
            TestTableFormat::NDJson,
        ] {
            let base = directory.join(format!("manifest-{}", format.extension()));
            write_test_table(&mut manifest.clone(), &base, format);
            let path = base.with_extension(format.extension());
            let observations = read_structure_observations(&path, "name", "file").unwrap();
            assert_eq!(
                observations
                    .iter()
                    .map(|observation| observation.id.as_str())
                    .collect::<Vec<_>>(),
                ["a", "b"]
            );
        }
    }

    #[test]
    fn manifest_structure_paths_must_be_regular_files() {
        let directory = std::env::temp_dir().join(format!(
            "arpeggia-non-regular-manifest-{}",
            std::process::id()
        ));
        let structure_directory = directory.join("not-a-file.pdb");
        std::fs::create_dir_all(&structure_directory).unwrap();
        let mut manifest = df!("id" => ["a"], "path" => ["not-a-file.pdb"]).unwrap();
        let manifest_path = directory.join("manifest");
        write_test_table(&mut manifest, &manifest_path, TestTableFormat::Csv);
        assert!(matches!(
            read_structure_observations(
                &manifest_path.with_extension("csv"),
                "id",
                "path"
            ),
            Err(ArpeggiaError::InvalidArgument(message))
                if message.contains("not a regular file")
        ));
    }

    #[test]
    fn pairwise_rmsd_runs_on_repeated_structure_files() {
        let source = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
        let directory = std::env::temp_dir().join(format!(
            "arpeggia-pairwise-directory-{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&directory).unwrap();
        for id in ["a", "b", "c"] {
            std::fs::copy(&source, directory.join(format!("{id}.pdb"))).unwrap();
        }
        let options = PairwiseRmsdOptions {
            residues: "A:1-20".into(),
            num_threads: 2,
            ..PairwiseRmsdOptions::default()
        };
        let observations = read_structure_observations(&directory, "id", "path").unwrap();
        let selector = ResidueSelector::parse(&options.residues).unwrap();
        let reference = load_observation(&observations[0], &options, &selector, None).unwrap();
        let other = load_observation(
            &observations[1],
            &options,
            &selector,
            reference.keys.as_deref(),
        )
        .unwrap();
        assert!(other.keys.is_none());
        let analysis = get_pairwise_rmsd(&directory, "id", "path", &options).unwrap();
        assert_eq!(analysis.value.height(), 3);
        assert!(
            analysis
                .value
                .column("rmsd")
                .unwrap()
                .f64()
                .unwrap()
                .into_no_null_iter()
                .all(|value| value < 1e-12)
        );
        let matrix = PairwiseRmsdMatrix::from_dataframe(&analysis.value, None, 3).unwrap();
        let clusters = cluster_pairwise_rmsd(
            &matrix,
            &ClusterOptions {
                max_clusters: Some(2),
                ..ClusterOptions::default()
            },
        )
        .unwrap()
        .value;
        assert!(
            clusters
                .column("cluster_id")
                .unwrap()
                .u32()
                .unwrap()
                .into_no_null_iter()
                .all(|cluster| cluster == 0)
        );
    }

    #[test]
    fn pairwise_matrix_requires_two_observations() {
        assert!(matches!(
            get_pairwise_rmsd_matrix(&[], &PairwiseRmsdOptions::default()),
            Err(ArpeggiaError::InvalidArgument(_))
        ));
    }

    #[test]
    fn pairwise_rmsd_is_identical_across_thread_counts() {
        let source = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
        let directory =
            std::env::temp_dir().join(format!("arpeggia-pairwise-threads-{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        for id in ["a", "b", "c", "d"] {
            std::fs::copy(&source, directory.join(format!("{id}.pdb"))).unwrap();
        }
        let mut observations = read_structure_observations(&directory, "id", "path").unwrap();
        observations.reverse();
        let options = PairwiseRmsdOptions {
            residues: "A:1-20".into(),
            num_threads: 1,
            ..PairwiseRmsdOptions::default()
        };
        let serial = get_pairwise_rmsd_matrix(&observations, &options)
            .unwrap()
            .value;
        let parallel = get_pairwise_rmsd_matrix(
            &observations,
            &PairwiseRmsdOptions {
                num_threads: 4,
                ..options.clone()
            },
        )
        .unwrap()
        .value;
        assert_eq!(serial.ids, parallel.ids);
        assert_eq!(serial.data, parallel.data);
        assert_eq!(serial.ids, ["a", "b", "c", "d"]);

        observations.push(observations[0].clone());
        assert!(get_pairwise_rmsd_matrix(&observations, &options).is_err());
    }

    #[test]
    fn requested_threads_are_capped_to_available_processors() {
        let available = std::thread::available_parallelism().map_or(1, usize::from);
        assert_eq!(effective_thread_count(0), available);
        assert_eq!(effective_thread_count(usize::MAX), available);
    }

    #[test]
    fn fixed_k_medoids_separates_two_groups() {
        let matrix = PairwiseRmsdMatrix {
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
        let matrix = PairwiseRmsdMatrix {
            ids: ["a", "b", "c", "d", "e", "f"]
                .into_iter()
                .map(str::to_string)
                .collect(),
            data: vec![
                0.1, 0.2, 0.1, 5.0, 5.1, 5.0, 5.1, 5.0, 5.1, 0.1, 5.0, 5.1, 5.0, 0.2, 0.1,
            ],
        };
        let scaled = PairwiseRmsdMatrix {
            ids: matrix.ids.clone(),
            data: matrix.data.iter().map(|value| value * 1e307).collect(),
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
        let matrix = PairwiseRmsdMatrix {
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
        let matrix = PairwiseRmsdMatrix {
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
        let matrix = PairwiseRmsdMatrix {
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
        let matrix = PairwiseRmsdMatrix {
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
        let matrix = PairwiseRmsdMatrix {
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
        let matrix = PairwiseRmsdMatrix {
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
    fn pairwise_matrix_get_checks_indices() {
        let matrix = PairwiseRmsdMatrix {
            ids: vec!["a".into(), "b".into(), "c".into()],
            data: vec![1.0, 2.0, 3.0],
        };
        assert_eq!(matrix.get(0, 2), Some(2.0));
        assert_eq!(matrix.get(3, 3), None);
        assert_eq!(matrix.get(0, 3), None);
    }

    #[test]
    fn fixed_k_medoids_is_deterministic_for_ties() {
        let matrix = PairwiseRmsdMatrix {
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
        let matrix = PairwiseRmsdMatrix {
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
        let matrix = PairwiseRmsdMatrix {
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
        let matrix = PairwiseRmsdMatrix {
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
        let matrix = PairwiseRmsdMatrix {
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
