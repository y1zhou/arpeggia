//! Ensemble structure input, pairwise RMSD, and clustering.

use crate::rmsd::{
    AtomIdentity, ResidueSelector, kabsch_centered_rmsd, prepare_centered_coordinates,
    select_coordinates, validate_correspondence,
};
use crate::{
    Analysis, AnalysisWarning, ArpeggiaError, ArpeggiaResult, AtomSubset, WarningCode, load_model,
};
use clap::ValueEnum;
use kmedoids::ArrayAdapter;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::path::{Path, PathBuf};

const FALLBACK_WARNING_BYTES: u64 = 8 * 1024 * 1024 * 1024;

/// One named structure participating in an ensemble calculation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StructureObservation {
    /// Case-sensitive identifier used in result tables.
    pub id: String,
    /// Canonical path to the structure file.
    pub path: PathBuf,
}

/// Options controlling pairwise structural superposition.
#[derive(Clone, Debug)]
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

impl Default for PairwiseRmsdOptions {
    fn default() -> Self {
        Self {
            model_num: 0,
            residues: String::new(),
            atoms: AtomSubset::Ca,
            num_threads: 0,
            bypass_mem_check: false,
        }
    }
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
    if observations.len() < 2 {
        return Err(ArpeggiaError::InvalidArgument(
            "pairwise RMSD requires at least two structures".into(),
        ));
    }
    get_pairwise_rmsd_matrix(&observations, options).and_then(|analysis| {
        analysis
            .value
            .to_dataframe()
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
    pub fn get(&self, left: usize, right: usize) -> f64 {
        if left == right {
            0.0
        } else {
            let (row, column) = if left > right {
                (left, right)
            } else {
                (right, left)
            };
            self.data[row * (row - 1) / 2 + column]
        }
    }

    /// Convert the packed matrix to one row per unordered pair.
    pub fn to_dataframe(&self) -> ArpeggiaResult<DataFrame> {
        let mut id_1 = Vec::with_capacity(self.data.len());
        let mut id_2 = Vec::with_capacity(self.data.len());
        for row in 1..self.ids.len() {
            for column in 0..row {
                id_1.push(self.ids[column].as_str());
                id_2.push(self.ids[row].as_str());
            }
        }
        df!("id_1" => id_1, "id_2" => id_2, "rmsd" => self.data.clone()).map_err(polars_error)
    }

    /// Validate a complete long pairwise table and pack it for clustering.
    pub fn from_dataframe(
        dataframe: &DataFrame,
        expected_ids: Option<&[String]>,
        minimum_ids: usize,
    ) -> ArpeggiaResult<Self> {
        let left = dataframe
            .column("id_1")
            .map_err(polars_error)?
            .str()
            .map_err(polars_error)?;
        let right = dataframe
            .column("id_2")
            .map_err(polars_error)?
            .str()
            .map_err(polars_error)?;
        let rmsd = dataframe
            .column("rmsd")
            .map_err(polars_error)?
            .cast(&DataType::Float64)
            .map_err(polars_error)?;
        let rmsd = rmsd.f64().map_err(polars_error)?;

        let mut ids = BTreeSet::new();
        for row in 0..dataframe.height() {
            let left = left.get(row).ok_or_else(|| {
                ArpeggiaError::InvalidArgument("pairwise id_1 contains null values".into())
            })?;
            let right = right.get(row).ok_or_else(|| {
                ArpeggiaError::InvalidArgument("pairwise id_2 contains null values".into())
            })?;
            if left.is_empty() || right.is_empty() {
                return Err(ArpeggiaError::InvalidArgument(
                    "pairwise IDs cannot be empty".into(),
                ));
            }
            ids.insert(left.to_string());
            ids.insert(right.to_string());
        }
        let ids = ids.into_iter().collect::<Vec<_>>();
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
            .collect::<BTreeMap<_, _>>();
        let pair_count = checked_pair_count(ids.len())?;
        if dataframe.height() != pair_count {
            return Err(ArpeggiaError::InvalidArgument(format!(
                "pairwise table has {} rows; {} IDs require exactly {pair_count}",
                dataframe.height(),
                ids.len()
            )));
        }
        let mut data = vec![f64::NAN; pair_count];
        for input_row in 0..dataframe.height() {
            let left_id = left.get(input_row).expect("null values were rejected");
            let right_id = right.get(input_row).expect("null values were rejected");
            let value = rmsd.get(input_row).ok_or_else(|| {
                ArpeggiaError::InvalidArgument("pairwise rmsd contains null values".into())
            })?;
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
            let (row, column) = if left_index > right_index {
                (left_index, right_index)
            } else {
                (right_index, left_index)
            };
            let output_index = row * (row - 1) / 2 + column;
            if !data[output_index].is_nan() {
                return Err(ArpeggiaError::InvalidArgument(format!(
                    "pairwise table contains duplicate pair {left_id}/{right_id}"
                )));
            }
            data[output_index] = value;
        }
        if data.iter().any(|value| value.is_nan()) {
            return Err(ArpeggiaError::InvalidArgument(
                "pairwise table does not contain every unordered ID pair".into(),
            ));
        }
        Ok(Self { ids, data })
    }
}

impl ArrayAdapter<f64> for PairwiseRmsdMatrix {
    fn len(&self) -> usize {
        self.ids.len()
    }

    fn is_square(&self) -> bool {
        checked_pair_count(self.ids.len()).is_ok_and(|count| count == self.data.len())
    }

    fn get(&self, left: usize, right: usize) -> f64 {
        PairwiseRmsdMatrix::get(self, left, right)
    }
}

/// Read and validate a cached pairwise RMSD table for the expected IDs.
pub fn read_pairwise_matrix(
    path: &Path,
    expected_ids: &[String],
) -> ArpeggiaResult<PairwiseRmsdMatrix> {
    let dataframe = read_dataframe(path)?;
    PairwiseRmsdMatrix::from_dataframe(&dataframe, Some(expected_ids), 2)
}

/// Calculate every unordered RMSD pair into a packed matrix.
pub fn get_pairwise_rmsd_matrix(
    observations: &[StructureObservation],
    options: &PairwiseRmsdOptions,
) -> ArpeggiaResult<Analysis<PairwiseRmsdMatrix>> {
    let pair_count = checked_pair_count(observations.len())?;
    let matrix_bytes = pair_count.checked_mul(size_of::<f64>()).ok_or_else(|| {
        ArpeggiaError::InvalidArgument("pairwise matrix size overflows usize".into())
    })? as u64;
    let mut warnings = check_memory(matrix_bytes, options.bypass_mem_check, "packed RMSD matrix")?;

    let selector = ResidueSelector::parse(&options.residues)?;
    let first_loaded = load_observation(&observations[0], options, &selector, None)?;
    let reference_keys = first_loaded.keys;
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
    let parse_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(parse_threads)
        .build()
        .map_err(|error| ArpeggiaError::Calculation(error.to_string()))?;
    let remaining = parse_pool.install(|| {
        observations[1..]
            .par_iter()
            .map(|observation| {
                load_observation(observation, options, &selector, Some(&reference_keys))
            })
            .collect::<Vec<_>>()
    });

    let mut coordinates = Vec::with_capacity(observations.len());
    coordinates.push(first_loaded.coordinates);
    for result in remaining {
        let prepared = result?;
        warnings.extend(prepared.warnings);
        coordinates.push(prepared.coordinates);
    }

    let pair_threads = effective_thread_count(options.num_threads);
    let chunk_size = pair_count
        .div_ceil(pair_threads.saturating_mul(8).max(1))
        .max(1);
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
                    *value = kabsch_centered_rmsd(&coordinates[column], &coordinates[row])?;
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
    if n < 3 {
        return Err(ArpeggiaError::InvalidArgument(
            "structure clustering requires at least three structures".into(),
        ));
    }
    if n > u32::MAX as usize {
        return Err(ArpeggiaError::InvalidArgument(
            "structure count exceeds UInt32 cluster output capacity".into(),
        ));
    }
    if options.max_iterations == 0 {
        return Err(ArpeggiaError::InvalidArgument(
            "max_iterations must be positive".into(),
        ));
    }

    let mut warnings = Vec::new();
    let requested_clusters = if let Some(k) = options.num_clusters {
        if options.max_clusters.is_some() {
            warnings.push(AnalysisWarning::new(
                WarningCode::ArgumentIgnored,
                "num_clusters takes priority; max_clusters was ignored",
            ));
        }
        ClusterCount::Fixed(k)
    } else if let Some(k) = options.max_clusters {
        ClusterCount::Automatic(k)
    } else {
        return Err(ArpeggiaError::InvalidArgument(
            "one of num_clusters or max_clusters is required".into(),
        ));
    };

    let (medoids, assignments) = match requested_clusters {
        ClusterCount::Fixed(k) => {
            validate_cluster_count(k, n, "num_clusters")?;
            run_fasterpam(matrix, k, options.max_iterations)?
        }
        ClusterCount::Automatic(max_k) => {
            if max_k < 2 || max_k > n {
                return Err(ArpeggiaError::InvalidArgument(format!(
                    "max_clusters must be between 2 and {n}"
                )));
            }
            if matrix.data.iter().all(|distance| *distance == 0.0) {
                (vec![0], vec![0; n])
            } else {
                run_dynmsc(matrix, max_k, options.max_iterations)?
            }
        }
    };
    let value = cluster_dataframe(matrix, &medoids, &assignments)?;
    Ok(Analysis::new(value, warnings))
}

enum ClusterCount {
    Fixed(usize),
    Automatic(usize),
}

fn validate_cluster_count(k: usize, n: usize, name: &str) -> ArpeggiaResult<()> {
    if k == 0 || k > n {
        Err(ArpeggiaError::InvalidArgument(format!(
            "{name} must be between 1 and {n}"
        )))
    } else {
        Ok(())
    }
}

fn run_fasterpam(
    matrix: &PairwiseRmsdMatrix,
    k: usize,
    max_iterations: usize,
) -> ArpeggiaResult<(Vec<usize>, Vec<usize>)> {
    let (_, _, mut medoids) = kmedoids::pam_build::<_, f64, f64>(matrix, k);
    let (_, assignments, iterations, _) =
        kmedoids::fasterpam::<_, f64, f64>(matrix, &mut medoids, max_iterations);
    if iterations >= max_iterations {
        return Err(ArpeggiaError::Calculation(format!(
            "k-medoids reached the {max_iterations}-iteration limit"
        )));
    }
    Ok((medoids, assignments))
}

fn run_dynmsc(
    matrix: &PairwiseRmsdMatrix,
    max_k: usize,
    max_iterations: usize,
) -> ArpeggiaResult<(Vec<usize>, Vec<usize>)> {
    let (_, _, initial_medoids) = kmedoids::pam_build::<_, f64, f64>(matrix, max_k);
    let (_, assignments, iterations, _, medoids, _) =
        kmedoids::dynmsc::<_, f64, f64>(matrix, &initial_medoids, 2, max_iterations);
    if iterations >= max_iterations {
        return Err(ArpeggiaError::Calculation(format!(
            "automatic k-medoids reached the {max_iterations}-iteration limit"
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
        rmsd_to_medoid.push(matrix.get(index, medoid));
    }
    df!(
        "id" => matrix.ids.iter().map(String::as_str).collect::<Vec<_>>(),
        "cluster_id" => output_cluster_ids,
        "medoid_id" => medoid_ids,
        "rmsd_to_medoid" => rmsd_to_medoid,
    )
    .map_err(polars_error)
}

struct PreparedObservation {
    keys: Vec<AtomIdentity>,
    coordinates: Vec<[f64; 3]>,
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
    if let Some(reference_keys) = reference_keys {
        validate_correspondence(reference_keys, &selected.keys)?;
    }
    let coordinates = prepare_centered_coordinates(selected.coordinates)?;
    let mut warnings = loaded.warnings;
    warnings.extend(selected.warnings);
    Ok(PreparedObservation {
        keys: selected.keys,
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
    let dataframe = read_dataframe(manifest)?;
    let ids = dataframe
        .column(id_column)
        .map_err(polars_error)?
        .str()
        .map_err(polars_error)?;
    let paths = dataframe
        .column(path_column)
        .map_err(polars_error)?
        .str()
        .map_err(polars_error)?;
    let base = manifest.parent().unwrap_or_else(|| Path::new("."));
    (0..dataframe.height())
        .map(|row| {
            let id = ids.get(row).ok_or_else(|| {
                ArpeggiaError::InvalidArgument(format!(
                    "manifest column {id_column} contains null at row {row}"
                ))
            })?;
            let path = paths.get(row).ok_or_else(|| {
                ArpeggiaError::InvalidArgument(format!(
                    "manifest column {path_column} contains null at row {row}"
                ))
            })?;
            if id.is_empty() || path.is_empty() {
                return Err(ArpeggiaError::InvalidArgument(format!(
                    "manifest row {row} contains an empty ID or path"
                )));
            }
            let path = Path::new(path);
            let path = if path.is_absolute() {
                path.to_path_buf()
            } else {
                base.join(path)
            };
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

pub(crate) fn read_dataframe(path: &Path) -> ArpeggiaResult<DataFrame> {
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
    let file = File::open(path)?;
    match extension.as_str() {
        "csv" => CsvReadOptions::default()
            .into_reader_with_file_handle(file)
            .finish()
            .map_err(polars_error),
        "parquet" => ParquetReader::new(file).finish().map_err(polars_error),
        "json" => JsonReader::new(file)
            .with_json_format(JsonFormat::Json)
            .finish()
            .map_err(polars_error),
        "ndjson" => JsonReader::new(file)
            .with_json_format(JsonFormat::JsonLines)
            .finish()
            .map_err(polars_error),
        _ => Err(ArpeggiaError::InvalidArgument(format!(
            "unsupported table format .{extension}; expected csv, parquet, json, or ndjson"
        ))),
    }
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
            matches!(
                extension.to_ascii_lowercase().as_str(),
                "pdb" | "cif" | "mmcif"
            )
        })
}

fn checked_pair_count(n: usize) -> ArpeggiaResult<usize> {
    n.checked_mul(n.saturating_sub(1))
        .map(|value| value / 2)
        .ok_or_else(|| ArpeggiaError::InvalidArgument("pair count overflows usize".into()))
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
    if requested == 0 {
        std::thread::available_parallelism().map_or(1, usize::from)
    } else {
        requested
    }
    .max(1)
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
        memory_warnings_or_error(estimated_bytes, effective_available_memory(), false, label)
    }
}

fn memory_warnings_or_error(
    estimated_bytes: u64,
    available_bytes: Option<u64>,
    bypass: bool,
    label: &str,
) -> ArpeggiaResult<Vec<AnalysisWarning>> {
    if bypass {
        return Ok(Vec::new());
    }
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

fn polars_error(error: PolarsError) -> ArpeggiaError {
    ArpeggiaError::Io(std::io::Error::other(error.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn memory_guard_has_a_bypass_and_fallback_warning() {
        assert!(memory_warnings_or_error(90, Some(100), false, "test").is_err());
        assert!(memory_warnings_or_error(90, Some(100), true, "test").is_ok());
        assert_eq!(
            memory_warnings_or_error(FALLBACK_WARNING_BYTES + 1, None, false, "test").unwrap()[0]
                .code,
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
        let analysis = get_pairwise_rmsd(
            &directory,
            "id",
            "path",
            &PairwiseRmsdOptions {
                residues: "A:1-20".into(),
                num_threads: 2,
                ..PairwiseRmsdOptions::default()
            },
        )
        .unwrap();
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
    fn automatic_k_medoids_short_circuits_identical_structures() {
        let matrix = PairwiseRmsdMatrix {
            ids: vec!["a".into(), "b".into(), "c".into()],
            data: vec![0.0; 3],
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
}
