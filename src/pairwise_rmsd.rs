//! Ensemble structure input and pairwise RMSD calculation.

use crate::rmsd::{
    AtomIdentity, PreparedCoordinates, ResidueSelector, kabsch_prepared_rmsd, prepare_coordinates,
    select_coordinates, validate_correspondence,
};
use crate::utils::polars_calculation_error;
use crate::{
    Analysis, AnalysisWarning, ArpeggiaError, ArpeggiaResult, AtomSubset, WarningCode, load_model,
};
use kmedoids::ArrayAdapter;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};

const FALLBACK_WARNING_BYTES: u64 = 8 * 1024 * 1024 * 1024;

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

/// Symmetric pairwise RMSD matrix stored as a packed lower triangle.
#[derive(Clone, Debug)]
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

    #[cfg(test)]
    pub(crate) fn from_test_data(ids: Vec<String>, data: Vec<f64>) -> Self {
        Self { ids, data }
    }

    #[cfg(test)]
    pub(crate) fn packed_data(&self) -> &[f64] {
        &self.data
    }

    /// RMSD between two matrix indices.
    pub fn get(&self, left: usize, right: usize) -> Option<f64> {
        (left < self.len() && right < self.len()).then(|| self.distance(left, right))
    }

    pub(crate) fn distance(&self, left: usize, right: usize) -> f64 {
        if left == right {
            0.0
        } else {
            let (row, column) = (left.max(right), left.min(right));
            self.data[row * (row - 1) / 2 + column]
        }
    }

    pub(crate) fn distance_range(&self) -> (f64, Option<f64>) {
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
            .into_reader_with_file_handle(std::fs::File::open(path)?)
            .finish()
            .map_err(polars_input_error),
        "parquet" => {
            let mut reader = ParquetReader::new(std::fs::File::open(path)?);
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
        "ndjson" => {
            LazyJsonLineReader::new(PlRefPath::try_from_path(path).map_err(polars_input_error)?)
                .with_n_rows(expected_rows.map(|rows| rows.saturating_add(1)))
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
        _ => Err(ArpeggiaError::InvalidArgument(format!(
            "unsupported table format .{extension}; expected csv, parquet, or ndjson"
        ))),
    }
}

fn polars_input_error(error: PolarsError) -> ArpeggiaError {
    ArpeggiaError::InvalidArgument(error.to_string())
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
    use crate::clustering::{ClusterOptions, cluster_pairwise_rmsd};

    #[derive(Clone, Copy)]
    enum TestTableFormat {
        Csv,
        Parquet,
        NDJson,
    }

    impl TestTableFormat {
        fn extension(self) -> &'static str {
            match self {
                Self::Csv => "csv",
                Self::Parquet => "parquet",
                Self::NDJson => "ndjson",
            }
        }
    }

    fn write_test_table(dataframe: &mut DataFrame, base: &Path, format: TestTableFormat) {
        let mut file = std::fs::File::create(base.with_extension(format.extension())).unwrap();
        match format {
            TestTableFormat::Csv => CsvWriter::new(&mut file).finish(dataframe).unwrap(),
            TestTableFormat::Parquet => {
                ParquetWriter::new(&mut file).finish(dataframe).unwrap();
            }
            TestTableFormat::NDJson => JsonWriter::new(&mut file)
                .with_json_format(JsonFormat::JsonLines)
                .finish(dataframe)
                .unwrap(),
        };
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
    fn cached_table_readers_reject_wrong_height() {
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
    fn pairwise_matrix_get_checks_indices() {
        let matrix = PairwiseRmsdMatrix {
            ids: vec!["a".into(), "b".into(), "c".into()],
            data: vec![1.0, 2.0, 3.0],
        };
        assert_eq!(matrix.get(0, 2), Some(2.0));
        assert_eq!(matrix.get(3, 3), None);
        assert_eq!(matrix.get(0, 3), None);
    }
}
