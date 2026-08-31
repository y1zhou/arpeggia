use crate::{ArpeggiaError, ArpeggiaResult};
use polars::prelude::*;
use std::path::{Path, PathBuf};

pub(crate) fn string_values<'a>(values: impl Iterator<Item = &'a str>) -> StringChunked {
    StringChunked::from_iter_values(PlSmallStr::EMPTY, values)
}

pub(crate) fn polars_calculation_error(error: PolarsError) -> ArpeggiaError {
    ArpeggiaError::Calculation(error.to_string())
}

pub(crate) fn polars_input_error(error: PolarsError) -> ArpeggiaError {
    ArpeggiaError::InvalidArgument(error.to_string())
}

/// Execute a parallel operation with the configured thread limit.
/// Uses rayon's thread pool with the specified number of threads.
pub fn run_with_threads<F, T>(num_threads: isize, f: F) -> T
where
    F: FnOnce() -> T + Send,
    T: Send,
{
    let n_threads = num_threads.max(0) as usize;
    if n_threads == 0 || n_threads >= rayon::current_num_threads() {
        // Use global pool directly when auto or enough threads
        f()
    } else {
        // Create a scoped pool with limited threads
        match rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build()
        {
            Ok(pool) => pool.install(f),
            Err(_) => f(),
        }
    }
}

/// Prepare output directory for `DataFrame` output
///
/// # Arguments
///
/// * `output` - Output directory or file path
/// * `filename` - Base filename for the output file
/// * `output_format` - File format for the output file
///
/// # Returns
///
/// The canonicalized output file path with the correct file extension.
pub fn prepare_df_output_dir(
    output: &PathBuf,
    filename: &str,
    output_format: DataFrameFileType,
) -> ArpeggiaResult<PathBuf> {
    let mut components = Path::new(filename).components();
    let valid_filename = matches!(components.next(), Some(std::path::Component::Normal(_)))
        && components.next().is_none();
    if !valid_filename {
        return Err(ArpeggiaError::InvalidArgument(
            "output filename must be one normal path component".into(),
        ));
    }
    let output_path = std::path::absolute(output)?;
    std::fs::create_dir_all(&output_path)?;
    Ok(output_path
        .join(filename)
        .with_extension(output_format.to_string()))
}

/// Write a `DataFrame` to a CSV file
///
pub fn write_df_to_file(
    df: &mut DataFrame,
    file_path: &Path,
    file_type: DataFrameFileType,
) -> ArpeggiaResult<()> {
    let file_suffix = file_type.to_string();
    let mut file = std::fs::File::create(file_path.with_extension(file_suffix))?;
    write_df(df, &mut file, file_type)
}

/// Write a `DataFrame` only when the destination does not already exist.
pub fn write_df_to_new_file(
    df: &mut DataFrame,
    file_path: &Path,
    file_type: DataFrameFileType,
) -> ArpeggiaResult<()> {
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(file_path.with_extension(file_type.to_string()))?;
    write_df(df, &mut file, file_type)
}

fn write_df(
    df: &mut DataFrame,
    file: &mut std::fs::File,
    file_type: DataFrameFileType,
) -> ArpeggiaResult<()> {
    match file_type {
        DataFrameFileType::Csv => {
            CsvWriter::new(&mut *file)
                .finish(df)
                .map_err(|error| ArpeggiaError::Io(error.into()))?;
        }
        DataFrameFileType::Parquet => {
            ParquetWriter::new(&mut *file)
                .finish(df)
                .map_err(|error| ArpeggiaError::Io(error.into()))?;
        }
        DataFrameFileType::Json => {
            JsonWriter::new(&mut *file)
                .with_json_format(JsonFormat::Json)
                .finish(df)
                .map_err(|error| ArpeggiaError::Io(error.into()))?;
        }
        DataFrameFileType::NDJson => {
            JsonWriter::new(&mut *file)
                .with_json_format(JsonFormat::JsonLines)
                .finish(df)
                .map_err(|error| ArpeggiaError::Io(error.into()))?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_filename_cannot_escape_the_output_directory() {
        let output = std::env::temp_dir().join("arpeggia-output-path-test");
        assert!(prepare_df_output_dir(&output, "../escape", DataFrameFileType::Csv).is_err());
        assert!(prepare_df_output_dir(&output, "nested/name", DataFrameFileType::Csv).is_err());
        assert_eq!(
            prepare_df_output_dir(&output, "contacts", DataFrameFileType::Csv)
                .unwrap()
                .file_name(),
            Some(std::ffi::OsStr::new("contacts.csv"))
        );
    }

    #[test]
    fn new_dataframe_output_does_not_overwrite() {
        let path = std::env::temp_dir().join(format!("arpeggia-no-clobber-{}", std::process::id()));
        let output = path.with_extension("csv");
        std::fs::write(&output, "existing").unwrap();
        let mut dataframe = df!("value" => [1]).unwrap();
        assert!(write_df_to_new_file(&mut dataframe, &path, DataFrameFileType::Csv).is_err());
        assert_eq!(std::fs::read_to_string(output).unwrap(), "existing");
    }
}

/// File format for writing `DataFrames`.
#[derive(clap::ValueEnum, Clone, Debug, Copy)]
pub enum DataFrameFileType {
    /// Comma-separated values
    Csv,
    /// Parquet columnar storage
    Parquet,
    /// Standard JSON
    Json,
    /// Newline-delimited JSON
    #[value(name = "ndjson")]
    NDJson,
}

impl std::fmt::Display for DataFrameFileType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            DataFrameFileType::Csv => write!(f, "csv"),
            DataFrameFileType::Parquet => write!(f, "parquet"),
            DataFrameFileType::Json => write!(f, "json"),
            DataFrameFileType::NDJson => write!(f, "ndjson"),
        }
    }
}
