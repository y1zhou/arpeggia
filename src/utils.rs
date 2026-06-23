use polars::prelude::*;
use std::path::{Path, PathBuf};

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
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build()
            .unwrap_or_else(|_| {
                // Fallback to global pool if creation fails
                rayon::ThreadPoolBuilder::new()
                    .build()
                    .expect("Failed to create thread pool")
            });
        pool.install(f)
    }
}

/// Sum the SASA column of a `DataFrame`.
///
/// # Arguments
///
/// * `df` - `DataFrame` containing a "sasa" column
///
/// # Returns
///
/// The sum of all SASA values, or 0.0 if the column is empty or doesn't exist.
pub fn sum_float_col(df: &DataFrame, colname: &str) -> f32 {
    df.column(colname)
        .unwrap()
        .f32()
        .unwrap()
        .sum()
        .unwrap_or(0.0)
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
) -> PathBuf {
    let output_path = match std::path::absolute(output) {
        Ok(path) => path,
        Err(e) => {
            panic!("Failed to resolve the output directory: {}", e);
        }
    };
    let _ = std::fs::create_dir_all(&output_path);

    if output_path.is_dir() {
        output_path.join(filename)
    } else {
        output_path.canonicalize().unwrap()
    }
    .with_extension(output_format.to_string())
}

/// Write a `DataFrame` to a CSV file
///
/// # Panics
/// This function will panic if the file cannot be created or written to.
pub fn write_df_to_file(df: &mut DataFrame, file_path: &Path, file_type: DataFrameFileType) {
    let file_suffix = file_type.to_string();
    let mut file = std::fs::File::create(file_path.with_extension(file_suffix)).unwrap();
    match file_type {
        DataFrameFileType::Csv => {
            CsvWriter::new(&mut file).finish(df).unwrap();
        }
        DataFrameFileType::Parquet => {
            ParquetWriter::new(&mut file).finish(df).unwrap();
        }
        DataFrameFileType::Json => {
            JsonWriter::new(&mut file)
                .with_json_format(JsonFormat::Json)
                .finish(df)
                .unwrap();
        }
        DataFrameFileType::NDJson => {
            JsonWriter::new(&mut file)
                .with_json_format(JsonFormat::JsonLines)
                .finish(df)
                .unwrap();
        }
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
