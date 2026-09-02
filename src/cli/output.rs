use arpeggia::{ArpeggiaError, ArpeggiaResult};
use clap::ValueEnum;
use polars::prelude::*;
use std::path::{Path, PathBuf};

#[derive(ValueEnum, Clone, Debug, Copy)]
pub(super) enum DataFrameFileType {
    Csv,
    Parquet,
    #[value(name = "ndjson")]
    NDJson,
}

impl std::fmt::Display for DataFrameFileType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str(match self {
            Self::Csv => "csv",
            Self::Parquet => "parquet",
            Self::NDJson => "ndjson",
        })
    }
}

pub(super) fn prepare_df_output_dir(
    output: &Path,
    filename: &str,
    output_format: DataFrameFileType,
) -> ArpeggiaResult<PathBuf> {
    let mut components = Path::new(filename).components();
    if !matches!(components.next(), Some(std::path::Component::Normal(_)))
        || components.next().is_some()
    {
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

pub(super) fn write_df_to_file(
    dataframe: &mut DataFrame,
    path: &Path,
    file_type: DataFrameFileType,
) -> ArpeggiaResult<()> {
    let mut file = std::fs::File::create(path.with_extension(file_type.to_string()))?;
    write_df(dataframe, &mut file, file_type)
}

pub(super) fn write_df_to_new_file(
    dataframe: &mut DataFrame,
    path: &Path,
    file_type: DataFrameFileType,
) -> ArpeggiaResult<()> {
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path.with_extension(file_type.to_string()))?;
    write_df(dataframe, &mut file, file_type)
}

fn write_df(
    dataframe: &mut DataFrame,
    file: &mut std::fs::File,
    file_type: DataFrameFileType,
) -> ArpeggiaResult<()> {
    match file_type {
        DataFrameFileType::Csv => CsvWriter::new(file),
        DataFrameFileType::Parquet => {
            ParquetWriter::new(file)
                .finish(dataframe)
                .map_err(|error| ArpeggiaError::Io(error.into()))?;
            return Ok(());
        }
        DataFrameFileType::NDJson => {
            JsonWriter::new(file)
                .with_json_format(JsonFormat::JsonLines)
                .finish(dataframe)
                .map_err(|error| ArpeggiaError::Io(error.into()))?;
            return Ok(());
        }
    }
    .finish(dataframe)
    .map_err(|error| ArpeggiaError::Io(error.into()))?;
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
