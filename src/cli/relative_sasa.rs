use arpeggia::{
    ArpeggiaError, ArpeggiaResult, DataFrameFileType, prepare_df_output_dir, run_with_threads,
    write_df_to_file,
};
use clap::Parser;
use polars::prelude::{ChunkCompareEq, Column};
use std::path::{Path, PathBuf};
use tracing::{debug, error, info, trace};

#[derive(Parser, Debug, Clone)]
#[command(version, about = "Calculate relative SASA (RSA) for each residue")]
pub(crate) struct Args {
    /// Path to the PDB or mmCIF file to be analyzed
    #[arg(short, long)]
    input: PathBuf,

    /// Output directory
    #[arg(short, long)]
    output: PathBuf,

    /// Name of the output file
    #[arg(short = 'f', long = "filename", default_value_t = String::from("relative_sasa"))]
    filename: String,

    /// Output file type
    #[arg(short = 't', long, default_value_t = DataFrameFileType::Csv)]
    output_format: DataFrameFileType,

    /// Model number to analyze (default: 0, the first model)
    #[arg(short = 'm', long = "model", default_value_t = 0)]
    model_num: usize,

    /// Probe radius r (smaller r detects more surface details and reports a larger surface)
    #[arg(short = 'r', long = "probe-radius", default_value_t = 1.4)]
    probe_radius: f32,

    /// Number of points on the sphere for sampling
    #[arg(short = 'n', long = "num-points", default_value_t = 100)]
    n_points: usize,

    /// Number of threads to use for parallel processing
    #[arg(short = 'j', long = "num-threads", default_value_t = 1)]
    num_threads: usize,

    /// Comma-separated chain IDs to include (e.g., "A,B,C"). If empty, includes all chains.
    #[arg(short = 'c', long = "chains", default_value_t = String::new())]
    chains: String,
}

pub(crate) fn run(args: &Args) -> ArpeggiaResult<()> {
    trace!("{args:?}");

    // Make sure `input` exists
    let pdb = super::load_input(Path::new(&args.input))?;

    // Calculate relative SASA
    let analysis = run_with_threads(args.num_threads as isize, || {
        debug!("Using {} thread(s)", rayon::current_num_threads());
        arpeggia::get_relative_sasa(
            &pdb,
            args.probe_radius,
            args.n_points,
            args.model_num,
            &args.chains,
        )
    })?;
    for warning in analysis.warnings {
        tracing::warn!("{warning}");
    }
    let mut df_relative_sasa = analysis.value;

    if df_relative_sasa.height() == 0 {
        error!(
            "No data found in the input file. Please check the provided arguments, especially the model number."
        );
        return Ok(());
    }

    // Log summary statistics
    if tracing::enabled!(tracing::Level::DEBUG) {
        let non_zero_sasa_mask = df_relative_sasa
            .column("sasa")
            .and_then(Column::f32)
            .map_err(|error| ArpeggiaError::Calculation(error.to_string()))?
            .not_equal(0.0);
        let df_sasa_nonzero = df_relative_sasa
            .filter(&non_zero_sasa_mask)
            .map_err(|error| ArpeggiaError::Calculation(error.to_string()))?;
        debug!(
            "Found {} residues with non-zero SASA\n{}",
            df_sasa_nonzero.height(),
            df_sasa_nonzero
        );
    }

    // Save results to file
    let output_file = prepare_df_output_dir(&args.output, &args.filename, args.output_format)?;
    write_df_to_file(&mut df_relative_sasa, &output_file, args.output_format)?;
    info!("Results saved to {}", output_file.display());
    Ok(())
}
