use arpeggia::{DataFrameFileType, load_model, get_num_threads, write_df_to_file};
use clap::Parser;
use polars::prelude::ChunkCompareEq;
use std::path::{Path, PathBuf};
use tracing::{debug, error, info, trace, warn};

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
}

pub(crate) fn run(args: &Args) {
    trace!("{args:?}");

    // Create Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.num_threads)
        .build_global()
        .unwrap();
    debug!("Using {} thread(s)", rayon::current_num_threads());

    // Make sure `input` exists
    let input_path = Path::new(&args.input).canonicalize().unwrap();
    let input_file: String = input_path.to_str().unwrap().parse().unwrap();

    // Load file as complex structure
    let (pdb, pdb_warnings) = load_model(&input_file);
    if !pdb_warnings.is_empty() {
        pdb_warnings.iter().for_each(|e| match e.level() {
            pdbtbx::ErrorLevel::BreakingError => error!("{e}"),
            pdbtbx::ErrorLevel::InvalidatingError => error!("{e}"),
            _ => warn!("{e}"),
        });
    }

    // Convert thread count to isize for rust-sasa
    let num_threads = get_num_threads(args.num_threads);

    // Calculate relative SASA
    let mut df_relative_sasa = arpeggia::get_relative_sasa(
        &pdb,
        args.probe_radius,
        args.n_points,
        args.model_num,
        num_threads,
    );

    if df_relative_sasa.is_empty() {
        error!(
            "No data found in the input file. Please check the provided arguments, especially the model number."
        );
        return;
    }

    // Prepare output directory
    let output_path = Path::new(&args.output).canonicalize().unwrap();
    let _ = std::fs::create_dir_all(output_path.clone());
    let output_file = match output_path.is_dir() {
        true => output_path.join(args.filename.clone()),
        false => output_path,
    }
    .with_extension(args.output_format.to_string());

    // Log summary statistics
    let non_zero_sasa_mask = df_relative_sasa
        .column("sasa")
        .unwrap()
        .as_materialized_series()
        .not_equal(0.0)
        .unwrap();
    let df_sasa_nonzero = df_relative_sasa.filter(&non_zero_sasa_mask).unwrap();
    debug!(
        "Found {} residues with non-zero SASA\n{}",
        df_sasa_nonzero.height(),
        df_sasa_nonzero
    );

    // Save results to file
    write_df_to_file(&mut df_relative_sasa, &output_file, args.output_format);
    let output_file_str = output_file.to_str().unwrap();
    info!("Results saved to {output_file_str}");
}
