use arpeggia::{DataFrameFileType, get_num_threads, load_model, write_df_to_file};
use clap::{Parser, ValueEnum};
use std::path::{Path, PathBuf};
use tracing::{debug, error, info, trace, warn};

/// Granularity level for SAP score calculation.
#[derive(ValueEnum, Clone, Debug, Copy, Default)]
pub enum SapLevel {
    /// Calculate SAP score for each individual atom
    Atom,
    /// Aggregate SAP score by residue
    #[default]
    Residue,
}

impl std::fmt::Display for SapLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SapLevel::Atom => write!(f, "atom"),
            SapLevel::Residue => write!(f, "residue"),
        }
    }
}

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub(crate) struct Args {
    /// Path to the PDB or mmCIF file to be analyzed
    #[arg(short, long)]
    input: PathBuf,

    /// Output directory
    #[arg(short, long)]
    output: PathBuf,

    /// Name of the output file
    #[arg(short = 'f', long = "filename", default_value_t = String::from("sap"))]
    filename: String,

    /// Output file type
    #[arg(short = 't', long, default_value_t = DataFrameFileType::Csv)]
    output_format: DataFrameFileType,

    /// Model number to analyze (default: 0, the first model)
    #[arg(short = 'm', long = "model", default_value_t = 0)]
    model_num: usize,

    /// Probe radius r for SASA calculation (smaller r detects more surface details)
    #[arg(short = 'r', long = "probe-radius", default_value_t = 1.4)]
    probe_radius: f32,

    /// Number of points on the sphere for SASA sampling
    #[arg(short = 'n', long = "num-points", default_value_t = 100)]
    n_points: usize,

    /// Radius for SAP neighbor search in Ångströms
    #[arg(short = 's', long = "sap-radius", default_value_t = 5.0)]
    sap_radius: f32,

    /// Number of threads to use for parallel processing
    #[arg(short = 'j', long = "num-threads", default_value_t = 1)]
    num_threads: usize,

    /// Granularity level for SAP calculation
    #[arg(short = 'l', long = "level", default_value_t = SapLevel::Residue)]
    level: SapLevel,

    /// Comma-separated chain IDs to include (e.g., "H,L"). If empty, includes all chains.
    #[arg(short = 'c', long = "chains", default_value_t = String::new())]
    chains: String,
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

    // Calculate SAP score based on the specified level
    let mut df_sap = match args.level {
        SapLevel::Atom => arpeggia::get_per_atom_sap_score(
            &pdb,
            args.probe_radius,
            args.n_points,
            args.model_num,
            args.sap_radius,
            num_threads,
            &args.chains,
        ),
        SapLevel::Residue => arpeggia::get_per_residue_sap_score(
            &pdb,
            args.probe_radius,
            args.n_points,
            args.model_num,
            args.sap_radius,
            num_threads,
            &args.chains,
        ),
    };

    if df_sap.is_empty() {
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

    // Log SAP score summary
    let entity_name = match args.level {
        SapLevel::Atom => "atoms",
        SapLevel::Residue => "residues",
    };
    debug!(
        "Calculated SAP score for {} {}\n{}",
        df_sap.height(),
        entity_name,
        df_sap
    );

    // Save results to file
    write_df_to_file(&mut df_sap, &output_file, args.output_format);
    let output_file_str = output_file.to_str().unwrap();
    info!("Results saved to {output_file_str}");
}
