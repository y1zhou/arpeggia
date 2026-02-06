use arpeggia::{DataFrameFileType, load_model, run_with_threads, write_df_to_file};
use clap::Parser;
use pdbtbx::*;
use polars::prelude::*;
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use tracing::{debug, error, info, trace, warn};

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub(crate) struct Args {
    /// Path to the PDB or mmCIF file to be analyzed
    #[arg(short, long)]
    input: PathBuf,

    /// Output directory
    #[arg(short, long)]
    output: PathBuf,

    /// Group chains for interactions:
    /// e.g. A,B/C,D
    /// where chains A and B are the "ligand" and C and D are the "receptor".
    /// Chains can exist on both sides, in which case intra-chain interactions will be calculated.
    /// If only one group is provided, all remaining chains will be considered as the other group.
    /// If no groups are provided ('/'), all inter- and intra-chain interactions will be calculated.
    #[arg(short, long, default_value_t = String::from("/"))]
    groups: String,

    /// Name of the output file
    #[arg(short = 'f', long = "filename", default_value_t = String::from("contacts"))]
    filename: String,

    /// Output file type
    #[arg(short = 't', long, default_value_t = DataFrameFileType::Csv)]
    output_format: DataFrameFileType,

    /// Compensation factor for VdW radii dependent interaction types
    #[arg(short = 'c', long = "vdw-comp", default_value_t = 0.1)]
    vdw_comp: f64,

    /// Distance cutoff when searching for neighboring atoms
    #[arg(short, long, default_value_t = 6.5)]
    dist_cutoff: f64,

    /// Number of threads to use for parallel processing. One thread should be sufficient unless the system is very large
    #[arg(short = 'j', long = "num-threads", default_value_t = 1)]
    num_threads: usize,

    /// Ignore atoms with zero occupancy
    #[arg(long = "ignore-zero-occupancy", default_value_t = false)]
    ignore_zero_occupancy: bool,
}

pub(crate) fn run(args: &Args) {
    trace!("{args:?}");

    // Make sure `input` exists
    let input_path = match Path::new(&args.input).canonicalize() {
        Ok(path) => path,
        Err(e) => {
            error!("Failed to retrieve input file: {}", e);
            return;
        }
    };
    let output_path = match std::path::absolute(&args.output) {
        Ok(path) => path,
        Err(e) => {
            error!("Failed to resolve the output directory: {}", e);
            return;
        }
    };
    let input_file: String = input_path.to_str().unwrap().parse().unwrap();

    // Load file as complex structure
    let (mut pdb, pdb_warnings) = load_model(&input_file);
    if !pdb_warnings.is_empty() {
        for e in &pdb_warnings {
            match e.level() {
                pdbtbx::ErrorLevel::BreakingError => error!("{e}"),
                pdbtbx::ErrorLevel::InvalidatingError => error!("{e}"),
                _ => warn!("{e}"),
            }
        }
    }

    // Filter out atoms with zero occupancy if requested
    if args.ignore_zero_occupancy {
        pdb.remove_atoms_by(|atom| atom.occupancy() == 0.0);
        debug!("Removed atoms with zero occupancy");
    }

    if pdb
        .par_atoms()
        .filter(|a| a.element().unwrap() == &Element::H)
        .count()
        == 0
    {
        warn!(
            "No hydrogen atoms found in the structure. This may affect the accuracy of the results."
        );
    }

    // Use the library function
    let mut df_contacts = run_with_threads(args.num_threads as isize, || {
        debug!("Using {} thread(s)", rayon::current_num_threads());
        arpeggia::get_contacts(&pdb, args.groups.as_str(), args.vdw_comp, args.dist_cutoff)
    });

    // Prepare output directory
    let _ = std::fs::create_dir_all(output_path.clone());
    let output_file = output_path
        .join(args.filename.clone())
        .with_extension(args.output_format.to_string());

    // Save results and log the identified interactions
    let df_clash = df_contacts
        .clone()
        .lazy()
        .filter(col("interaction").eq(lit("StericClash")))
        .collect()
        .unwrap();
    if df_clash.height() > 0 {
        warn!(
            "Found {} steric {}\n{}",
            df_clash.height(),
            match df_clash.height() {
                1 => "clash",
                _ => "clashes",
            },
            df_clash
        );
    }

    // Save res to CSV files
    write_df_to_file(&mut df_contacts, &output_file, args.output_format);
    let output_file_str = output_file.to_str().unwrap();
    info!("Results saved to {output_file_str}");
}
