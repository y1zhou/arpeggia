#![warn(missing_docs)]
#![allow(clippy::type_complexity)]
#![doc = include_str!("../README.md")]

mod chains;
mod core;
mod interactions;
mod residues;
mod utils;

use crate::chains::ChainExt;

use clap::Parser;
use polars::prelude::*;
use std::path::{Path, PathBuf};
use tracing::{debug, error, info, trace, warn};

#[derive(Parser, Debug)]
#[command(version, about)]
struct Args {
    /// Path to the PDB or mmCIF file to be analyzed
    #[arg(short, long)]
    input: PathBuf,

    /// Output directory
    #[arg(short, long)]
    output: PathBuf,

    /// Group chains for interactions:
    /// e.g. A,B/C,D
    /// where chains A and B are the "ligand" and C and D are the "receptor"
    #[arg(short, long)]
    groups: String,

    /// Compensation factor for VdW radii dependent interaction types
    #[arg(short = 'c', long = "vdw-comp", default_value_t = 0.1)]
    vdw_comp: f64,

    /// Distance cutoff when searching for neighboring atoms
    #[arg(short, long, default_value_t = 4.5)]
    dist_cutoff: f64,

    /// Number of threads to use for parallel processing
    #[arg(short = 'j', long = "num-threads", default_value_t = 0)]
    num_threads: usize,
}

/// Entry to the CLI tool. Verbosity can be controlled with the `RUST_LOG` environment variable.
fn main() {
    let args = Args::parse();
    tracing_subscriber::fmt::init();
    trace!("{args:?}");

    // Create Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.num_threads)
        .build_global()
        .unwrap();
    debug!("Using {} thread(s)", rayon::current_num_threads());

    let run_result = core::core(
        args.input,
        args.groups.as_str(),
        args.vdw_comp,
        args.dist_cutoff,
    );
    match run_result {
        Err(e) => {
            e.iter().for_each(|e| match e.level() {
                pdbtbx::ErrorLevel::BreakingError => error!("{e}"),
                pdbtbx::ErrorLevel::InvalidatingError => error!("{e}"),
                _ => warn!("{e}"),
            });
            std::process::exit(1);
        }
        Ok((mut df_atomic, mut df_ring, i_complex, pdb_warnings)) => {
            // Notify of any PDB warnings
            if !pdb_warnings.is_empty() {
                pdb_warnings.iter().for_each(|e| warn!("{e}"));
            }

            // Information on the sequence of the chains in the model
            info!(
                "Loaded {} {}",
                i_complex.model.chain_count(),
                match i_complex.model.chain_count() {
                    1 => "chain",
                    _ => "chains",
                }
            );
            for chain in i_complex.model.chains() {
                debug!(">{}: {}", chain.id(), chain.pdb_seq().join(""));
            }

            debug!(
                "Parsed ligand chains {lig:?}; receptor chains {receptor:?}",
                lig = i_complex.ligand,
                receptor = i_complex.receptor
            );

            // Prepare output directory
            // let file_id = input_path.file_stem().unwrap().to_str().unwrap();
            let output_path = Path::new(&args.output);
            let _ = std::fs::create_dir_all(output_path);
            let output_dir = output_path.to_str().unwrap();

            debug!("Results will be saved to {output_dir}");

            // Save results and log the identified interactions
            config_polars_output();
            info!(
                "Found {} atom-atom contacts\n{}",
                df_atomic.shape().0,
                df_atomic
            );
            let df_clash = df_atomic
                .clone()
                .lazy()
                .filter(col("interaction").eq(lit("StericClash")))
                .collect()
                .unwrap();
            if df_clash.height() > 0 {
                warn!("Found {} steric clashes\n{}", df_clash.shape().0, df_clash);
            }
            info!("Found {} ring contacts\n{}", df_ring.shape().0, df_ring);

            // Save results to CSV files
            write_df_to_csv(&mut df_atomic, output_path.join("atomic_contacts.csv"));
            write_df_to_csv(&mut df_ring, output_path.join("ring_contacts.csv"));
        }
    }
}

fn config_polars_output() {
    std::env::set_var("POLARS_FMT_TABLE_HIDE_DATAFRAME_SHAPE_INFORMATION", "1");
    std::env::set_var("POLARS_FMT_TABLE_HIDE_COLUMN_DATA_TYPES", "1");
    std::env::set_var("POLARS_FMT_TABLE_ROUNDED_CORNERS", "1");
    std::env::set_var("POLARS_FMT_MAX_COLS", "14");
}

fn write_df_to_csv(df: &mut DataFrame, file_path: PathBuf) {
    let mut file = std::fs::File::create(file_path).unwrap();
    CsvWriter::new(&mut file).finish(df).unwrap();
}
