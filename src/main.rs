#![warn(missing_docs)]
#![allow(clippy::type_complexity)]
#![doc = include_str!("../README.md")]

mod chains;
mod cli;
mod interactions;
mod residues;
mod utils;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(version, about, author)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Clone)]
enum Commands {
    /// Analyze atomic and ring contacts in a PDB or mmCIF file
    Contacts(crate::cli::contacts::Args),
    /// Calculate the solvent accessible surface area (SASA) of each atom in a PDB or mmCIF file
    Sasa(crate::cli::sasa::Args),
}

/// Entry to the CLI tool. Verbosity can be controlled with the `RUST_LOG` environment variable.
fn main() {
    tracing_subscriber::fmt::init();
    config_polars_output();

    let cli = Cli::parse();
    match &cli.command {
        Commands::Contacts(args) => {
            crate::cli::contacts::run(args);
        }
        Commands::Sasa(args) => {
            crate::cli::sasa::run(args);
        }
    }
}

fn config_polars_output() {
    std::env::set_var("POLARS_FMT_TABLE_HIDE_DATAFRAME_SHAPE_INFORMATION", "1");
    std::env::set_var("POLARS_FMT_TABLE_HIDE_COLUMN_DATA_TYPES", "1");
    std::env::set_var("POLARS_FMT_TABLE_ROUNDED_CORNERS", "1");
    std::env::set_var("POLARS_FMT_MAX_COLS", "14");
}
