#![warn(missing_docs)]
#![doc = include_str!("../README.md")]

mod cli;

use arpeggia::ArpeggiaResult;
use clap::{Parser, Subcommand};
use std::process::ExitCode;

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
    /// Calculate the buried surface area (dSASA) at the interface between chain groups
    Dsasa(crate::cli::dsasa::Args),
    /// Calculate relative SASA (RSA) for each residue, normalized by Tien et al. (2013) MaxASA values
    RelativeSasa(crate::cli::relative_sasa::Args),
    /// Calculate Spatial Aggregation Propensity (SAP) score for aggregation-prone region prediction
    Sap(crate::cli::sap::Args),
    /// Calculate Shape Complementarity (SC) between two chain groups
    Sc(crate::cli::sc::Args),
    /// Print the sequences of all chains in a PDB or mmCIF file
    Seq(crate::cli::pdb2seq::Args),
    /// Print declared SEQRES/entity-polymer sequences
    Seqres(crate::cli::pdb2seq::Args),
}

/// Entry to the CLI tool. Verbosity can be controlled with the `RUST_LOG` environment variable.
fn main() -> ExitCode {
    tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
        .init();
    config_polars_output();

    let cli = Cli::parse();
    let result: ArpeggiaResult<()> = match &cli.command {
        Commands::Contacts(args) => crate::cli::contacts::run(args),
        Commands::Sasa(args) => crate::cli::sasa::run(args),
        Commands::Dsasa(args) => crate::cli::dsasa::run(args),
        Commands::RelativeSasa(args) => crate::cli::relative_sasa::run(args),
        Commands::Sap(args) => crate::cli::sap::run(args),
        Commands::Sc(args) => crate::cli::sc::run(args),
        Commands::Seq(args) => crate::cli::pdb2seq::run(args),
        Commands::Seqres(args) => crate::cli::pdb2seq::run_declared(args),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            tracing::error!("{error}");
            ExitCode::FAILURE
        }
    }
}

fn config_polars_output() {
    unsafe {
        std::env::set_var("POLARS_FMT_TABLE_HIDE_DATAFRAME_SHAPE_INFORMATION", "1");
        std::env::set_var("POLARS_FMT_TABLE_HIDE_COLUMN_DATA_TYPES", "1");
        std::env::set_var("POLARS_FMT_TABLE_ROUNDED_CORNERS", "1");
        std::env::set_var("POLARS_FMT_MAX_COLS", "14");
    }
}
