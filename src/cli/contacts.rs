use super::{DataFrameFileType, prepare_df_output_dir, write_df_to_file};
use arpeggia::{ArpeggiaResult, ProtonationMode, run_with_threads};
use clap::Parser;
use pdbtbx::*;
use polars::prelude::*;
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use tracing::{debug, info, trace, warn};

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

    /// How unresolved histidine protonation affects ionic-contact typing
    #[arg(long, value_enum, default_value_t = ProtonationMode::AllCharged)]
    protonation: ProtonationMode,

    /// pH used only by the heuristic protonation mode
    #[arg(long, default_value_t = 7.4)]
    ph: f64,
}

pub(crate) fn run(args: &Args) -> ArpeggiaResult<()> {
    trace!("{args:?}");

    // Make sure `input` exists
    let output_path = std::path::absolute(&args.output)?;

    // Load file as complex structure
    let mut pdb = super::load_input(Path::new(&args.input))?;
    let metadata = arpeggia::read_metadata(&args.input)?;
    for warning in metadata.warnings {
        warn!("{warning}");
    }
    let metadata = metadata.value;

    // Filter out atoms with zero occupancy if requested
    if args.ignore_zero_occupancy {
        pdb.remove_atoms_by(|atom| atom.occupancy() == 0.0);
        debug!("Removed atoms with zero occupancy");
    }

    if pdb
        .par_atoms()
        .filter(|atom| atom.element() == Some(&Element::H))
        .count()
        == 0
    {
        warn!(
            "No hydrogen atoms found in the structure. This may affect the accuracy of the results."
        );
    }

    // Use the library function
    let analysis = run_with_threads(args.num_threads as isize, || {
        debug!("Using {} thread(s)", rayon::current_num_threads());
        arpeggia::analyze_contacts(
            &pdb,
            Some(&metadata),
            args.groups.as_str(),
            args.vdw_comp,
            args.dist_cutoff,
            args.protonation,
            args.ph,
        )
    })?;
    for warning in analysis.warnings {
        warn!("{warning}");
    }
    let mut df_contacts = analysis.value;

    // Save results and log the identified interactions
    let clash_mask = df_contacts
        .column("interaction")
        .unwrap()
        .str()
        .unwrap()
        .equal("StericClash");
    let df_clash = df_contacts.filter(&clash_mask).unwrap();
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

    let output_file = prepare_df_output_dir(&output_path, &args.filename, args.output_format)?;
    write_df_to_file(&mut df_contacts, &output_file, args.output_format)?;
    info!("Results saved to {}", output_file.display());
    Ok(())
}
