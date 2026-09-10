use arpeggia::{ArpeggiaResult, get_dsasa_components, run_with_threads};
use clap::Parser;
use std::path::{Path, PathBuf};
use tracing::{debug, info, trace};

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub(crate) struct Args {
    /// Path to the PDB or mmCIF file to be analyzed
    #[arg(short, long)]
    input: PathBuf,

    /// Disjoint, non-empty interface groups, e.g. A,B/C,D or A/ (A vs remaining chains).
    /// Two-sided area in Å² = SASA(group1) + SASA(group2) - SASA(complex).
    #[arg(short, long)]
    groups: String,

    /// Model serial to analyze (0 selects the first model)
    #[arg(short = 'm', long = "model", default_value_t = 0)]
    model_num: usize,

    /// Solvent probe radius in Ångströms. Smaller probes access narrower crevices;
    /// larger probes exclude them. Total SASA changes depend on the structure.
    #[arg(short = 'r', long = "probe-radius", default_value_t = 1.4)]
    probe_radius: f32,

    /// Positive number of sample points per atomic sphere
    #[arg(short = 'n', long = "num-points", default_value_t = 100)]
    n_points: usize,

    /// Worker count (0 uses available processors)
    #[arg(short = 'j', long = "num-threads", default_value_t = 1)]
    num_threads: usize,
}

pub(crate) fn run(args: &Args) -> ArpeggiaResult<()> {
    trace!("{args:?}");

    // Make sure `input` exists
    let pdb = super::load_input(Path::new(&args.input))?;

    // Convert thread count to isize for rust-sasa
    let dsasa = run_with_threads(args.num_threads as isize, || {
        debug!("Using {} thread(s)", rayon::current_num_threads());
        get_dsasa_components(
            &pdb,
            &args.groups,
            args.probe_radius,
            args.n_points,
            args.model_num,
        )
    })?;
    for warning in dsasa.warnings {
        tracing::warn!("{warning}");
    }
    let dsasa = dsasa.value;
    info!(
        "Buried surface area (dSASA) at the interface between chains [{}]: {:.2} Å²",
        args.groups, dsasa.dsasa
    );
    info!(
        "polar={:.2} Å² hydrophobic={:.2} Å² unclassified={:.2} Å²",
        dsasa.polar_dsasa, dsasa.hydrophobic_dsasa, dsasa.unclassified_dsasa
    );
    Ok(())
}
