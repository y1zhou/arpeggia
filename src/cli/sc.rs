use arpeggia::load_model;
use clap::Parser;
use std::path::{Path, PathBuf};
use tracing::{debug, error, info, trace, warn};

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub(crate) struct Args {
    /// Path to the PDB or mmCIF file to be analyzed
    #[arg(short, long)]
    input: PathBuf,

    /// Chain groups for SC calculation:
    /// e.g. A,B/C,D where chains A and B form one surface and C and D form the other.
    /// Both groups must be specified for SC calculation.
    #[arg(short, long)]
    groups: String,

    /// Model number to analyze (default: 0, the first model)
    #[arg(short = 'm', long = "model", default_value_t = 0)]
    model_num: usize,

    /// Probe radius for surface generation (default: 1.7 Å, Connolly 1983)
    #[arg(short = 'r', long = "probe-radius", default_value_t = 1.7)]
    probe_radius: f64,

    /// Dot density for surface sampling (dots per Ų)
    #[arg(short = 'd', long = "density", default_value_t = 15.0)]
    density: f64,

    /// Peripheral band width for trimming (Å)
    #[arg(short = 'b', long = "band", default_value_t = 1.5)]
    band: f64,

    /// Separation cutoff for interface atom classification (Å)
    #[arg(short = 's', long = "separation", default_value_t = 8.0)]
    separation: f64,

    /// Gaussian weight for distance-weighted SC calculation (Å^-2)
    #[arg(short = 'w', long = "weight", default_value_t = 0.5)]
    weight: f64,

    /// Output format: 'text' for human-readable or 'json' for machine-readable
    #[arg(short = 'f', long = "format", default_value_t = String::from("text"))]
    output_format: String,
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
    let input_file: String = input_path.to_str().unwrap().parse().unwrap();

    // Validate groups
    if !args.groups.contains('/') {
        error!("Groups must be specified as 'A,B/C,D' with both surfaces defined");
        return;
    }

    // Load file as complex structure
    let (pdb, pdb_warnings) = load_model(&input_file);
    if !pdb_warnings.is_empty() {
        pdb_warnings.iter().for_each(|e| match e.level() {
            pdbtbx::ErrorLevel::BreakingError => error!("{e}"),
            pdbtbx::ErrorLevel::InvalidatingError => error!("{e}"),
            _ => warn!("{e}"),
        });
    }

    // Create settings
    let settings = arpeggia::ScSettings {
        probe_radius: args.probe_radius,
        density: args.density,
        band: args.band,
        separation: args.separation,
        weight: args.weight,
    };

    // Calculate SC
    let result = arpeggia::get_sc(&pdb, &args.groups, args.model_num, Some(settings));

    if !result.valid {
        error!("SC calculation failed. Check that both chain groups have atoms at the interface.");
        return;
    }

    // Output results
    if args.output_format == "json" {
        println!(
            r#"{{"sc": {:.6}, "distance": {:.4}, "area": {:.2}, "surface1": {{"atoms": {}, "buried_atoms": {}, "dots": {}, "trimmed_dots": {}, "trimmed_area": {:.2}, "s_median": {:.4}, "d_median": {:.4}}}, "surface2": {{"atoms": {}, "buried_atoms": {}, "dots": {}, "trimmed_dots": {}, "trimmed_area": {:.2}, "s_median": {:.4}, "d_median": {:.4}}} }}"#,
            result.sc,
            result.distance,
            result.area,
            result.surfaces[0].n_atoms,
            result.surfaces[0].n_buried_atoms,
            result.surfaces[0].n_dots,
            result.surfaces[0].n_trimmed_dots,
            result.surfaces[0].trimmed_area,
            result.surfaces[0].s_median,
            result.surfaces[0].d_median,
            result.surfaces[1].n_atoms,
            result.surfaces[1].n_buried_atoms,
            result.surfaces[1].n_dots,
            result.surfaces[1].n_trimmed_dots,
            result.surfaces[1].trimmed_area,
            result.surfaces[1].s_median,
            result.surfaces[1].d_median,
        );
    } else {
        info!("Shape Complementarity Results:");
        info!("  SC Score: {:.4}", result.sc);
        info!("  Median Distance: {:.4} Å", result.distance);
        info!("  Total Trimmed Area: {:.2} Ų", result.area);
        debug!("Surface 1:");
        debug!("  Atoms: {} ({} at interface)", result.surfaces[0].n_atoms, result.surfaces[0].n_buried_atoms);
        debug!("  Dots: {} total, {} trimmed", result.surfaces[0].n_dots, result.surfaces[0].n_trimmed_dots);
        debug!("  Trimmed Area: {:.2} Ų", result.surfaces[0].trimmed_area);
        debug!("  SC Median: {:.4}, Distance Median: {:.4} Å", result.surfaces[0].s_median, result.surfaces[0].d_median);
        debug!("Surface 2:");
        debug!("  Atoms: {} ({} at interface)", result.surfaces[1].n_atoms, result.surfaces[1].n_buried_atoms);
        debug!("  Dots: {} total, {} trimmed", result.surfaces[1].n_dots, result.surfaces[1].n_trimmed_dots);
        debug!("  Trimmed Area: {:.2} Ų", result.surfaces[1].trimmed_area);
        debug!("  SC Median: {:.4}, Distance Median: {:.4} Å", result.surfaces[1].s_median, result.surfaces[1].d_median);
    }
}
