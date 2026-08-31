use arpeggia::{
    ArpeggiaError, ArpeggiaResult, AtomSubset, ClusterOptions, ClusteringMethod, DataFrameFileType,
    PairwiseRmsdOptions, cluster_pairwise_rmsd, get_pairwise_rmsd_matrix, prepare_df_output_dir,
    read_pairwise_matrix, read_structure_observations, validate_residue_selection,
    write_df_to_file, write_df_to_new_file,
};
use clap::Parser;
use std::path::PathBuf;
use tracing::{debug, info, warn};

#[derive(Parser, Debug, Clone)]
#[command(version, about)]
pub(crate) struct Args {
    /// Non-recursive directory of PDB or mmCIF structures
    #[arg(short, long)]
    input: PathBuf,

    /// Output directory
    #[arg(short, long)]
    output: PathBuf,

    /// Name of the cluster table
    #[arg(short = 'f', long = "filename", default_value = "clusters")]
    filename: String,

    /// Output table format
    #[arg(short = 't', long, default_value_t = DataFrameFileType::Csv)]
    output_format: DataFrameFileType,

    /// Also save and reuse the pairwise RMSD table
    #[arg(long)]
    pairwise_rmsd: bool,

    /// Name of the optional pairwise RMSD table
    #[arg(long, default_value = "pairwise_rmsd")]
    pairwise_filename: String,

    /// Clustering algorithm
    #[arg(long, default_value = "k-medoids")]
    method: ClusteringMethod,

    /// Exact cluster count; takes priority over --max-clusters
    #[arg(long)]
    num_clusters: Option<usize>,

    /// Largest cluster count considered automatically
    #[arg(long)]
    max_clusters: Option<usize>,

    /// Maximum clustering iterations
    #[arg(long, default_value_t = 100)]
    max_iterations: usize,

    /// Model number to select (0 selects the first model)
    #[arg(short = 'm', long = "model", default_value_t = 0)]
    model_num: usize,

    /// Comma-separated chain and author-residue ranges
    #[arg(short = 'r', long, default_value_t = String::new())]
    residues: String,

    /// Atom population used for fitting and RMSD
    #[arg(short = 'a', long, default_value = "ca")]
    atoms: AtomSubset,

    /// Worker count (0 uses available processors)
    #[arg(short = 'j', long, default_value_t = 0)]
    num_threads: usize,

    /// Skip heuristic pairwise-matrix memory checks
    #[arg(long)]
    bypass_mem_check: bool,
}

pub(crate) fn run(args: &Args) -> ArpeggiaResult<()> {
    let cluster_options = ClusterOptions {
        method: args.method,
        num_clusters: args.num_clusters,
        max_clusters: args.max_clusters,
        max_iterations: args.max_iterations,
    };
    cluster_options.validate_without_structure_count()?;
    validate_residue_selection(&args.residues)?;

    let metadata = args.input.metadata().map_err(|error| {
        ArpeggiaError::Io(std::io::Error::new(
            error.kind(),
            format!("cannot inspect input {}: {error}", args.input.display()),
        ))
    })?;
    if !metadata.is_dir() {
        return Err(ArpeggiaError::InvalidArgument(
            "cluster-structs CLI input must be a directory".into(),
        ));
    }
    let observations = read_structure_observations(&args.input, "id", "path")?;

    cluster_options.validate(observations.len())?;

    let output_path = prepare_df_output_dir(&args.output, &args.filename, args.output_format)?;
    let pairwise_path = args
        .pairwise_rmsd
        .then(|| prepare_df_output_dir(&args.output, &args.pairwise_filename, args.output_format))
        .transpose()?;
    if pairwise_path.as_ref() == Some(&output_path) {
        return Err(ArpeggiaError::InvalidArgument(
            "cluster and pairwise RMSD output paths must differ".into(),
        ));
    }
    let pairwise_options = PairwiseRmsdOptions {
        model_num: args.model_num,
        residues: args.residues.clone(),
        atoms: args.atoms,
        num_threads: args.num_threads,
        bypass_mem_check: args.bypass_mem_check,
    };
    let matrix =
        if let Some(path) = pairwise_path.as_deref().filter(|path| path.exists()) {
            let expected_ids = observations
                .iter()
                .map(|observation| observation.id().to_string())
                .collect::<Vec<_>>();
            let analysis = read_pairwise_matrix(path, &expected_ids, args.bypass_mem_check)
                .map_err(|error| match error {
                    ArpeggiaError::Calculation(_) => error,
                    _ => ArpeggiaError::InvalidArgument(format!(
                        "invalid pairwise RMSD cache {}; remove it to recalculate: {error}",
                        path.display()
                    )),
                })?;
            for warning in analysis.warnings {
                warn!("{warning}");
            }
            debug!("Reusing pairwise RMSD cache {}", path.display());
            analysis.value
        } else {
            let analysis = get_pairwise_rmsd_matrix(&observations, &pairwise_options)?;
            for warning in analysis.warnings {
                warn!("{warning}");
            }
            if let Some(path) = &pairwise_path {
                let mut table = analysis.value.to_dataframe()?;
                write_df_to_new_file(&mut table, path, args.output_format)?;
                info!("Pairwise RMSD saved to {}", path.display());
            }
            analysis.value
        };

    let analysis = cluster_pairwise_rmsd(&matrix, &cluster_options)?;
    for warning in analysis.warnings {
        warn!("{warning}");
    }
    let mut clusters = analysis.value;
    write_df_to_file(&mut clusters, &output_path, args.output_format)?;
    info!("Clusters saved to {}", output_path.display());
    Ok(())
}
