//! Python bindings for the arpeggia library using PyO3.
//!
//! This module provides Python-friendly wrappers around the core Rust functions,
//! converting Polars DataFrames to Python using pyo3-polars for efficient zero-copy data transfer.

use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use std::ffi::CString;

mod ab_numbering;

fn value_enum<T: clap::ValueEnum>(value: &str, error: &'static str) -> PyResult<T> {
    T::from_str(&value.replace('_', "-"), true)
        .map_err(|_| pyo3::exceptions::PyValueError::new_err(error))
}

fn python_error(error: crate::ArpeggiaError) -> PyErr {
    match error {
        crate::ArpeggiaError::Io(error) => pyo3::exceptions::PyOSError::new_err(error.to_string()),
        crate::ArpeggiaError::Parse(message) | crate::ArpeggiaError::InvalidArgument(message) => {
            pyo3::exceptions::PyValueError::new_err(message)
        }
        crate::ArpeggiaError::Calculation(message) => {
            pyo3::exceptions::PyRuntimeError::new_err(message)
        }
    }
}

fn load_for_python(py: Python<'_>, input_file: &str) -> PyResult<pdbtbx::PDB> {
    let analysis = py
        .detach(|| crate::load_model(input_file))
        .map_err(python_error)?;
    emit_python_warnings(py, analysis.warnings)?;
    Ok(analysis.value)
}

fn emit_python_warnings(
    py: Python<'_>,
    warnings: impl IntoIterator<Item = impl std::fmt::Display>,
) -> PyResult<()> {
    let category = py.get_type::<pyo3::exceptions::PyUserWarning>();
    for warning in warnings {
        let message = CString::new(warning.to_string())
            .map_err(|_| pyo3::exceptions::PyValueError::new_err("warning contains a NUL byte"))?;
        PyErr::warn(py, &category, &message, 1)?;
    }
    Ok(())
}

fn protonation_mode(value: &str) -> PyResult<crate::ProtonationMode> {
    value_enum(
        value,
        "protonation must be 'all-charged', 'heuristic', or 'explicit-only'",
    )
}

fn atom_subset(value: &str) -> PyResult<crate::AtomSubset> {
    value_enum(value, "atoms must be 'ca', 'backbone', 'heavy', or 'all'")
}

fn clustering_method(value: &str) -> PyResult<crate::ClusteringMethod> {
    value_enum(value, "method must be 'k-medoids'")
}

fn display_options(py: Python<'_>, width: Option<usize>, color: &str) -> PyResult<(usize, bool)> {
    let color: crate::AlignmentColor =
        value_enum(color, "color must be 'auto', 'always', or 'never'")?;
    let stdout = py.import("sys")?.getattr("stdout")?;
    let terminal = stdout
        .call_method0("isatty")
        .and_then(|v| v.extract::<bool>())
        .unwrap_or(false);
    let width = match width {
        Some(width) => width,
        None => py
            .import("shutil")?
            .call_method0("get_terminal_size")?
            .getattr("columns")?
            .extract()?,
    };
    Ok((width, crate::seq_alignment::color_enabled(color, terminal)))
}

#[pymethods]
impl crate::SeqAlignment {
    /// Format statistics and three alignment rows, optionally with position rulers.
    ///
    /// Args:
    ///     width (int | None): Total visible columns including labels. None
    ///         detects terminal width or uses 80 when unavailable.
    ///     color (str): "auto", "always", or "never". Auto (default) inspects
    ///         sys.stdout and respects NO_COLOR.
    ///     rulers (bool): Show position rulers by default. False hides ticks
    ///         but retains start/end numbers.
    ///
    /// Returns:
    ///     str: Statistics and wrapped alignment text. Color escapes appear
    ///         only when enabled; stored data fields always remain plain.
    ///
    /// Raises:
    ///     ValueError: The width cannot fit labels and one residue, or the
    ///         color policy is unsupported.
    #[pyo3(name = "format", signature = (width=None, color="auto", rulers=true))]
    fn format_python(
        &self,
        py: Python<'_>,
        width: Option<usize>,
        color: &str,
        rulers: bool,
    ) -> PyResult<String> {
        let (width, color) = display_options(py, width, color)?;
        self.render(width, color, rulers).map_err(python_error)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.format_python(py, None, "auto", true)
    }
}

#[pymethods]
impl crate::RmsdResult {
    fn __repr__(&self) -> String {
        format!(
            "RmsdResult(rmsd={:.6}, core_rmsd={:.6}, fit_atoms={}/{}, evaluation_atoms={}, cycles={}, chains={})",
            self.rmsd,
            self.core_rmsd,
            self.retained_fit_atoms,
            self.initial_fit_atoms,
            self.evaluation_atoms,
            self.cycles,
            self.chain_alignments.len()
        )
    }
}

/// Align two unaligned amino-acid strings and return a read-only SeqAlignment.
///
/// Printing the result shows colored, wrapped sequences and position rulers when
/// the terminal permits. Use result.format(width=60, color="never", rulers=False)
/// for plain output without rulers. Data fields contain no color escapes.
///
/// Args:
///     reference (str): Reference input sequence. Lowercase is
///         normalized; standard amino acids and B/Z/X/U/O are accepted. Inputs
///         must be non-empty and contain no whitespace, gaps, or stop symbols.
///     query (str): Query input sequence, with the same alphabet restrictions.
///     mode (str): "global" (default) consumes both inputs; "local" selects
///         best-scoring subsequences; "semi-global" consumes the entire query
///         sequence with free reference terminal overhangs.
///     gap_open (float): Opening cost ≥ 0.01, default 10, with at most two
///         decimals and no smaller than gap_extend. A gap of length L costs
///         gap_open + (L - 1) * gap_extend.
///     gap_extend (float): Cost per additional gap residue ≥ 0.01, default
///         0.5, with at most two decimals and no larger than gap_open.
///     reference_name (str): Keyword-only reference display label, default
///         "Reference". Stored on the result and in JSON; does not affect scoring.
///     query_name (str): Keyword-only query display label, default "Query".
///
/// Returns:
///     SeqAlignment: Plain gapped aligned_reference/aligned_query strings,
///         reference-to-query operations (space for match, + insertion, - deletion,
///         : positive-score substitution, x other substitution), zero-based half-open
///         input spans, score, matches, mismatches, gap_residues, and gap_runs.
///         identity_alignment/identity_shorter and coverage_alignment/coverage_shorter
///         are ratios using alignment-column and shorter-full-input denominators.
///         edit_distance is full-input Levenshtein distance, independent of score.
///         Empty local results have score zero and None alignment-length ratios.
///         U/O score as C/K with a warning but retain distinct identity.
///         Positive-score substitutions are blue in displays and remain mismatches.
///
/// Examples:
///     >>> alignment = arpeggia.align_seqs("GGACDEFGHIKGG", "ACDEFGHIK", mode="semi-global")
///     >>> alignment.reference_span
///     (2, 11)
#[pyfunction]
#[pyo3(signature = (reference, query, mode="global", gap_open=10.0, gap_extend=0.5, *, reference_name="Reference", query_name="Query"))]
#[allow(clippy::too_many_arguments)]
fn align_seqs(
    py: Python<'_>,
    reference: &str,
    query: &str,
    mode: &str,
    gap_open: f64,
    gap_extend: f64,
    reference_name: &str,
    query_name: &str,
) -> PyResult<crate::SeqAlignment> {
    let options = crate::SeqAlignOptions {
        mode: value_enum(mode, "mode must be 'global', 'local', or 'semi-global'")?,
        gap_open,
        gap_extend,
    };
    let analysis = py
        .detach(|| crate::align_seqs(reference, query, &options))
        .map_err(python_error)?;
    emit_python_warnings(py, analysis.warnings)?;
    Ok(crate::SeqAlignment {
        reference_name: reference_name.into(),
        query_name: query_name.into(),
        ..analysis.value
    })
}

/// Superpose two PDB/mmCIF structures; return a read-only RmsdResult in Ångströms.
///
/// Args:
///     reference (str): Path to the reference PDB/mmCIF structure.
///     query (str): Path to the query PDB/mmCIF structure to superpose.
///     model_num (int): Model serial; 0 independently selects each first model.
///     superpose_residues (str): Fit selection on the reference structure, using
///         reference chain IDs and author residue numbers. With align_seqs=True,
///         corresponding query residues come from the sequence alignment. With
///         align_seqs=False, apply the same selection to the query and require
///         exact atom identities. Empty selects all eligible residues.
///         Use "A" for a whole chain or "A:1-100,A:110-120,B"
///         for a union of inclusive author-number ranges and chains. Repeat the
///         chain in each clause. Negative numbers ("A:-5--1") and insertion
///         codes ("B:10A-20") are accepted; "A:10" includes all insertion variants.
///     rmsd_residues (str): Evaluation selection on the reference structure,
///         using the same syntax and query-correspondence rules as the fit
///         selection. Empty independently means all; it does not inherit
///         superpose_residues. Evaluation never determines or refits the transform.
///     atoms (str): "ca" (default), "backbone" (N/CA/C/O/OXT), "heavy", or "all".
///     align_seqs (bool): Align complete observed chains before atom selection.
///         When True, both residue selectors use reference author numbering;
///         otherwise selections apply to both structures with exact correspondence.
///     chain_map (dict[str, str] | None): Explicit reference-to-query chain IDs,
///         e.g. {"A": "H"}. Requires align_seqs=True, covers exactly selected
///         reference chains, and assigns unique query chains. None infers a
///         unique maximum-score assignment; ambiguous homomers require a map.
///     alignment_mode (str): "global" (default), "local", or "semi-global" for
///         final chain alignment. Semi-global consumes the complete query chain.
///         Chain inference independently scores shorter against longer semi-globally.
///     gap_open (float): Alignment opening cost ≥ 0.01, default 10; at most
///         two decimals and no smaller than gap_extend. Gap cost is open+(L-1)*extend.
///     gap_extend (float): Cost per additional gap residue ≥ 0.01, default 0.5;
///         at most two decimals and no larger than gap_open.
///     refine_cycles (int): Maximum rejection/refit rounds after the initial fit;
///         0 (default) performs no rejection. Stops early when unchanged or exact.
///     refine_cutoff (float): Positive rejection multiplier (default 2) of current
///         fitting RMSD. Rejection is permanent and works with or without alignment.
///
/// Returns:
///     RmsdResult: rmsd evaluates every mapped selected evaluation pair, including
///         rejected fitting pairs. core_rmsd measures surviving fitting pairs;
///         initial_rmsd is before rejection. Counts, cycles, parameters, and
///         chain_alignments retain correspondence. Fitting needs at least three
///         non-collinear atom pairs; evaluation needs one. Invalid survivors fail.
///         Sequence mismatches can pair backbone atoms; side chains require matching
///         residue types, atom names, and elements. Omitted endpoints emit warnings.
///
/// Examples:
///     >>> result = arpeggia.rmsd("reference.cif", "query.cif", align_seqs=True,
///     ...     superpose_residues="A:1-100", rmsd_residues="A:101-120")
///     >>> print(result.rmsd, result.core_rmsd)
///
/// Selection examples and output schemas:
/// https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md
#[pyfunction]
#[pyo3(signature = (reference, query, model_num=0, superpose_residues="", rmsd_residues="", atoms="ca", *, align_seqs=false, chain_map=None, alignment_mode="global", gap_open=10.0, gap_extend=0.5, refine_cycles=0, refine_cutoff=2.0))]
#[allow(clippy::too_many_arguments)]
fn rmsd(
    py: Python<'_>,
    reference: String,
    query: String,
    model_num: usize,
    superpose_residues: &str,
    rmsd_residues: &str,
    atoms: &str,
    align_seqs: bool,
    chain_map: Option<std::collections::BTreeMap<String, String>>,
    alignment_mode: &str,
    gap_open: f64,
    gap_extend: f64,
    refine_cycles: usize,
    refine_cutoff: f64,
) -> PyResult<crate::RmsdResult> {
    let options = crate::RmsdOptions {
        model_num,
        superpose_residues: superpose_residues.into(),
        rmsd_residues: rmsd_residues.into(),
        atoms: atom_subset(atoms)?,
        align_seqs,
        chain_map: chain_map.unwrap_or_default(),
        alignment: crate::SeqAlignOptions {
            mode: value_enum(
                alignment_mode,
                "alignment_mode must be 'global', 'local', or 'semi-global'",
            )?,
            gap_open,
            gap_extend,
        },
        refine_cycles,
        refine_cutoff,
    };
    crate::validate_rmsd_selections(superpose_residues, rmsd_residues).map_err(python_error)?;
    let reference = load_for_python(py, &reference)?;
    let query = load_for_python(py, &query)?;
    let analysis = py
        .detach(move || crate::get_rmsd(reference, query, &options))
        .map_err(python_error)?;
    emit_python_warnings(py, analysis.warnings)?;
    Ok(analysis.value)
}

/// Calculate every unordered pairwise RMSD using exact atom correspondence.
///
/// Args:
///     input (str): Non-recursive PDB/mmCIF directory, or CSV, Parquet, or NDJSON
///         manifest. At least two structures are required. Directory IDs are
///         filename stems; manifest relative paths resolve against the manifest.
///     id_col (str): Manifest ID column, default "id". IDs must be unique.
///     path_col (str): Manifest path column, default "path". Resolved file paths
///         must be unique.
///     model_num (int): Model serial; 0 selects each structure's first model.
///     superpose_residues (str): Fit selection, applied to every structure using
///         chain IDs and author numbering. Empty selects all eligible residues.
///         Use "A:1-100,B" for inclusive ranges/whole chains; repeat chain IDs in
///         comma-separated clauses. Insertion codes and negative numbers are valid.
///     rmsd_residues (str): Evaluation selection with the same syntax, applied
///         to every structure. Empty independently selects all eligible residues;
///         it does not inherit superpose_residues.
///     atoms (str): "ca" (default), "backbone" (N/CA/C/O/OXT), "heavy", or "all".
///     num_threads (int): Worker limit; 0 uses available processors.
///     bypass_mem_check (bool): Skip heuristic memory checks (default False).
///         Pair storage is quadratic in structure count. Checks are estimates,
///         not a guarantee against running out of memory.
///
/// Returns:
///     polars.DataFrame: id_1 (String), id_2 (String), rmsd (Float64, Ångströms),
///         one row per unordered pair. Atom identities must agree across the
///         collection; sequence alignment and refinement are not supported here.
///
/// Examples:
///     >>> pairs = arpeggia.pairwise_rmsd("structures/", superpose_residues="A",
///     ...     rmsd_residues="B", num_threads=8)
#[pyfunction]
#[pyo3(signature = (input, id_col="id", path_col="path", model_num=0, superpose_residues="", rmsd_residues="", atoms="ca", num_threads=0, bypass_mem_check=false))]
#[allow(clippy::too_many_arguments)]
fn pairwise_rmsd(
    py: Python<'_>,
    input: String,
    id_col: &str,
    path_col: &str,
    model_num: usize,
    superpose_residues: &str,
    rmsd_residues: &str,
    atoms: &str,
    num_threads: usize,
    bypass_mem_check: bool,
) -> PyResult<PyDataFrame> {
    let options = crate::PairwiseRmsdOptions {
        model_num,
        superpose_residues: superpose_residues.to_string(),
        rmsd_residues: rmsd_residues.to_string(),
        atoms: atom_subset(atoms)?,
        num_threads,
        bypass_mem_check,
    };
    let analysis = py
        .detach(|| {
            crate::get_pairwise_rmsd(std::path::Path::new(&input), id_col, path_col, &options)
        })
        .map_err(python_error)?;
    emit_python_warnings(py, analysis.warnings)?;
    Ok(PyDataFrame(analysis.value))
}

/// Cluster at least three structures with k-medoids and return representatives.
///
/// Args:
///     input (str | None): Directory or manifest accepted by pairwise_rmsd().
///     pairwise_rmsd (polars.DataFrame | None): Complete unordered pair table with
///         id_1, id_2, rmsd columns. Supply exactly one of input or pairwise_rmsd.
///     id_col (str): Manifest ID column, default "id".
///     path_col (str): Manifest path column, default "path".
///     method (str): "k-medoids", currently the only supported method.
///     num_clusters (int | None): Fixed count from 1 through the structure count.
///         Takes precedence over max_clusters with a warning when both are given.
///     max_clusters (int | None): Upper bound for automatic selection, at least 2
///         and less than the structure count. Supply this or num_clusters.
///         An effectively identical ensemble forms one cluster automatically.
///     max_iterations (int): Iteration budget, default 100; nonconvergence fails.
///     model_num (int): Model serial; 0 selects each first model for input structures.
///     superpose_residues (str): Fit selection on every input structure. Empty
///         selects all eligible residues. Use "A:1-100,B" for a comma union of
///         inclusive author-number ranges/whole chains; insertion codes and
///         negative numbers are valid. Exact atom correspondence is required.
///     rmsd_residues (str): Evaluation selection on every input structure, using
///         the same syntax. Empty independently selects all eligible residues;
///         it does not inherit superpose_residues.
///     atoms (str): "ca" (default), "backbone" (N/CA/C/O/OXT), "heavy", or "all".
///     num_threads (int): Pairwise worker limit; 0 uses available processors.
///     bypass_mem_check (bool): Skip heuristic memory checks, default False.
///         Structure-selection options do not alter a supplied pairwise table.
///
/// Returns:
///     polars.DataFrame: id (String), cluster_id (UInt32, zero-based), medoid_id
///         (String), rmsd_to_medoid (Float64, Ångströms).
///
/// Examples:
///     >>> clusters = arpeggia.cluster_structs(input="structures/", num_clusters=3)
///
/// Schemas, CLI cache behavior, and memory use:
/// https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md
#[pyfunction]
#[pyo3(signature = (input=None, pairwise_rmsd=None, id_col="id", path_col="path", method="k-medoids", num_clusters=None, max_clusters=None, max_iterations=100, model_num=0, superpose_residues="", rmsd_residues="", atoms="ca", num_threads=0, bypass_mem_check=false))]
#[allow(clippy::too_many_arguments)]
fn cluster_structs(
    py: Python<'_>,
    input: Option<String>,
    pairwise_rmsd: Option<PyDataFrame>,
    id_col: &str,
    path_col: &str,
    method: &str,
    num_clusters: Option<usize>,
    max_clusters: Option<usize>,
    max_iterations: usize,
    model_num: usize,
    superpose_residues: &str,
    rmsd_residues: &str,
    atoms: &str,
    num_threads: usize,
    bypass_mem_check: bool,
) -> PyResult<PyDataFrame> {
    if input.is_some() == pairwise_rmsd.is_some() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "exactly one of input or pairwise_rmsd is required",
        ));
    }
    let cluster_options = crate::ClusterOptions {
        method: clustering_method(method)?,
        num_clusters,
        max_clusters,
        max_iterations,
    };
    cluster_options
        .validate_without_structure_count()
        .map_err(python_error)?;
    let matrix = match (input, pairwise_rmsd) {
        (Some(input), None) => {
            let atoms = atom_subset(atoms)?;
            crate::validate_rmsd_selections(superpose_residues, rmsd_residues)
                .map_err(python_error)?;
            let observations = py
                .detach(|| {
                    crate::read_structure_observations(
                        std::path::Path::new(&input),
                        id_col,
                        path_col,
                    )
                })
                .map_err(python_error)?;
            cluster_options
                .validate(observations.len())
                .map_err(python_error)?;
            let options = crate::PairwiseRmsdOptions {
                model_num,
                superpose_residues: superpose_residues.to_string(),
                rmsd_residues: rmsd_residues.to_string(),
                atoms,
                num_threads,
                bypass_mem_check,
            };
            let analysis = py
                .detach(|| crate::get_pairwise_rmsd_matrix(&observations, &options))
                .map_err(python_error)?;
            emit_python_warnings(py, analysis.warnings)?;
            analysis.value
        }
        (None, Some(dataframe)) => {
            let warnings = crate::rmsd::pairwise::check_packed_matrix_memory(
                dataframe.0.height(),
                bypass_mem_check,
            )
            .map_err(python_error)?;
            emit_python_warnings(py, warnings)?;
            py.detach(move || crate::PairwiseRmsdMatrix::from_dataframe(&dataframe.0, None, 3))
                .map_err(python_error)?
        }
        _ => unreachable!("exactly one clustering input was validated"),
    };
    let analysis = py
        .detach(|| crate::cluster_pairwise_rmsd(&matrix, &cluster_options))
        .map_err(python_error)?;
    emit_python_warnings(py, analysis.warnings)?;
    Ok(PyDataFrame(analysis.value))
}

/// Load a PDB or mmCIF file and calculate atomic and ring contacts.
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     groups (str, optional): Chain groups specification. Defaults to "/" (all-to-all).
///         Examples: "A,B/C,D" for chains A,B vs C,D; "A/" for chain A vs all others.
///     vdw_comp (float, optional): VdW distance tolerance in Ångströms. Defaults to 0.1.
///     dist_cutoff (float, optional): Distance cutoff for neighbor searches in Ångströms. Defaults to 6.5.
///     ignore_zero_occupancy (bool, optional): If True, ignore atoms with zero occupancy. Defaults to False.
///     protonation (str, optional): Histidine policy: "all-charged" (default),
///         "heuristic", or "explicit-only".
///     ph (float, optional): pH used by heuristic protonation. Defaults to 7.4.
///     num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.
///
/// Returns:
///     polars.DataFrame: A DataFrame containing all identified contacts with columns:
///         - model, interaction, distance
///         - from_chain, from_resn, from_resi, from_insertion, from_altloc, from_atomn, from_atomi
///         - to_chain, to_resn, to_resi, to_insertion, to_altloc, to_atomn, to_atomi
///         - sc_centroid_dist, sc_dihedral, sc_centroid_angle
///
/// Examples:
///     >>> import arpeggia
///     >>> contacts = arpeggia.contacts("structure.pdb", groups="/", vdw_comp=0.1)
///     >>> print(f"Found {len(contacts)} contacts")
#[pyfunction]
#[pyo3(signature = (input_file, groups="/", vdw_comp=0.1, dist_cutoff=6.5, ignore_zero_occupancy=false, protonation="all-charged", ph=7.4, num_threads=1))]
#[allow(clippy::too_many_arguments)]
fn contacts(
    py: Python<'_>,
    input_file: String,
    groups: &str,
    vdw_comp: f64,
    dist_cutoff: f64,
    ignore_zero_occupancy: bool,
    protonation: &str,
    ph: f64,
    num_threads: usize,
) -> PyResult<PyDataFrame> {
    // Load the PDB file
    let mut pdb = load_for_python(py, &input_file)?;
    let metadata = py
        .detach(|| crate::read_metadata(&input_file))
        .map_err(python_error)?;
    emit_python_warnings(py, metadata.warnings)?;
    let metadata = metadata.value;

    // Filter out atoms with zero occupancy if requested
    if ignore_zero_occupancy {
        pdb.remove_atoms_by(|atom| atom.occupancy() == 0.0);
    }

    // Get contacts
    let protonation = protonation_mode(protonation)?;
    let analysis = py
        .detach(|| {
            crate::run_with_threads(num_threads as isize, || {
                crate::analyze_contacts(
                    &pdb,
                    Some(&metadata),
                    groups,
                    vdw_comp,
                    dist_cutoff,
                    protonation,
                    ph,
                )
            })
        })
        .map_err(python_error)?;
    emit_python_warnings(py, analysis.warnings)?;

    // Convert to PyDataFrame for Python
    Ok(PyDataFrame(analysis.value))
}

/// Load a PDB or mmCIF file and calculate solvent accessible surface area (SASA).
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     level (str, optional): Aggregation level for SASA calculation. Options:
///         - "atom": Calculate SASA for each atom (default)
///         - "residue": Aggregate SASA by residue
///         - "chain": Aggregate SASA by chain
///     probe_radius (float, optional): Probe radius in Ångströms. Defaults to 1.4.
///         Smaller probes access narrower crevices; larger probes exclude them.
///         Total SASA changes depend on the structure.
///     n_points (int, optional): Number of points for surface calculation. Defaults to 100.
///     model_num (int, optional): Model serial to analyze (0 for first model). Defaults to 0.
///     chains (str, optional): Comma-separated chain IDs to include (e.g., "A,B,C").
///         If empty, includes all chains. Defaults to "".
///     num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.
///
/// Returns:
///     polars.DataFrame: A DataFrame with SASA areas in Å². Columns depend on the level:
///         - atom: atomi, sasa, polarity, chain, resn, resi, insertion, altloc, atomn
///         - residue: chain, resn, resi, insertion, sasa, polar_sasa,
///           hydrophobic_sasa, unclassified_sasa
///         - chain: chain, sasa, polar_sasa, hydrophobic_sasa, unclassified_sasa
///
/// Examples:
///     >>> import arpeggia
///     >>> # Atom-level SASA for all chains
///     >>> sasa_df = arpeggia.sasa("structure.pdb", level="atom")
///     >>> print(f"Calculated SASA for {len(sasa_df)} atoms")
///     >>>
///     >>> # Residue-level SASA for only chains A and B
///     >>> residue_sasa = arpeggia.sasa("structure.pdb", level="residue", chains="A,B")
///     >>> print(f"Calculated SASA for {len(residue_sasa)} residues")
///     >>>
///     >>> # Chain-level SASA
///     >>> chain_sasa = arpeggia.sasa("structure.pdb", level="chain")
///     >>> print(f"Calculated SASA for {len(chain_sasa)} chains")
#[pyfunction]
#[pyo3(signature = (input_file, level="atom", probe_radius=1.4, n_points=100, model_num=0, chains="", num_threads=1))]
#[allow(clippy::too_many_arguments)]
fn sasa(
    py: Python<'_>,
    input_file: String,
    level: &str,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
    num_threads: usize,
) -> PyResult<PyDataFrame> {
    // Load the PDB file
    let pdb = load_for_python(py, &input_file)?;

    let level = level.to_ascii_lowercase();
    if !matches!(level.as_str(), "atom" | "residue" | "chain") {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Invalid level '{level}'. Must be one of: 'atom', 'residue', 'chain'"
        )));
    }
    let analysis = py
        .detach(|| {
            crate::run_with_threads(num_threads as isize, || match level.as_str() {
                "atom" => crate::get_atom_sasa(&pdb, probe_radius, n_points, model_num, chains),
                "residue" => {
                    crate::get_residue_sasa(&pdb, probe_radius, n_points, model_num, chains)
                }
                "chain" => crate::get_chain_sasa(&pdb, probe_radius, n_points, model_num, chains),
                _ => unreachable!(),
            })
        })
        .map_err(python_error)?;
    emit_python_warnings(py, analysis.warnings)?;
    Ok(PyDataFrame(analysis.value))
}

/// Load a PDB or mmCIF file and calculate buried surface area at the interface between chain groups.
///
/// The two-sided buried surface area (dSASA) is calculated as:
/// dSASA = SASA_group1 + SASA_group2 - SASA_complex
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     groups (str): Chain groups specification for interface calculation.
///         Format: "A,B/C,D" where chains A,B form one side and C,D form the other.
///         Groups must be disjoint and non-empty; "A/" selects A vs all remaining chains.
///     probe_radius (float, optional): Probe radius in Ångströms. Defaults to 1.4.
///         Smaller probes access narrower crevices; larger probes exclude them.
///         Total SASA changes depend on the structure.
///     n_points (int, optional): Number of points for surface calculation. Defaults to 100.
///     model_num (int, optional): Model serial to analyze (0 for first model). Defaults to 0.
///     num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.
///
/// Returns:
///     float: The buried surface area at the interface in square Ångströms.
///
/// Examples:
///     >>> import arpeggia
///     >>> bsa = arpeggia.dsasa("structure.pdb", groups="A,B/C,D")
///     >>> print(f"Buried surface area: {bsa:.2f} Å²")
#[pyfunction]
#[pyo3(signature = (input_file, groups, probe_radius=1.4, n_points=100, model_num=0, num_threads=1))]
fn dsasa(
    py: Python<'_>,
    input_file: String,
    groups: &str,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    num_threads: usize,
) -> PyResult<f32> {
    Ok(dsasa_components(
        py,
        input_file,
        groups,
        probe_radius,
        n_points,
        model_num,
        num_threads,
    )?
    .0)
}

/// Return (total, polar, hydrophobic, unclassified) two-sided dSASA in Å².
///
/// Args:
///     input_file (str): PDB or mmCIF path.
///     groups (str): Disjoint non-empty groups, e.g. "A,B/C" or "A/" (A vs rest).
///     probe_radius (float): Solvent probe radius in Å, default 1.4.
///         Smaller probes access narrower crevices; larger probes exclude them.
///         Total SASA changes depend on the structure.
///     n_points (int): Positive surface sample count per sphere, default 100.
///     model_num (int): Model serial, or 0 for the first model.
///     num_threads (int): Worker limit, default 1; 0 uses available processors.
///
/// Returns:
///     tuple[float, float, float, float]: (total, polar, hydrophobic, unclassified)
///         areas in Å². Total is SASA(group1)+SASA(group2)-SASA(complex), and the
///         partitions sum to that total.
///
/// Examples:
///     >>> total, polar, hydrophobic, unknown = arpeggia.dsasa_components(
///     ...     "complex.cif", groups="H,L/A")
#[pyfunction]
#[pyo3(signature = (input_file, groups, probe_radius=1.4, n_points=100, model_num=0, num_threads=1))]
fn dsasa_components(
    py: Python<'_>,
    input_file: String,
    groups: &str,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    num_threads: usize,
) -> PyResult<(f32, f32, f32, f32)> {
    let pdb = load_for_python(py, &input_file)?;
    let result = py
        .detach(|| {
            crate::run_with_threads(num_threads as isize, || {
                crate::get_dsasa_components(&pdb, groups, probe_radius, n_points, model_num)
            })
        })
        .map_err(python_error)?;
    emit_python_warnings(py, result.warnings)?;
    let result = result.value;
    Ok((
        result.dsasa,
        result.polar_dsasa,
        result.hydrophobic_dsasa,
        result.unclassified_dsasa,
    ))
}

/// Extract coordinate-observed protein sequences for all chains in a selected model.
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file.
///     model_num (int): Model serial; 0 selects the first model (default).
///
/// Returns:
///     list[tuple[str, str]]: (chain ID, observed protein sequence) tuples.
///         Declared residues without coordinates are absent; use seqres() for those.
///
/// Examples:
///     >>> import arpeggia
///     >>> sequences = arpeggia.seq("structure.pdb")
///     >>> for chain_id, seq in sequences:
///     ...     print(f"Chain {chain_id}: {seq}")
#[pyfunction]
#[pyo3(signature = (input_file, model_num=0))]
fn seq(
    py: Python<'_>,
    input_file: String,
    model_num: usize,
) -> PyResult<std::vec::Vec<(String, String)>> {
    // Load the PDB file
    let pdb = load_for_python(py, &input_file)?;

    // Get sequences
    let seqs = py
        .detach(|| crate::get_sequences(&pdb, model_num))
        .map_err(python_error)?;

    Ok(seqs)
}

/// Return declared PDB SEQRES or mmCIF entity-polymer sequences by chain.
///
/// No coordinate model is selected. Use seq() for observed protein sequences.
///
/// Args:
///     input_file (str): Path to a PDB or mmCIF file.
///
/// Returns:
///     list[tuple[str, str]]: (chain ID, declared sequence) tuples, including
///         declared residues without coordinates.
///
/// Examples:
///     >>> declared = arpeggia.seqres("structure.cif")
#[pyfunction]
fn seqres(py: Python<'_>, input_file: String) -> PyResult<Vec<(String, String)>> {
    let analysis = py
        .detach(|| crate::get_seqres(input_file))
        .map_err(python_error)?;
    emit_python_warnings(py, analysis.warnings)?;
    Ok(analysis.value)
}

/// Load a PDB or mmCIF file and calculate relative SASA (RSA) for each residue.
///
/// RSA is calculated as the ratio of observed SASA to the maximum possible SASA
/// for each amino acid type, based on Tien et al. (2013) theoretical values.
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     probe_radius (float, optional): Probe radius in Ångströms. Defaults to 1.4.
///         Smaller probes access narrower crevices; larger probes exclude them.
///         Total SASA changes depend on the structure.
///     n_points (int, optional): Number of points for surface calculation. Defaults to 100.
///     model_num (int, optional): Model serial to analyze (0 for first model). Defaults to 0.
///     chains (str, optional): Comma-separated chain IDs to include (e.g., "A,B,C").
///         If empty, includes all chains. Defaults to "".
///     num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.
///
/// Returns:
///     polars.DataFrame: A DataFrame with areas in Å² and dimensionless relative_sasa ratios
///         (not percentages), with columns:
///         - chain, resn, resi, insertion, sasa, polar_sasa,
///           hydrophobic_sasa, unclassified_sasa, relative_sasa
///
/// Examples:
///     >>> import arpeggia
///     >>> # RSA for all chains
///     >>> rsa = arpeggia.relative_sasa("structure.pdb", probe_radius=1.4)
///     >>> print(f"Calculated RSA for {len(rsa)} residues")
///     >>>
///     >>> # RSA for only chain A
///     >>> rsa_a = arpeggia.relative_sasa("structure.pdb", chains="A")
#[pyfunction]
#[pyo3(signature = (input_file, probe_radius=1.4, n_points=100, model_num=0, chains="", num_threads=1))]
fn relative_sasa(
    py: Python<'_>,
    input_file: String,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
    num_threads: usize,
) -> PyResult<PyDataFrame> {
    // Load the PDB file
    let pdb = load_for_python(py, &input_file)?;

    // Get relative SASA
    let analysis = py
        .detach(|| {
            crate::run_with_threads(num_threads as isize, || {
                crate::get_relative_sasa(&pdb, probe_radius, n_points, model_num, chains)
            })
        })
        .map_err(python_error)?;
    emit_python_warnings(py, analysis.warnings)?;

    // Convert to PyDataFrame for Python
    Ok(PyDataFrame(analysis.value))
}

/// Load a PDB or mmCIF file and calculate Spatial Aggregation Propensity (SAP) scores.
///
/// The SAP score quantifies the aggregation propensity by combining the solvent-accessible
/// hydrophobic surface area of neighboring residues. It was developed by Chennamsetty et al.
/// and is described in "Developability Index: A Rapid In Silico Tool for the Screening of
/// Antibody Aggregation Propensity" (J Pharm Sci, 2012).
///
/// The formula is:
/// SAP(i) = Σ{j ∈ neighbors(i, R)} [ Hydrophobicity(j) × (SASA(j) / SASA_max(j)) ]
///
/// Where:
/// - Neighbors are atoms/residues within radius R of atom/residue i
/// - Hydrophobicity uses Rosetta's Black & Mould-derived constants
/// - SASA is the side-chain solvent accessible surface area
/// - SASA_max is the maximum SASA for that residue type
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     level (str, optional): Aggregation level for SAP calculation. Options:
///         - "atom": Calculate SAP for each atom
///         - "residue": Aggregate SAP by residue (default)
///     probe_radius (float, optional): Probe radius in Ångströms for SASA calculation. Defaults to 1.1.
///         Smaller probes access narrower crevices; larger probes exclude them.
///         Total SASA changes depend on the structure.
///     n_points (int, optional): Number of points for SASA surface calculation. Defaults to 100.
///     model_num (int, optional): Model serial to analyze (0 for first model). Defaults to 0.
///     sap_radius (float, optional): Radius in Ångströms for neighbor search. Defaults to 5.0.
///     chains (str, optional): Comma-separated chain IDs to include (e.g., "H,L").
///         If empty, includes all chains. Defaults to "".
///     num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.
///
/// Returns:
///     polars.DataFrame: A DataFrame with SAP scores. Columns depend on the level:
///         - atom: chain, resn, resi, insertion, atomn, atomi, sasa, sap_score
///         - residue: chain, resn, resi, insertion, sc_sasa, sap_score,
///           max_sc_asa, relative_sc_sasa
///
/// Examples:
///     >>> import arpeggia
///     >>> # Residue-level SAP scores for all chains
///     >>> residue_sap = arpeggia.sap_score("structure.pdb")
///     >>> print(f"Calculated SAP for {len(residue_sap)} residues")
///     >>>
///     >>> # SAP scores for only antibody heavy and light chains
///     >>> sap_hl = arpeggia.sap_score("antibody.pdb", chains="H,L")
///     >>> print(f"Calculated SAP for {len(sap_hl)} residues")
#[pyfunction]
#[pyo3(signature = (input_file, level="residue", probe_radius=1.1, n_points=100, model_num=0, sap_radius=5.0, chains="", num_threads=1))]
#[allow(clippy::too_many_arguments)]
fn sap_score(
    py: Python<'_>,
    input_file: String,
    level: &str,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    sap_radius: f32,
    chains: &str,
    num_threads: usize,
) -> PyResult<PyDataFrame> {
    // Load the PDB file
    let pdb = load_for_python(py, &input_file)?;

    let level = level.to_ascii_lowercase();
    if !matches!(level.as_str(), "atom" | "residue") {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Invalid level '{level}'. Must be one of: 'atom', 'residue'"
        )));
    }
    let analysis = py
        .detach(|| {
            crate::run_with_threads(num_threads as isize, || match level.as_str() {
                "atom" => crate::get_per_atom_sap_score(
                    &pdb,
                    probe_radius,
                    n_points,
                    model_num,
                    sap_radius,
                    chains,
                ),
                "residue" => crate::get_per_residue_sap_score(
                    &pdb,
                    probe_radius,
                    n_points,
                    model_num,
                    sap_radius,
                    chains,
                ),
                _ => unreachable!(),
            })
        })
        .map_err(python_error)?;
    emit_python_warnings(py, analysis.warnings)?;
    Ok(PyDataFrame(analysis.value))
}

/// Calculate Shape Complementarity (SC) between two chain groups.
///
/// Shape complementarity measures how well two molecular surfaces fit together,
/// following Lawrence & Colman (1993) "Shape Complementarity at Protein/Protein Interfaces".
/// Higher SC values (closer to 1.0) indicate better geometric fit between surfaces.
/// Typical protein-protein interfaces have SC values between 0.5 and 0.7.
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     groups (str): Chain groups specification, e.g., "H,L/A" for chains H,L vs chain A.
///         Groups must be disjoint and non-empty, separated by "/"; "H,L/"
///         compares H,L against all remaining chains.
///     model_num (int, optional): Model serial to analyze (0 for first model). Defaults to 0.
///     num_threads (int, optional): Number of threads for parallel calculations (0 for auto). Defaults to 0.
///
/// Returns:
///     float: The shape complementarity score (approximately -1 to 1).
///
/// Examples:
///     >>> import arpeggia
///     >>> sc = arpeggia.sc("antibody_antigen.pdb", groups="H,L/A")
///     >>> print(f"SC Score: {sc:.3f}")
#[pyfunction]
#[pyo3(signature = (input_file, groups, model_num=0, num_threads=0))]
fn sc(
    py: Python<'_>,
    input_file: String,
    groups: &str,
    model_num: usize,
    num_threads: usize,
) -> PyResult<f64> {
    // Load the PDB file; SC applies shared structure preparation internally.
    let pdb = load_for_python(py, &input_file)?;

    // Calculate SC
    let analysis = py
        .detach(|| {
            crate::run_with_threads(num_threads as isize, || {
                crate::get_sc_details(&pdb, groups, model_num)
            })
        })
        .map_err(python_error)?;
    emit_python_warnings(py, analysis.warnings)?;
    Ok(analysis.value.sc)
}

/// Python module for protein structure analysis.
///
/// This module provides functions for analyzing protein structures from PDB and mmCIF files,
/// including contact detection, SASA calculation, SAP score calculation, and sequence extraction.
#[pymodule]
fn arpeggia(m: &Bound<'_, PyModule>) -> PyResult<()> {
    ab_numbering::register(m)?;
    m.add_class::<crate::SeqAlignment>()?;
    m.add_class::<crate::RmsdResult>()?;
    m.add_class::<crate::ChainAlignment>()?;
    m.add_class::<crate::ResiduePair>()?;
    m.add_function(wrap_pyfunction!(align_seqs, m)?)?;
    m.add_function(wrap_pyfunction!(rmsd, m)?)?;
    m.add_function(wrap_pyfunction!(pairwise_rmsd, m)?)?;
    m.add_function(wrap_pyfunction!(cluster_structs, m)?)?;
    m.add_function(wrap_pyfunction!(contacts, m)?)?;
    m.add_function(wrap_pyfunction!(sasa, m)?)?;
    m.add_function(wrap_pyfunction!(dsasa, m)?)?;
    m.add_function(wrap_pyfunction!(dsasa_components, m)?)?;
    m.add_function(wrap_pyfunction!(relative_sasa, m)?)?;
    m.add_function(wrap_pyfunction!(sap_score, m)?)?;
    m.add_function(wrap_pyfunction!(sc, m)?)?;
    m.add_function(wrap_pyfunction!(seq, m)?)?;
    m.add_function(wrap_pyfunction!(seqres, m)?)?;
    Ok(())
}
