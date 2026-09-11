"""Type stubs for the arpeggia Rust module."""

from collections.abc import Sequence
from typing import Literal

import polars as pl

from ._contract import (
    AlignmentMode,
    AtomSubset,
    CdrDefinition,
    ClusteringMethod,
    DsasaComponents,
    GermlineSpecies,
    NumberingScheme,
    ProtonationMode,
    SapLevel,
    SasaLevel,
    SequenceList,
)

def align_seqs(
    reference: str,
    query: str,
    mode: AlignmentMode = ...,
    gap_open: float = ...,
    gap_extend: float = ...,
    *,
    reference_name: str = "Reference",
    query_name: str = "Query",
) -> SeqAlignment:
    """Align two unaligned amino-acid strings and return a read-only SeqAlignment.

    Printing the result shows colored, wrapped sequences and position rulers when
    the terminal permits. Use result.format(width=60, color="never", rulers=False)
    for plain output without rulers. Data fields contain no color escapes.

    Args:
        reference (str): Reference input sequence. Lowercase is
            normalized; standard amino acids and B/Z/X/U/O are accepted. Inputs
            must be non-empty and contain no whitespace, gaps, or stop symbols.
        query (str): Query input sequence, with the same alphabet restrictions.
        mode (str): "global" (default) consumes both inputs; "local" selects
            best-scoring subsequences; "semi-global" consumes the entire query
            sequence with free reference terminal overhangs.
        gap_open (float): Positive opening cost, default 10, with at most two
            decimals and no smaller than gap_extend. A gap of length L costs
            gap_open + (L - 1) * gap_extend.
        gap_extend (float): Positive cost per additional gap residue, default
            0.5, with at most two decimals and no larger than gap_open.
        reference_name (str): Keyword-only reference display label, default
            "Reference". Stored on the result and in JSON; does not affect scoring.
        query_name (str): Keyword-only query display label, default "Query".

    Returns:
        SeqAlignment: Plain gapped aligned_reference/aligned_query strings,
            reference-to-query operations (space for match, + insertion, - deletion,
            : positive-score substitution, x other substitution), zero-based half-open
            input spans, score, matches, mismatches, gap_residues, and gap_runs.
            identity_alignment/identity_shorter and coverage_alignment/coverage_shorter
            are ratios using alignment-column and shorter-full-input denominators.
            edit_distance is full-input Levenshtein distance, independent of score.
            Empty local results have score zero and None alignment-length ratios.
            U/O score as C/K with a warning but retain distinct identity.
            Positive-score substitutions are blue in displays and remain mismatches.

    Examples:
        >>> alignment = arpeggia.align_seqs("GGACDEFGHIKGG", "ACDEFGHIK", mode="semi-global")
        >>> alignment.reference_span
        (2, 11)
    """

def rmsd(
    reference: str,
    query: str,
    model_num: int = ...,
    superpose_residues: str = ...,
    rmsd_residues: str = ...,
    atoms: AtomSubset = ...,
    *,
    align_seqs: bool = ...,
    chain_map: dict[str, str] | None = ...,
    alignment_mode: AlignmentMode = ...,
    gap_open: float = ...,
    gap_extend: float = ...,
    refine_cycles: int = ...,
    refine_cutoff: float = ...,
) -> RmsdResult:
    """Superpose two PDB/mmCIF structures; return a read-only RmsdResult in Ångströms.

    Args:
        reference (str): Path to the reference PDB/mmCIF structure.
        query (str): Path to the query PDB/mmCIF structure to superpose.
        model_num (int): Model serial; 0 independently selects each first model.
        superpose_residues (str): Fit selection on the reference structure, using
            reference chain IDs and author residue numbers. With align_seqs=True,
            corresponding query residues come from the sequence alignment. With
            align_seqs=False, apply the same selection to the query and require
            exact atom identities. Empty selects all eligible residues.
            Use "A" for a whole chain or "A:1-100,A:110-120,B"
            for a union of inclusive author-number ranges and chains. Repeat the
            chain in each clause. Negative numbers ("A:-5--1") and insertion
            codes ("B:10A-20") are accepted; "A:10" includes all insertion variants.
        rmsd_residues (str): Evaluation selection on the reference structure,
            using the same syntax and query-correspondence rules as the fit
            selection. Empty independently means all; it does not inherit
            superpose_residues. Evaluation never determines or refits the transform.
        atoms (str): "ca" (default), "backbone" (N/CA/C/O/OXT), "heavy", or "all".
        align_seqs (bool): Align complete observed chains before atom selection.
            When True, both residue selectors use reference author numbering;
            otherwise selections apply to both structures with exact correspondence.
        chain_map (dict[str, str] | None): Explicit reference-to-query chain IDs,
            e.g. {"A": "H"}. Requires align_seqs=True, covers exactly selected
            reference chains, and assigns unique query chains. None infers a
            unique maximum-score assignment; ambiguous homomers require a map.
        alignment_mode (str): "global" (default), "local", or "semi-global" for
            final chain alignment. Semi-global consumes the complete query chain.
            Chain inference independently scores shorter against longer semi-globally.
        gap_open (float): Positive alignment opening cost, default 10; at most
            two decimals and no smaller than gap_extend. Gap cost is open+(L-1)*extend.
        gap_extend (float): Positive cost per additional gap residue, default 0.5;
            at most two decimals and no larger than gap_open.
        refine_cycles (int): Maximum rejection/refit rounds after the initial fit;
            0 (default) performs no rejection. Stops early when unchanged or exact.
        refine_cutoff (float): Positive rejection multiplier (default 2) of current
            fitting RMSD. Rejection is permanent and works with or without alignment.

    Returns:
        RmsdResult: rmsd evaluates every mapped selected evaluation pair, including
            rejected fitting pairs. core_rmsd measures surviving fitting pairs;
            initial_rmsd is before rejection. Counts, cycles, parameters, and
            chain_alignments retain correspondence. Fitting needs at least three
            non-collinear atom pairs; evaluation needs one. Invalid survivors fail.
            Sequence mismatches can pair backbone atoms; side chains require matching
            residue types, atom names, and elements. Omitted endpoints emit warnings.

    Examples:
        >>> result = arpeggia.rmsd("reference.cif", "query.cif", align_seqs=True,
        ...     superpose_residues="A:1-100", rmsd_residues="A:101-120")
        >>> print(result.rmsd, result.core_rmsd)

    Selection examples and output schemas:
    https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md
    """

def pairwise_rmsd(
    input: str,
    id_col: str = ...,
    path_col: str = ...,
    model_num: int = ...,
    superpose_residues: str = ...,
    rmsd_residues: str = ...,
    atoms: AtomSubset = ...,
    num_threads: int = ...,
    bypass_mem_check: bool = ...,
) -> pl.DataFrame:
    """Calculate every unordered pairwise RMSD using exact atom correspondence.

    Args:
        input (str): Non-recursive PDB/mmCIF directory, or CSV, Parquet, or NDJSON
            manifest. At least two structures are required. Directory IDs are
            filename stems; manifest relative paths resolve against the manifest.
        id_col (str): Manifest ID column, default "id". IDs must be unique.
        path_col (str): Manifest path column, default "path". Resolved file paths
            must be unique.
        model_num (int): Model serial; 0 selects each structure's first model.
        superpose_residues (str): Fit selection, applied to every structure using
            chain IDs and author numbering. Empty selects all eligible residues.
            Use "A:1-100,B" for inclusive ranges/whole chains; repeat chain IDs in
            comma-separated clauses. Insertion codes and negative numbers are valid.
        rmsd_residues (str): Evaluation selection with the same syntax, applied
            to every structure. Empty independently selects all eligible residues;
            it does not inherit superpose_residues.
        atoms (str): "ca" (default), "backbone" (N/CA/C/O/OXT), "heavy", or "all".
        num_threads (int): Worker limit; 0 uses available processors.
        bypass_mem_check (bool): Skip heuristic memory checks (default False).
            Pair storage is quadratic in structure count. Checks are estimates,
            not a guarantee against running out of memory.

    Returns:
        polars.DataFrame: id_1 (String), id_2 (String), rmsd (Float64, Ångströms),
            one row per unordered pair. Atom identities must agree across the
            collection; sequence alignment and refinement are not supported here.

    Examples:
        >>> pairs = arpeggia.pairwise_rmsd("structures/", superpose_residues="A",
        ...     rmsd_residues="B", num_threads=8)
    """

def cluster_structs(
    input: str | None = ...,
    pairwise_rmsd: pl.DataFrame | None = ...,
    id_col: str = ...,
    path_col: str = ...,
    method: ClusteringMethod = ...,
    num_clusters: int | None = ...,
    max_clusters: int | None = ...,
    max_iterations: int = ...,
    model_num: int = ...,
    superpose_residues: str = ...,
    rmsd_residues: str = ...,
    atoms: AtomSubset = ...,
    num_threads: int = ...,
    bypass_mem_check: bool = ...,
) -> pl.DataFrame:
    """Cluster at least three structures with k-medoids and return representatives.

    Args:
        input (str | None): Directory or manifest accepted by pairwise_rmsd().
        pairwise_rmsd (polars.DataFrame | None): Complete unordered pair table with
            id_1, id_2, rmsd columns. Supply exactly one of input or pairwise_rmsd.
        id_col (str): Manifest ID column, default "id".
        path_col (str): Manifest path column, default "path".
        method (str): "k-medoids", currently the only supported method.
        num_clusters (int | None): Fixed count from 1 through the structure count.
            Takes precedence over max_clusters with a warning when both are given.
        max_clusters (int | None): Upper bound for automatic selection, at least 2
            and less than the structure count. Supply this or num_clusters.
            An effectively identical ensemble forms one cluster automatically.
        max_iterations (int): Iteration budget, default 100; nonconvergence fails.
        model_num (int): Model serial; 0 selects each first model for input structures.
        superpose_residues (str): Fit selection on every input structure. Empty
            selects all eligible residues. Use "A:1-100,B" for a comma union of
            inclusive author-number ranges/whole chains; insertion codes and
            negative numbers are valid. Exact atom correspondence is required.
        rmsd_residues (str): Evaluation selection on every input structure, using
            the same syntax. Empty independently selects all eligible residues;
            it does not inherit superpose_residues.
        atoms (str): "ca" (default), "backbone" (N/CA/C/O/OXT), "heavy", or "all".
        num_threads (int): Pairwise worker limit; 0 uses available processors.
        bypass_mem_check (bool): Skip heuristic memory checks, default False.
            Structure-selection options do not alter a supplied pairwise table.

    Returns:
        polars.DataFrame: id (String), cluster_id (UInt32, zero-based), medoid_id
            (String), rmsd_to_medoid (Float64, Ångströms).

    Examples:
        >>> clusters = arpeggia.cluster_structs(input="structures/", num_clusters=3)

    Schemas, CLI cache behavior, and memory use:
    https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md
    """

def contacts(
    input_file: str,
    groups: str = ...,
    vdw_comp: float = ...,
    dist_cutoff: float = ...,
    ignore_zero_occupancy: bool = ...,
    protonation: ProtonationMode = ...,
    ph: float = ...,
    num_threads: int = ...,
) -> pl.DataFrame:
    """Load a PDB or mmCIF file and calculate atomic and ring contacts.

    Args:
        input_file (str): Path to the PDB or mmCIF file
        groups (str, optional): Chain groups specification. Defaults to "/" (all-to-all).
            Examples: "A,B/C,D" for chains A,B vs C,D; "A/" for chain A vs all others.
        vdw_comp (float, optional): VdW distance tolerance in Ångströms. Defaults to 0.1.
        dist_cutoff (float, optional): Distance cutoff for neighbor searches in Ångströms. Defaults to 6.5.
        ignore_zero_occupancy (bool, optional): If True, ignore atoms with zero occupancy. Defaults to False.
        protonation (str, optional): Histidine policy: "all-charged" (default),
            "heuristic", or "explicit-only".
        ph (float, optional): pH used by heuristic protonation. Defaults to 7.4.
        num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.

    Returns:
        polars.DataFrame: A DataFrame containing all identified contacts with columns:
            - model, interaction, distance
            - from_chain, from_resn, from_resi, from_insertion, from_altloc, from_atomn, from_atomi
            - to_chain, to_resn, to_resi, to_insertion, to_altloc, to_atomn, to_atomi
            - sc_centroid_dist, sc_dihedral, sc_centroid_angle

    Examples:
        >>> import arpeggia
        >>> contacts = arpeggia.contacts("structure.pdb", groups="/", vdw_comp=0.1)
        >>> print(f"Found {len(contacts)} contacts")
    """

def sasa(
    input_file: str,
    level: SasaLevel = ...,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    chains: str = ...,
    num_threads: int = ...,
) -> pl.DataFrame:
    """Load a PDB or mmCIF file and calculate solvent accessible surface area (SASA).

    Args:
        input_file (str): Path to the PDB or mmCIF file
        level (str, optional): Aggregation level for SASA calculation. Options:
            - "atom": Calculate SASA for each atom (default)
            - "residue": Aggregate SASA by residue
            - "chain": Aggregate SASA by chain
        probe_radius (float, optional): Probe radius in Ångströms. Defaults to 1.4.
            Smaller probes access narrower crevices; larger probes exclude them.
            Total SASA changes depend on the structure.
        n_points (int, optional): Number of points for surface calculation. Defaults to 100.
        model_num (int, optional): Model serial to analyze (0 for first model). Defaults to 0.
        chains (str, optional): Comma-separated chain IDs to include (e.g., "A,B,C").
            If empty, includes all chains. Defaults to "".
        num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.

    Returns:
        polars.DataFrame: A DataFrame with SASA areas in Å². Columns depend on the level:
            - atom: atomi, sasa, polarity, chain, resn, resi, insertion, altloc, atomn
            - residue: chain, resn, resi, insertion, sasa, polar_sasa,
              hydrophobic_sasa, unclassified_sasa
            - chain: chain, sasa, polar_sasa, hydrophobic_sasa, unclassified_sasa

    Examples:
        >>> import arpeggia
        >>> # Atom-level SASA for all chains
        >>> sasa_df = arpeggia.sasa("structure.pdb", level="atom")
        >>> print(f"Calculated SASA for {len(sasa_df)} atoms")
        >>>
        >>> # Residue-level SASA for only chains A and B
        >>> residue_sasa = arpeggia.sasa("structure.pdb", level="residue", chains="A,B")
        >>> print(f"Calculated SASA for {len(residue_sasa)} residues")
        >>>
        >>> # Chain-level SASA
        >>> chain_sasa = arpeggia.sasa("structure.pdb", level="chain")
        >>> print(f"Calculated SASA for {len(chain_sasa)} chains")
    """

def dsasa(
    input_file: str,
    groups: str,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    num_threads: int = ...,
) -> float:
    """Load a PDB or mmCIF file and calculate buried surface area at the interface between chain groups.

    The two-sided buried surface area (dSASA) is calculated as:
    dSASA = SASA_group1 + SASA_group2 - SASA_complex

    Args:
        input_file (str): Path to the PDB or mmCIF file
        groups (str): Chain groups specification for interface calculation.
            Format: "A,B/C,D" where chains A,B form one side and C,D form the other.
            Groups must be disjoint and non-empty; "A/" selects A vs all remaining chains.
        probe_radius (float, optional): Probe radius in Ångströms. Defaults to 1.4.
            Smaller probes access narrower crevices; larger probes exclude them.
            Total SASA changes depend on the structure.
        n_points (int, optional): Number of points for surface calculation. Defaults to 100.
        model_num (int, optional): Model serial to analyze (0 for first model). Defaults to 0.
        num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.

    Returns:
        float: The buried surface area at the interface in square Ångströms.

    Examples:
        >>> import arpeggia
        >>> bsa = arpeggia.dsasa("structure.pdb", groups="A,B/C,D")
        >>> print(f"Buried surface area: {bsa:.2f} Å²")
    """

def dsasa_components(
    input_file: str,
    groups: str,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    num_threads: int = ...,
) -> DsasaComponents:
    """Return (total, polar, hydrophobic, unclassified) two-sided dSASA in Å².

    Args:
        input_file (str): PDB or mmCIF path.
        groups (str): Disjoint non-empty groups, e.g. "A,B/C" or "A/" (A vs rest).
        probe_radius (float): Solvent probe radius in Å, default 1.4.
            Smaller probes access narrower crevices; larger probes exclude them.
            Total SASA changes depend on the structure.
        n_points (int): Positive surface sample count per sphere, default 100.
        model_num (int): Model serial, or 0 for the first model.
        num_threads (int): Worker limit, default 1; 0 uses available processors.

    Returns:
        tuple[float, float, float, float]: (total, polar, hydrophobic, unclassified)
            areas in Å². Total is SASA(group1)+SASA(group2)-SASA(complex), and the
            partitions sum to that total.

    Examples:
        >>> total, polar, hydrophobic, unknown = arpeggia.dsasa_components(
        ...     "complex.cif", groups="H,L/A")
    """

def relative_sasa(
    input_file: str,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    chains: str = ...,
    num_threads: int = ...,
) -> pl.DataFrame:
    """Load a PDB or mmCIF file and calculate relative SASA (RSA) for each residue.

    RSA is calculated as the ratio of observed SASA to the maximum possible SASA
    for each amino acid type, based on Tien et al. (2013) theoretical values.

    Args:
        input_file (str): Path to the PDB or mmCIF file
        probe_radius (float, optional): Probe radius in Ångströms. Defaults to 1.4.
            Smaller probes access narrower crevices; larger probes exclude them.
            Total SASA changes depend on the structure.
        n_points (int, optional): Number of points for surface calculation. Defaults to 100.
        model_num (int, optional): Model serial to analyze (0 for first model). Defaults to 0.
        chains (str, optional): Comma-separated chain IDs to include (e.g., "A,B,C").
            If empty, includes all chains. Defaults to "".
        num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.

    Returns:
        polars.DataFrame: A DataFrame with areas in Å² and dimensionless relative_sasa ratios
            (not percentages), with columns:
            - chain, resn, resi, insertion, sasa, polar_sasa,
              hydrophobic_sasa, unclassified_sasa, relative_sasa

    Examples:
        >>> import arpeggia
        >>> # RSA for all chains
        >>> rsa = arpeggia.relative_sasa("structure.pdb", probe_radius=1.4)
        >>> print(f"Calculated RSA for {len(rsa)} residues")
        >>>
        >>> # RSA for only chain A
        >>> rsa_a = arpeggia.relative_sasa("structure.pdb", chains="A")
    """

def sap_score(
    input_file: str,
    level: SapLevel = ...,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    sap_radius: float = ...,
    chains: str = ...,
    num_threads: int = ...,
) -> pl.DataFrame:
    """Load a PDB or mmCIF file and calculate Spatial Aggregation Propensity (SAP) scores.

    The SAP score quantifies the aggregation propensity by combining the solvent-accessible
    hydrophobic surface area of neighboring residues. It was developed by Chennamsetty et al.
    and is described in "Developability Index: A Rapid In Silico Tool for the Screening of
    Antibody Aggregation Propensity" (J Pharm Sci, 2012).

    The formula is:
    SAP(i) = Σ{j ∈ neighbors(i, R)} [ Hydrophobicity(j) × (SASA(j) / SASA_max(j)) ]

    Where:
    - Neighbors are atoms/residues within radius R of atom/residue i
    - Hydrophobicity uses Rosetta's Black & Mould-derived constants
    - SASA is the side-chain solvent accessible surface area
    - SASA_max is the maximum SASA for that residue type

    Args:
        input_file (str): Path to the PDB or mmCIF file
        level (str, optional): Aggregation level for SAP calculation. Options:
            - "atom": Calculate SAP for each atom
            - "residue": Aggregate SAP by residue (default)
        probe_radius (float, optional): Probe radius in Ångströms for SASA calculation. Defaults to 1.1.
            Smaller probes access narrower crevices; larger probes exclude them.
            Total SASA changes depend on the structure.
        n_points (int, optional): Number of points for SASA surface calculation. Defaults to 100.
        model_num (int, optional): Model serial to analyze (0 for first model). Defaults to 0.
        sap_radius (float, optional): Radius in Ångströms for neighbor search. Defaults to 5.0.
        chains (str, optional): Comma-separated chain IDs to include (e.g., "H,L").
            If empty, includes all chains. Defaults to "".
        num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.

    Returns:
        polars.DataFrame: A DataFrame with SAP scores. Columns depend on the level:
            - atom: chain, resn, resi, insertion, atomn, atomi, sasa, sap_score
            - residue: chain, resn, resi, insertion, sc_sasa, sap_score,
              max_sc_asa, relative_sc_sasa

    Examples:
        >>> import arpeggia
        >>> # Residue-level SAP scores for all chains
        >>> residue_sap = arpeggia.sap_score("structure.pdb")
        >>> print(f"Calculated SAP for {len(residue_sap)} residues")
        >>>
        >>> # SAP scores for only antibody heavy and light chains
        >>> sap_hl = arpeggia.sap_score("antibody.pdb", chains="H,L")
        >>> print(f"Calculated SAP for {len(sap_hl)} residues")
    """

def seq(input_file: str, model_num: int = ...) -> SequenceList:
    """Extract coordinate-observed protein sequences for all chains in a selected model.

    Args:
        input_file (str): Path to the PDB or mmCIF file.
        model_num (int): Model serial; 0 selects the first model (default).

    Returns:
        list[tuple[str, str]]: (chain ID, observed protein sequence) tuples.
            Declared residues without coordinates are absent; use seqres() for those.

    Examples:
        >>> import arpeggia
        >>> sequences = arpeggia.seq("structure.pdb")
        >>> for chain_id, seq in sequences:
        ...     print(f"Chain {chain_id}: {seq}")
    """

def seqres(input_file: str) -> SequenceList:
    """Return declared PDB SEQRES or mmCIF entity-polymer sequences by chain.

    No coordinate model is selected. Use seq() for observed protein sequences.

    Args:
        input_file (str): Path to a PDB or mmCIF file.

    Returns:
        list[tuple[str, str]]: (chain ID, declared sequence) tuples, including
            declared residues without coordinates.

    Examples:
        >>> declared = arpeggia.seqres("structure.cif")
    """

def sc(
    input_file: str,
    groups: str,
    model_num: int = ...,
    num_threads: int = ...,
) -> float:
    """Calculate Shape Complementarity (SC) between two chain groups.

    Shape complementarity measures how well two molecular surfaces fit together,
    following Lawrence & Colman (1993) "Shape Complementarity at Protein/Protein Interfaces".
    Higher SC values (closer to 1.0) indicate better geometric fit between surfaces.
    Typical protein-protein interfaces have SC values between 0.5 and 0.7.

    Args:
        input_file (str): Path to the PDB or mmCIF file
        groups (str): Chain groups specification, e.g., "H,L/A" for chains H,L vs chain A.
            Groups must be disjoint and non-empty, separated by "/"; "H,L/"
            compares H,L against all remaining chains.
        model_num (int, optional): Model serial to analyze (0 for first model). Defaults to 0.
        num_threads (int, optional): Number of threads for parallel calculations (0 for auto). Defaults to 0.

    Returns:
        float: The shape complementarity score (approximately -1 to 1).

    Examples:
        >>> import arpeggia
        >>> sc = arpeggia.sc("antibody_antigen.pdb", groups="H,L/A")
        >>> print(f"SC Score: {sc:.3f}")
    """

class SeqAlignment:
    """Read-only alignment: gapped strings, input spans, scoring and statistics."""
    def format(
        self,
        width: int | None = ...,
        color: Literal["auto", "always", "never"] = ...,
        rulers: bool = ...,
    ) -> str:
        """Format statistics and three alignment rows, optionally with position rulers.

        Args:
            width (int | None): Total visible columns including labels. None
                detects terminal width or uses 80 when unavailable.
            color (str): "auto", "always", or "never". Auto (default) inspects
                sys.stdout and respects NO_COLOR.
            rulers (bool): Show position rulers by default. False hides ticks
                but retains start/end numbers.

        Returns:
            str: Statistics and wrapped alignment text. Color escapes appear
                only when enabled; stored data fields always remain plain.

        Raises:
            ValueError: The width cannot fit labels and one residue, or the
                color policy is unsupported.
        """
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...
    @property
    def reference_name(self) -> str: ...
    @property
    def query_name(self) -> str: ...
    @property
    def reference(self) -> str: ...
    @property
    def query(self) -> str: ...
    @property
    def mode(self) -> str: ...
    @property
    def matrix(self) -> str: ...
    @property
    def gap_open(self) -> float: ...
    @property
    def gap_extend(self) -> float: ...
    @property
    def score(self) -> float: ...
    @property
    def reference_span(self) -> tuple[int, int]: ...
    @property
    def query_span(self) -> tuple[int, int]: ...
    @property
    def aligned_reference(self) -> str: ...
    @property
    def aligned_query(self) -> str: ...
    @property
    def operations(self) -> str: ...
    @property
    def alignment_length(self) -> int: ...
    @property
    def shorter_length(self) -> int: ...
    @property
    def matches(self) -> int: ...
    @property
    def mismatches(self) -> int: ...
    @property
    def paired_residues(self) -> int: ...
    @property
    def gap_residues(self) -> int: ...
    @property
    def gap_runs(self) -> int: ...
    @property
    def edit_distance(self) -> int: ...
    @property
    def identity_alignment(self) -> float | None: ...
    @property
    def identity_shorter(self) -> float: ...
    @property
    def coverage_alignment(self) -> float | None: ...
    @property
    def coverage_shorter(self) -> float: ...

class ResiduePair:
    """Corresponding coordinate residues, retaining author numbers, insertions and names."""
    @property
    def reference_number(self) -> int: ...
    @property
    def reference_insertion(self) -> str: ...
    @property
    def reference_name(self) -> str: ...
    @property
    def query_number(self) -> int: ...
    @property
    def query_insertion(self) -> str: ...
    @property
    def query_name(self) -> str: ...

class ChainAlignment:
    """Reference/query chain IDs, optional sequence alignment and residue pairs."""
    @property
    def reference_chain(self) -> str: ...
    @property
    def query_chain(self) -> str: ...
    @property
    def alignment(self) -> SeqAlignment | None: ...
    @property
    def residue_pairs(self) -> list[ResiduePair]: ...

class RmsdResult:
    """Read-only RMSDs in Å: full evaluation, retained fit core, and initial fit.

    Counts refer to atom pairs. cycles counts rejection inspections, including an
    unchanged final pass; refine_cycles is the requested maximum. chain_alignments
    records correspondence. reference_model/query_model are selected serials.
    """
    @property
    def rmsd(self) -> float: ...
    @property
    def core_rmsd(self) -> float: ...
    @property
    def initial_rmsd(self) -> float: ...
    @property
    def refine_cutoff(self) -> float: ...
    @property
    def initial_fit_atoms(self) -> int: ...
    @property
    def retained_fit_atoms(self) -> int: ...
    @property
    def evaluation_atoms(self) -> int: ...
    @property
    def cycles(self) -> int: ...
    @property
    def refine_cycles(self) -> int: ...
    @property
    def reference_model(self) -> int: ...
    @property
    def query_model(self) -> int: ...
    @property
    def align_seqs(self) -> bool: ...
    @property
    def atoms(self) -> str: ...
    @property
    def superpose_residues(self) -> str: ...
    @property
    def rmsd_residues(self) -> str: ...
    @property
    def chain_alignments(self) -> list[ChainAlignment]: ...

def number_antibody(
    sequence: str,
    *,
    name: str = "Seq001",
    scheme: NumberingScheme | None = None,
    cdr_definition: CdrDefinition = "auto",
    species: GermlineSpecies | Sequence[GermlineSpecies] | None = None,
) -> NumberedAntibody:
    """Number one variable domain and compute separate, tied V/J similarities.

    Args:
        sequence (str): Unaligned amino-acid input; tags and constant tails are
            retained. Only terminal FR1/FR4 truncations are supported.
        name (str): Display name, default Seq001.
        scheme (str | None): IMGT by default; Martin, AHo and Kabat are supported.
            Chothia aliases Martin numbering.
        cdr_definition (str): Auto follows the scheme. An explicit definition
            requires an explicit scheme; Chothia uses distinct consensus regions.
        species (str | Sequence[str] | None): Human, mouse, alpaca, or a sequence
            of these names. None searches all bundled references.

    Returns:
        NumberedAntibody: Read-only residue correspondence, region sequences,
            confidence, diagnostics and V/J matches; missing matches are None.

    Raises:
        ValueError: Invalid input or options.
        RuntimeError: Unsupported domain, truncation, or insertion length.

    Conventions and citations:
    https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md
    IMGT: https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
    Martin/AbM: https://pmc.ncbi.nlm.nih.gov/articles/PMC10939163/
    AHo: https://pubs.rsc.org/en/content/articlehtml/2019/me/c9me00021f
    """

def align_antibodies(
    antibodies: Sequence[NumberedAntibody], *, reference_index: int = 0
) -> AntibodyAlignment:
    """Align compatible antibodies by numbered positions.

    Args:
        antibodies (Sequence[NumberedAntibody]): One or more antibodies in one
            scheme, all heavy or all light. K/L mixtures are supported.
        reference_index (int): Zero-based comparison reference, default 0.

    Returns:
        AntibodyAlignment: Antibodies and gapped rows stored in input order.
            Display places the reference first and hides germlines.

    Raises:
        ValueError: Empty/incompatible inputs or invalid reference index.
    """

class NumberedAntibody:
    """Read-only NumberedAntibody result; constructed by the antibody analysis functions."""

    @property
    def name(self) -> str:
        """Display name."""

    @property
    def input_sequence(self) -> str:
        """Complete normalized input, including unnumbered tails."""

    @property
    def domain_span(self) -> tuple[int, int]:
        """Zero-based, half-open numbered span in the original input."""

    @property
    def chain(self) -> str:
        """Detected H, K or L chain class."""

    @property
    def scheme(self) -> str:
        """Resolved numbering convention."""

    @property
    def cdr_definition(self) -> str:
        """Resolved CDR convention."""

    @property
    def residues(self) -> list[NumberedResidue]:
        """Ordered residue correspondence."""

    @property
    def confidence(self) -> float:
        """Profile confidence heuristic; not a probability of correct numbering."""

    @property
    def matched_profile_positions(self) -> int:
        """Distinct profile positions with matched input residues, excluding insertions."""

    @property
    def diagnostics(self) -> list[str]:
        """Recoverable limitations of this annotation."""

    @property
    def v_match(self) -> GermlineMatch | None:
        """Best qualifying V similarities, or none when reference evidence is insufficient."""

    @property
    def j_match(self) -> GermlineMatch | None:
        """Best qualifying J similarities, or none when reference evidence is insufficient."""

    @property
    def sequence(self) -> str:
        """Derived sequence sequence under this antibody's CDR definition."""

    @property
    def fr1(self) -> str:
        """Derived fr1 sequence under this antibody's CDR definition."""

    @property
    def cdr1(self) -> str:
        """Derived cdr1 sequence under this antibody's CDR definition."""

    @property
    def fr2(self) -> str:
        """Derived fr2 sequence under this antibody's CDR definition."""

    @property
    def cdr2(self) -> str:
        """Derived cdr2 sequence under this antibody's CDR definition."""

    @property
    def fr3(self) -> str:
        """Derived fr3 sequence under this antibody's CDR definition."""

    @property
    def cdr3(self) -> str:
        """Derived cdr3 sequence under this antibody's CDR definition."""

    @property
    def fr4(self) -> str:
        """Derived fr4 sequence under this antibody's CDR definition."""

    def impute(
        self, *, v_reference: str | None = None, j_reference: str | None = None
    ) -> NumberedAntibody:
        """Fill supported terminal FR1/FR4 residues in a new object.

        Args:
            v_reference (str | None): Exact tied V reference ID; None requires
                all tied references to agree on presence and a known residue.
            j_reference (str | None): Equivalent selector for J references.

        Returns:
            NumberedAntibody: Original input/span preserved, with source IDs
                and input_index=None for added residues. Internal gaps and
                unknown input residues remain unchanged.

        Raises:
            ValueError: Unknown reference selector.
        """

    def format(
        self,
        width: int | None = None,
        color: Literal["auto", "always", "never"] = "auto",
        rulers: bool = True,
    ) -> str:
        """Render the input above its combined V/J germlines, with CDR highlighting.

        Each block contains one CDR marker row, input ruler and sequence, germline
        ruler and sequence, then operations relative to the input. CDR1/2/3 bands
        are gray/pink/cyan across every row except operations. Yellow backgrounds mark
        imputed residues and their summary count. All tied reference names appear
        in the summary; only representative V/J sequences are shown.

        Args:
            width (int | None): Total columns including labels; None detects
                terminal width with an 80-column fallback.
            color (str): Auto uses terminal support and NO_COLOR; always/never
                override it. Stored fields remain plain.
            rulers (bool): Show one-based residue positions at every tenth residue.
                Supplied-input coordinates survive imputation; imputed positions
                are blank. Germline counts continue from V through J, ignoring gaps. False
                hides rulers but retains endpoint numbers and CDR markers.

        Returns:
            str: Wrapped sequence rows and separate V/J similarities for the supplied
                input. Outer germline padding is blank; unknown junction hyphens
                are gray with blank operations. A '+' marks a germline insertion,
                '-' a deletion, ':' a positive-BLOSUM62 substitution, and 'x' any
                other mismatch. Matches have blank operations.

        Raises:
            ValueError: Invalid display options or insufficient width.
        """

class AntibodyAlignment:
    """Read-only AntibodyAlignment result; constructed by the antibody analysis functions."""

    @property
    def antibodies(self) -> list[NumberedAntibody]:
        """Original numbered antibodies, retaining their CDR definitions and germlines."""

    @property
    def positions(self) -> list[NumberedPosition]:
        """Ordered union of numbered positions."""

    @property
    def aligned_sequences(self) -> list[str]:
        """One plain gapped sequence per input, in original row order."""

    @property
    def reference_index(self) -> int:
        """Zero-based input row used as the display comparison reference."""

    def format(
        self,
        width: int | None = None,
        color: Literal["auto", "always", "never"] = "auto",
        rulers: bool = True,
        *,
        reference_index: int | None = None,
    ) -> str:
        """Render comparisons with the selected reference first and germlines hidden.

        The selected reference defines the shared CDR bands and the convention
        summary. Each row retains its own original input coordinates. Imputed
        residues have yellow backgrounds; the summary count totals all antibodies.

        Args:
            width (int | None): Total columns including labels; None detects
                terminal width with an 80-column fallback.
            color (str): Auto uses terminal support and NO_COLOR; always/never
                override it. Stored fields remain plain.
            rulers (bool): Show one-based original input positions at every tenth
                residue. Imputed positions are blank. False hides rulers while
                retaining endpoint numbers and the single CDR marker row per block.
            reference_index (int | None): Zero-based reference override for
                this rendering, preserving stored row order and reference.

        Returns:
            str: Styled or plain alignment text.

        Raises:
            ValueError: Invalid display options or insufficient width.
        """

class NumberedPosition:
    """Read-only NumberedPosition result; constructed by the antibody analysis functions."""

    @property
    def number(self) -> int:
        """Position number in the selected scheme."""

    @property
    def insertion(self) -> str | None:
        """Insertion letter, when present."""

class NumberedResidue:
    """Read-only NumberedResidue result; constructed by the antibody analysis functions."""

    @property
    def position(self) -> NumberedPosition:
        """Numbered position."""

    @property
    def amino_acid(self) -> str:
        """Supplied amino-acid symbol."""

    @property
    def input_index(self) -> int | None:
        """Zero-based index in the original input; absent for imputed residues."""

    @property
    def region(self) -> str:
        """FR1, CDR1, FR2, CDR2, FR3, CDR3 or FR4."""

    @property
    def imputed_from(self) -> list[str]:
        """Source reference IDs for an imputed residue; empty for supplied residues."""

class GermlineReference:
    """Read-only GermlineReference result; constructed by the antibody analysis functions."""

    @property
    def id(self) -> str:
        """Stable reference selector: accession, gene/allele, species and source span."""

    @property
    def species(self) -> str:
        """Full source species name, including any strain/subspecies suffix."""

    @property
    def gene(self) -> str:
        """Gene name without allele suffix."""

    @property
    def allele(self) -> str:
        """Allele suffix."""

    @property
    def accession(self) -> str:
        """Source sequence accession."""

class GermlineHit:
    """Read-only GermlineHit result; constructed by the antibody analysis functions."""

    @property
    def references(self) -> list[GermlineReference]:
        """All references sharing this sequence and coverage, sorted by source identity."""

    @property
    def alignment(self) -> SeqAlignment:
        """Optimal local reference-to-input alignment with BLOSUM62 and gap costs 10/0.5."""

    @property
    def query_input_start(self) -> int:
        """Input index corresponding to query index zero in the segment alignment."""

    @property
    def imgt_positions(self) -> list[NumberedPosition | None]:
        """Internal IMGT positions per reference residue. J junction residues have none."""

    @property
    def imgt_span(self) -> tuple[int, int]:
        """One-based, half-open span of available IMGT reference coverage."""

    @property
    def known_pairs(self) -> int:
        """Nongap pairs where both residues are among the 20 standard amino acids."""

    @property
    def known_matches(self) -> int:
        """Identical standard-amino-acid pairs."""

    @property
    def known_fr4_pairs(self) -> int:
        """Known pairs in the input's IMGT FR4 (118–128)."""

    @property
    def known_identity(self) -> float:
        """Known matches divided by known pairs."""

    @property
    def reference_coverage(self) -> float:
        """Known pairs divided by all standard residues in the reference segment."""

    @property
    def query_coverage(self) -> float:
        """Known pairs divided by all standard residues in the input segment."""

class GermlineMatch:
    """Read-only GermlineMatch result; constructed by the antibody analysis functions."""

    @property
    def score(self) -> float:
        """Maximum qualifying local alignment score."""

    @property
    def hits(self) -> list[GermlineHit]:
        """Tied sequence/coverage groups; the first is the display representative."""
