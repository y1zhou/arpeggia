"""Type stubs for the arpeggia Rust module."""

from typing import Literal

import polars as pl

from ._contract import (
    AlignmentMode,
    AtomSubset,
    ClusteringMethod,
    DsasaComponents,
    ProtonationMode,
    SapLevel,
    SasaLevel,
    SequenceList,
)

def align_seqs(
    reference: str,
    mobile: str,
    mode: AlignmentMode = ...,
    gap_open: float = ...,
    gap_extend: float = ...,
) -> SeqAlignment:
    """Align two non-empty unaligned protein strings with BLOSUM62.

    mode is global (both complete), local, or semi-global (mobile complete,
    reference tails free). Gap costs are positive, at most two decimals, with
    open >= extend; a gap costs open + (length-1)*extend. Returns read-only
    gapped strings, operations, zero-based half-open spans, counts and ratios.
    edit_distance compares full inputs independently of the protein score.
    Example: align_seqs("GGACDEFGHIKGG", "ACDEFGHIK", mode="semi-global").
    """

def rmsd(
    reference: str,
    mobile: str,
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
    """Fit two PDB/mmCIF structures and return separate evaluation/core RMSDs in Å.

    superpose_residues sets the fit; rmsd_residues sets evaluation under that
    transform. Each independently defaults to all eligible residues. Syntax:
    "A:1-100,A:110-120,B" (inclusive author-number ranges or whole chains),
    "A:-5--1" (negative numbers), "A:10A" (insertion), "A:10" (all insertions).
    Repeat chain IDs in comma clauses. model_num=0 selects each first model.
    atoms: ca, backbone (N/CA/C/O/OXT), heavy, all.

    align_seqs=True establishes observed-chain correspondence before selection,
    using reference numbering. chain_map={"A": "H"} explicitly maps every
    selected reference chain to a unique mobile chain; otherwise infer a unique
    maximum-score assignment. alignment_mode/gap costs follow align_seqs().
    refine_cycles=0 means no rejection; positive values reject/refit at most N
    times using refine_cutoff (default 2) times current fit RMSD.

    result.rmsd includes all mapped evaluation pairs, even rejected fitting
    pairs; core_rmsd uses retained fitting pairs. Degenerate fits fail.
    Example: rmsd("ref.cif", "mob.cif", superpose_residues="A", rmsd_residues="B").
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
    """Return id_1, id_2, rmsd (Å) columns for every unordered structure pair.

    input is a non-recursive PDB/mmCIF directory or CSV/Parquet/NDJSON manifest;
    id_col/path_col name manifest columns, whose paths resolve against the
    manifest. At least two structures and exact atom correspondence are required.
    Selections use rmsd() syntax and independently default to all; alignment and
    refinement are unavailable. model_num=0 selects first models, num_threads=0
    uses available processors, bypass_mem_check skips heuristic memory checks.
    Example: pairwise_rmsd("structures/", superpose_residues="A", rmsd_residues="B").
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
    """Return id, cluster_id, medoid_id, rmsd_to_medoid (Å) columns.

    Supply exactly one input directory/manifest or complete pairwise_rmsd table
    (id_1, id_2, rmsd). At least three structures are required. Use a fixed
    num_clusters from 1 through n, or max_clusters from 2 through n-1 for automatic
    selection; fixed count wins with a warning. method is k-medoids.
    max_iterations bounds convergence. Structure options follow pairwise_rmsd()
    and do not change a supplied pair table. Cluster IDs are zero-based.
    Example: cluster_structs(input="structures/", num_clusters=3).
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
    """Return atomic/aromatic contacts as a Polars DataFrame.

    input_file is PDB/mmCIF. groups="/" means all pairs, "A,B/C" compares groups,
    "A/" compares A to remaining chains. vdw_comp and dist_cutoff are in Å.
    ignore_zero_occupancy removes occupancy-zero atoms. protonation is
    all-charged, heuristic (uses ph), or explicit-only. num_threads=0 is automatic.
    Example: contacts("complex.cif", groups="H,L/A").
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
    """Return solvent accessible areas in Å² per atom, residue, or chain.

    input_file is PDB/mmCIF; level is atom, residue, or chain. probe_radius is
    in Å; n_points is the positive sample count per sphere. model_num=0 selects
    the first model. chains="" means all; "A,B" filters chains. num_threads=0
    uses available processors. Example: sasa("complex.cif", level="residue").
    """

def dsasa(
    input_file: str,
    groups: str,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    num_threads: int = ...,
) -> float:
    """Return two-sided buried area in Å²: SASA(group1)+SASA(group2)-SASA(complex).

    PDB/mmCIF input; groups must be disjoint and non-empty, e.g. "H,L/A" or
    "A/" (A vs rest). probe_radius is in Å; n_points is samples per sphere;
    model_num=0 selects the first model and num_threads=0 is automatic.
    Example: dsasa("complex.cif", groups="H,L/A").
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

    Arguments and group restrictions match dsasa(); partitions sum to total.
    Example: dsasa_components("complex.cif", groups="H,L/A").
    """

def relative_sasa(
    input_file: str,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    chains: str = ...,
    num_threads: int = ...,
) -> pl.DataFrame:
    """Return per-residue SASA (Å²) and dimensionless relative_sasa ratios.

    Ratios use Tien reference maxima, not percentages. PDB/mmCIF input;
    probe_radius, n_points, model_num, chains and num_threads follow sasa().
    Example: relative_sasa("complex.cif", chains="A").
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
    """Return Spatial Aggregation Propensity per atom or residue (default).

    PDB/mmCIF input; probe_radius and sap_radius (neighbor search) are in Å.
    n_points, model_num, chains and num_threads follow sasa(). SASA columns
    contain Å²; sap_score and relative_sc_sasa are dimensionless.
    Example: sap_score("antibody.cif", chains="H,L").
    """

def seq(input_file: str, model_num: int = ...) -> SequenceList:
    """Return (chain ID, observed protein sequence) tuples from PDB/mmCIF.

    model_num is a model serial, or 0 for the first model. Residues without
    coordinates are absent. Example: seq("structure.cif", model_num=0).
    """

def seqres(input_file: str) -> SequenceList:
    """Return (chain ID, declared sequence) tuples from PDB SEQRES/mmCIF metadata.

    Includes residues without coordinates; use seq() for observed proteins.
    Example: seqres("structure.cif").
    """

def sc(
    input_file: str,
    groups: str,
    model_num: int = ...,
    num_threads: int = ...,
) -> float:
    """Return dimensionless shape complementarity (approximately -1 to 1).

    PDB/mmCIF input; groups must be disjoint and non-empty, e.g. "H,L/A" or
    "H,L/" (H,L vs remaining chains). Higher scores mean better geometric fit.
    model_num=0 selects the first model; num_threads=0 uses available processors.
    Example: sc("complex.cif", groups="H,L/A").
    """

class SeqAlignment:
    """Read-only alignment: gapped strings, input spans, scoring and statistics."""
    def format(
        self,
        width: int | None = ...,
        color: Literal["auto", "always", "never"] = ...,
        rulers: bool = ...,
    ) -> str:
        """Wrap to total width; auto color follows the terminal; rulers=False hides ticks."""
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...
    @property
    def reference(self) -> str: ...
    @property
    def mobile(self) -> str: ...
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
    def mobile_span(self) -> tuple[int, int]: ...
    @property
    def aligned_reference(self) -> str: ...
    @property
    def aligned_mobile(self) -> str: ...
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
    def mobile_number(self) -> int: ...
    @property
    def mobile_insertion(self) -> str: ...
    @property
    def mobile_name(self) -> str: ...

class ChainAlignment:
    """Reference/mobile chain IDs, optional sequence alignment and residue pairs."""
    @property
    def reference_chain(self) -> str: ...
    @property
    def mobile_chain(self) -> str: ...
    @property
    def alignment(self) -> SeqAlignment | None: ...
    @property
    def residue_pairs(self) -> list[ResiduePair]: ...

class RmsdResult:
    """Read-only RMSDs in Å: full evaluation, retained fit core, and initial fit.

    Counts refer to atom pairs. cycles counts rejection inspections, including an
    unchanged final pass; refine_cycles is the requested maximum. chain_alignments
    records correspondence. reference_model/mobile_model are selected serials.
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
    def mobile_model(self) -> int: ...
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
