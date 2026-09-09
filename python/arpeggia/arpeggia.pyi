"""Type stubs for the arpeggia Rust module."""

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
) -> SeqAlignment: ...
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
) -> RmsdResult: ...
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
) -> pl.DataFrame: ...
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
) -> pl.DataFrame: ...
def contacts(
    input_file: str,
    groups: str = ...,
    vdw_comp: float = ...,
    dist_cutoff: float = ...,
    ignore_zero_occupancy: bool = ...,
    protonation: ProtonationMode = ...,
    ph: float = ...,
    num_threads: int = ...,
) -> pl.DataFrame: ...
def sasa(
    input_file: str,
    level: SasaLevel = ...,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    chains: str = ...,
    num_threads: int = ...,
) -> pl.DataFrame: ...
def dsasa(
    input_file: str,
    groups: str,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    num_threads: int = ...,
) -> float: ...
def dsasa_components(
    input_file: str,
    groups: str,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    num_threads: int = ...,
) -> DsasaComponents: ...
def relative_sasa(
    input_file: str,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    chains: str = ...,
    num_threads: int = ...,
) -> pl.DataFrame: ...
def sap_score(
    input_file: str,
    level: SapLevel = ...,
    probe_radius: float = ...,
    n_points: int = ...,
    model_num: int = ...,
    sap_radius: float = ...,
    chains: str = ...,
    num_threads: int = ...,
) -> pl.DataFrame: ...
def seq(input_file: str, model_num: int = ...) -> SequenceList: ...
def seqres(input_file: str) -> SequenceList: ...
def sc(
    input_file: str,
    groups: str,
    model_num: int = ...,
    num_threads: int = ...,
) -> float: ...

class SeqAlignment:
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
    def columns(self) -> list[tuple[int | None, int | None]]: ...
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
    @property
    def reference_chain(self) -> str: ...
    @property
    def mobile_chain(self) -> str: ...
    @property
    def alignment(self) -> SeqAlignment | None: ...
    @property
    def residue_pairs(self) -> list[ResiduePair]: ...

class RmsdResult:
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
