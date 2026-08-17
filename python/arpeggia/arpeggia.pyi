"""Type stubs for the arpeggia Rust module."""

import polars as pl

from ._contract import (
    DsasaComponents,
    ProtonationMode,
    SapLevel,
    SasaLevel,
    SequenceList,
)

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
def seq(input_file: str) -> SequenceList: ...
def seqres(input_file: str) -> SequenceList: ...
def sc(
    input_file: str,
    groups: str,
    model_num: int = ...,
    num_threads: int = ...,
) -> float: ...
