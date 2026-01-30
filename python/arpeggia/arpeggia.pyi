"""Type stubs for arpeggia Rust module."""

from typing import Literal

import polars as pl

def contacts(
    input_file: str,
    groups: str = "/",
    vdw_comp: float = 0.1,
    dist_cutoff: float = 6.5,
    ignore_zero_occupancy: bool = False,
) -> pl.DataFrame:
    """Load a PDB or mmCIF file and calculate atomic and ring contacts.

    Args:
        input_file: Path to the PDB or mmCIF file
        groups: Chain groups specification. Defaults to "/" (all-to-all).
            Examples: "A,B/C,D" for chains A,B vs C,D; "A/" for chain A vs all others.
        vdw_comp: VdW radii compensation factor. Defaults to 0.1.
        dist_cutoff: Distance cutoff for neighbor searches in Ångströms. Defaults to 6.5.
        ignore_zero_occupancy: If True, ignore atoms with zero occupancy. Defaults to False.

    Returns:
        A Polars DataFrame containing all identified contacts with columns:
        - model, interaction, distance
        - from_chain, from_resn, from_resi, from_insertion, from_altloc, from_atomn, from_atomi
        - to_chain, to_resn, to_resi, to_insertion, to_altloc, to_atomn, to_atomi
        - sc_centroid_dist, sc_dihedral, sc_centroid_angle
    """
    ...

def sasa(
    input_file: str,
    level: Literal["atom", "residue", "chain"] = "atom",
    probe_radius: float = 1.4,
    n_points: int = 100,
    model_num: int = 0,
    num_threads: int = 1,
    chains: str = "",
) -> pl.DataFrame:
    """Load a PDB or mmCIF file and calculate solvent accessible surface area (SASA).

    Args:
        input_file: Path to the PDB or mmCIF file
        level: Aggregation level for SASA calculation. Options:
            - "atom": Calculate SASA for each atom (default)
            - "residue": Aggregate SASA by residue
            - "chain": Aggregate SASA by chain
        probe_radius: Probe radius in Ångströms. Defaults to 1.4.
        n_points: Number of points for surface calculation. Defaults to 100.
        model_num: Model number to analyze (0 for first model). Defaults to 0.
        num_threads: Number of threads for parallel processing (0 for all cores). Defaults to 1.
        chains: Comma-separated chain IDs to include (e.g., "A,B,C").
            If empty, includes all chains. Defaults to "".

    Returns:
        A Polars DataFrame with SASA values. Columns depend on the level:
        - atom: atomi, sasa, chain, resn, resi, insertion, altloc, atomn
        - residue: chain, resn, resi, insertion, altloc, sasa, is_polar
        - chain: chain, sasa
    """
    ...

def dsasa(
    input_file: str,
    groups: str,
    probe_radius: float = 1.4,
    n_points: int = 100,
    model_num: int = 0,
    num_threads: int = 1,
) -> float:
    """Load a PDB or mmCIF file and calculate buried surface area at the interface.

    The buried surface area (dSASA) is calculated as:
    dSASA = SASA_group1 + SASA_group2 - SASA_complex

    Args:
        input_file: Path to the PDB or mmCIF file
        groups: Chain groups specification for interface calculation.
            Format: "A,B/C,D" where chains A,B form one side and C,D form the other side.
        probe_radius: Probe radius in Ångströms. Defaults to 1.4.
        n_points: Number of points for surface calculation. Defaults to 100.
        model_num: Model number to analyze (0 for first model). Defaults to 0.
        num_threads: Number of threads for parallel processing (0 for all cores). Defaults to 1.

    Returns:
        The buried surface area at the interface in square Ångströms.
    """
    ...

def relative_sasa(
    input_file: str,
    probe_radius: float = 1.4,
    n_points: int = 100,
    model_num: int = 0,
    num_threads: int = 1,
    chains: str = "",
) -> pl.DataFrame:
    """Load a PDB or mmCIF file and calculate relative SASA (RSA) for each residue.

    RSA is calculated as the ratio of observed SASA to the maximum possible SASA
    for each amino acid type, based on Tien et al. (2013) theoretical values.

    Args:
        input_file: Path to the PDB or mmCIF file
        probe_radius: Probe radius in Ångströms. Defaults to 1.4.
        n_points: Number of points for surface calculation. Defaults to 100.
        model_num: Model number to analyze (0 for first model). Defaults to 0.
        num_threads: Number of threads for parallel processing (0 for all cores). Defaults to 1.
        chains: Comma-separated chain IDs to include (e.g., "A,B,C").
            If empty, includes all chains. Defaults to "".

    Returns:
        A Polars DataFrame with relative SASA values for each residue with columns:
        - chain, resn, resi, insertion, altloc, sasa, is_polar, relative_sasa
    """
    ...

def sap_score(
    input_file: str,
    level: Literal["atom", "residue"] = "residue",
    probe_radius: float = 1.4,
    n_points: int = 100,
    model_num: int = 0,
    sap_radius: float = 5.0,
    num_threads: int = 1,
    chains: str = "",
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
    - Hydrophobicity uses the Black & Mould (1991) scale, normalized so glycine = 0
    - SASA is the side-chain solvent accessible surface area
    - SASA_max is the maximum SASA for that residue type

    Args:
        input_file: Path to the PDB or mmCIF file
        level: Aggregation level for SAP calculation. Options:
            - "atom": Calculate SAP for each atom
            - "residue": Aggregate SAP by residue (default)
        probe_radius: Probe radius in Ångströms for SASA calculation. Defaults to 1.4.
        n_points: Number of points for SASA surface calculation. Defaults to 100.
        model_num: Model number to analyze (0 for first model). Defaults to 0.
        sap_radius: Radius in Ångströms for neighbor search. Defaults to 5.0.
        num_threads: Number of threads for parallel processing (0 for all cores). Defaults to 1.
        chains: Comma-separated chain IDs to include (e.g., "H,L").
            If empty, includes all chains. Defaults to "".

    Returns:
        A Polars DataFrame with SAP scores. Columns depend on the level:
        - atom: chain, resn, resi, insertion, atomn, atomi, sasa, sap_score
        - residue: chain, resn, resi, insertion, sc_sasa, sap_score
    """
    ...

def pdb2seq(input_file: str) -> dict[str, str]:
    """Load a PDB or mmCIF file and extract sequences for all chains.

    Args:
        input_file: Path to the PDB or mmCIF file

    Returns:
        A dictionary mapping chain IDs to their sequences
    """
    ...
