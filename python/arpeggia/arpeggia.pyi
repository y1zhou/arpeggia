"""Type stubs for arpeggia Rust module."""

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
    probe_radius: float = 1.4,
    n_points: int = 100,
    model_num: int = 0,
    num_threads: int = 1,
) -> pl.DataFrame:
    """Load a PDB or mmCIF file and calculate solvent accessible surface area (SASA) for each atom.

    Args:
        input_file: Path to the PDB or mmCIF file
        probe_radius: Probe radius in Ångströms. Defaults to 1.4.
        n_points: Number of points for surface calculation. Defaults to 100.
        model_num: Model number to analyze (0 for first model). Defaults to 0.
        num_threads: Number of threads for parallel processing (0 for all cores). Defaults to 1.

    Returns:
        A Polars DataFrame with SASA values for each atom with columns:
        - atomi, sasa
        - chain, resn, resi, insertion, altloc, atomn
    """
    ...

def residue_sasa(
    input_file: str,
    probe_radius: float = 1.4,
    n_points: int = 100,
    model_num: int = 0,
    num_threads: int = 1,
) -> pl.DataFrame:
    """Load a PDB or mmCIF file and calculate solvent accessible surface area (SASA) aggregated by residue.

    Args:
        input_file: Path to the PDB or mmCIF file
        probe_radius: Probe radius in Ångströms. Defaults to 1.4.
        n_points: Number of points for surface calculation. Defaults to 100.
        model_num: Model number to analyze (0 for first model). Defaults to 0.
        num_threads: Number of threads for parallel processing (0 for all cores). Defaults to 1.

    Returns:
        A Polars DataFrame with SASA values for each residue with columns:
        - chain, resn, resi, insertion, altloc, sasa, is_polar
    """
    ...

def chain_sasa(
    input_file: str,
    probe_radius: float = 1.4,
    n_points: int = 100,
    model_num: int = 0,
    num_threads: int = 1,
) -> pl.DataFrame:
    """Load a PDB or mmCIF file and calculate solvent accessible surface area (SASA) aggregated by chain.

    Args:
        input_file: Path to the PDB or mmCIF file
        probe_radius: Probe radius in Ångströms. Defaults to 1.4.
        n_points: Number of points for surface calculation. Defaults to 100.
        model_num: Model number to analyze (0 for first model). Defaults to 0.
        num_threads: Number of threads for parallel processing (0 for all cores). Defaults to 1.

    Returns:
        A Polars DataFrame with SASA values for each chain with columns:
        - chain, sasa
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
