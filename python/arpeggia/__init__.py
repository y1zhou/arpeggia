"""Arpeggia: Protein structure analysis library.

This module provides functions for analyzing protein structures from PDB and mmCIF files.
It includes functionality for:
- Detecting atomic and ring contacts between protein chains
- Calculating solvent accessible surface area (SASA) for atoms, residues, and chains
- Extracting protein sequences from structures

The module is built on Rust using PyO3 and returns results as Polars DataFrames for
efficient data manipulation.

Example:
    >>> import arpeggia
    >>> # Analyze contacts in a protein structure
    >>> contacts_df = arpeggia.contacts("structure.pdb", groups="/", vdw_comp=0.1, ignore_zero_occupancy=False)
    >>> print(f"Found {len(contacts_df)} contacts")
    >>>
    >>> # Calculate SASA for all atoms
    >>> sasa_df = arpeggia.sasa("structure.pdb", probe_radius=1.4)
    >>> print(f"Calculated SASA for {len(sasa_df)} atoms")
    >>>
    >>> # Calculate SASA aggregated by residue
    >>> residue_sasa_df = arpeggia.residue_sasa("structure.pdb", probe_radius=1.4)
    >>> print(f"Calculated SASA for {len(residue_sasa_df)} residues")
    >>>
    >>> # Calculate SASA aggregated by chain
    >>> chain_sasa_df = arpeggia.chain_sasa("structure.pdb", probe_radius=1.4)
    >>> print(f"Calculated SASA for {len(chain_sasa_df)} chains")
    >>>
    >>> # Extract sequences
    >>> sequences = arpeggia.pdb2seq("structure.pdb")
    >>> for chain_id, seq in sequences.items():
    ...     print(f"Chain {chain_id}: {seq}")
"""

from importlib.metadata import version

from arpeggia.arpeggia import chain_sasa, contacts, pdb2seq, residue_sasa, sasa

__version__ = version("arpeggia")
__all__ = ["contacts", "sasa", "residue_sasa", "chain_sasa", "pdb2seq"]
