"""Arpeggia: Protein structure analysis library.

This module provides functions for analyzing protein structures from PDB and mmCIF files.
It includes functionality for:
- Detecting atomic and ring contacts between protein chains
- Calculating solvent accessible surface area (SASA) for atoms
- Extracting protein sequences from structures

The module is built on Rust using PyO3 and returns results as Polars DataFrames for
efficient data manipulation.

Example:
    >>> import arpeggia
    >>> # Analyze contacts in a protein structure
    >>> contacts_df = arpeggia.contacts("structure.pdb", groups="/", vdw_comp=0.1)
    >>> print(f"Found {len(contacts_df)} contacts")
    >>>
    >>> # Calculate SASA for all atoms
    >>> sasa_df = arpeggia.sasa("structure.pdb", probe_radius=1.4)
    >>> print(f"Calculated SASA for {len(sasa_df)} atoms")
    >>>
    >>> # Extract sequences
    >>> sequences = arpeggia.pdb2seq("structure.pdb")
    >>> for chain_id, seq in sequences.items():
    ...     print(f"Chain {chain_id}: {seq}")
"""

from arpeggia.arpeggia import contacts, pdb2seq, sasa

__version__ = "0.5.0"
__all__ = ["contacts", "sasa", "pdb2seq"]
