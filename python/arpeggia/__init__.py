"""Arpeggia: Protein structure analysis library.

This module provides functions for analyzing protein structures from PDB and mmCIF files.
It includes functionality for:
- Detecting atomic and ring contacts between protein chains
- Calculating solvent accessible surface area (SASA) for atoms, residues, and chains
- Calculating relative SASA (RSA) for residues normalized by Tien et al. (2013) MaxASA values
- Calculating buried surface area (dSASA) at interfaces between chain groups
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
    >>> sasa_df = arpeggia.sasa("structure.pdb", level="atom")
    >>> print(f"Calculated SASA for {len(sasa_df)} atoms")
    >>>
    >>> # Calculate SASA aggregated by residue
    >>> residue_sasa_df = arpeggia.sasa("structure.pdb", level="residue")
    >>> print(f"Calculated SASA for {len(residue_sasa_df)} residues")
    >>>
    >>> # Calculate SASA aggregated by chain
    >>> chain_sasa_df = arpeggia.sasa("structure.pdb", level="chain")
    >>> print(f"Calculated SASA for {len(chain_sasa_df)} chains")
    >>>
    >>> # Calculate relative SASA (RSA) for each residue
    >>> rsa_df = arpeggia.relative_sasa("structure.pdb")
    >>> print(f"Calculated RSA for {len(rsa_df)} residues")
    >>>
    >>> # Calculate buried surface area at interface
    >>> bsa = arpeggia.dsasa("structure.pdb", groups="A,B/C,D")
    >>> print(f"Buried surface area: {bsa:.2f} Å²")
    >>>
    >>> # Extract sequences
    >>> sequences = arpeggia.pdb2seq("structure.pdb")
    >>> for chain_id, seq in sequences.items():
    ...     print(f"Chain {chain_id}: {seq}")
"""

from importlib.metadata import version

from arpeggia.arpeggia import contacts, dsasa, pdb2seq, relative_sasa, sasa

__version__ = version("arpeggia")
__all__ = ["contacts", "sasa", "relative_sasa", "dsasa", "pdb2seq"]
