"""Arpeggia: Protein structure analysis library.

This module provides functions for analyzing protein structures from PDB and mmCIF files.
It includes functionality for:
- Detecting atomic and ring contacts between protein chains
- Calculating solvent accessible surface area (SASA) for atoms, residues, and chains
- Calculating relative SASA (RSA) for residues normalized by Tien et al. (2013) MaxASA values
- Calculating Spatial Aggregation Propensity (SAP) scores for predicting aggregation-prone regions
- Calculating buried surface area (dSASA) at interfaces between chain groups
- Calculating Shape Complementarity (SC) between chain groups
- Extracting coordinate-observed and declared polymer sequences

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
    >>> # Calculate SASA for specific chains only
    >>> sasa_ab = arpeggia.sasa("structure.pdb", level="residue", chains="A,B")
    >>> print(f"Calculated SASA for {len(sasa_ab)} residues in chains A and B")
    >>>
    >>> # Calculate relative SASA (RSA) for each residue
    >>> rsa_df = arpeggia.relative_sasa("structure.pdb")
    >>> print(f"Calculated RSA for {len(rsa_df)} residues")
    >>>
    >>> # Calculate Spatial Aggregation Propensity (SAP) scores
    >>> sap_df = arpeggia.sap_score("antibody.pdb", level="residue")
    >>> print(f"Calculated SAP for {len(sap_df)} residues")
    >>>
    >>> # Calculate SAP for specific chains (e.g., antibody H and L chains)
    >>> sap_hl = arpeggia.sap_score("antibody.pdb", chains="H,L")
    >>> print(f"Calculated SAP for {len(sap_hl)} residues in H and L chains")
    >>>
    >>> # Calculate buried surface area at interface
    >>> bsa = arpeggia.dsasa("structure.pdb", groups="A,B/C,D")
    >>> print(f"Buried surface area: {bsa:.2f} Å²")
    >>> components = arpeggia.dsasa_components("structure.pdb", groups="A,B/C,D")
    >>>
    >>> # Calculate shape complementarity at interface
    >>> sc = arpeggia.sc("antibody_antigen.pdb", groups="H,L/A")
    >>> print(f"Shape complementarity: {sc:.3f}")
    >>>
    >>> # Extract sequences
    >>> sequences = arpeggia.seq("structure.pdb")
    >>> for chain_id, seq in sequences:
    ...     print(f"Chain {chain_id}: {seq}")
    >>> declared_sequences = arpeggia.seqres("structure.pdb")
"""

from importlib.metadata import version

from arpeggia._contract import EXPORTED_FUNCTIONS
from arpeggia.arpeggia import (  # noqa: F401
    cluster_structs,
    contacts,
    dsasa,
    dsasa_components,
    pairwise_rmsd,
    relative_sasa,
    rmsd,
    sap_score,
    sasa,
    sc,
    seq,
    seqres,
)

__version__ = version("arpeggia")
__all__ = list(EXPORTED_FUNCTIONS)
