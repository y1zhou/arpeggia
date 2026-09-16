"""Protein structure analysis, sequence alignment and antibody numbering.

Tabular analyses return Polars DataFrames; sequence alignment, RMSD and antibody
APIs return structured result objects. Scientific diagnostics emit UserWarning.
See https://github.com/y1zhou/arpeggia/blob/master/docs/examples.md for Python
and CLI examples and the individual functions for arguments and return types.
"""

from importlib.metadata import version

from arpeggia._contract import EXPORTED_CLASSES, EXPORTED_FUNCTIONS
from arpeggia.arpeggia import (  # noqa: F401
    AntibodyAlignment,
    ChainAlignment,
    GermlineHit,
    GermlineMatch,
    GermlineReference,
    NumberedAntibody,
    NumberedPosition,
    NumberedResidue,
    ResiduePair,
    RmsdResult,
    SeqAlignment,
    align_antibodies,
    align_seqs,
    cluster_structs,
    contacts,
    dsasa,
    dsasa_components,
    number_antibody,
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
__all__ = list(EXPORTED_FUNCTIONS + EXPORTED_CLASSES)
