"""Public stub usage checked by ty; this function is not executed."""

import arpeggia
import polars as pl
from arpeggia._contract import AtomSubset, ProtonationMode, SapLevel, SasaLevel


def public_api_types(
    path: str,
    protonation: ProtonationMode,
    sasa_level: SasaLevel,
    sap_level: SapLevel,
    atoms: AtomSubset,
) -> tuple[
    pl.DataFrame, pl.DataFrame, pl.DataFrame, arpeggia.RmsdResult, list[tuple[str, str]]
]:
    """Check accepted selections and concrete returns through public imports."""
    return (
        arpeggia.contacts(path, protonation=protonation),
        arpeggia.sasa(path, level=sasa_level),
        arpeggia.sap_score(path, level=sap_level),
        arpeggia.rmsd(path, path, atoms=atoms),
        arpeggia.seq(path, model_num=0),
    )


def alignment_api_types() -> tuple[arpeggia.SeqAlignment, float | None]:
    """Check the typed alignment object and its empty-alignment ratio."""
    result = arpeggia.align_seqs("ACDE", "ACD", mode="semi-global")
    return result, result.identity_alignment


def antibody_api_types(
    sequence: str,
) -> tuple[arpeggia.NumberedAntibody, arpeggia.AntibodyAlignment, str]:
    """Check antibody classes, region views, imputation and numbered alignments."""
    antibody = arpeggia.number_antibody(
        sequence, scheme="imgt", species=["human", "alpaca", "rat", "rabbit"]
    )
    completed = antibody.impute()
    comparison = arpeggia.align_antibodies([antibody, completed], reference_index=1)
    return antibody, comparison, comparison.format(color="never", reference_index=0)
