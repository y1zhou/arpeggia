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
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, float, list[tuple[str, str]]]:
    """Check accepted selections and concrete returns through public imports."""
    return (
        arpeggia.contacts(path, protonation=protonation),
        arpeggia.sasa(path, level=sasa_level),
        arpeggia.sap_score(path, level=sap_level),
        arpeggia.rmsd(path, path, atoms=atoms),
        arpeggia.seq(path, model_num=0),
    )
