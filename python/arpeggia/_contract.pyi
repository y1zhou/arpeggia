"""Type-level Python exposure contract for arpeggia."""

from typing import Final, Literal, TypeAlias

SasaLevel: TypeAlias = Literal["atom", "residue", "chain"]
SapLevel: TypeAlias = Literal["atom", "residue"]
SequenceList: TypeAlias = list[tuple[str, str]]

EXPORTED_FUNCTIONS: Final[tuple[str, ...]]
SASA_LEVELS: Final[tuple[SasaLevel, ...]]
SAP_LEVELS: Final[tuple[SapLevel, ...]]
DEFAULTS: Final[dict[str, dict[str, object]]]
CONTACT_COLUMNS: Final[tuple[str, ...]]
SASA_COLUMNS: Final[dict[SasaLevel, tuple[str, ...]]]
RELATIVE_SASA_COLUMNS: Final[tuple[str, ...]]
SAP_COLUMNS: Final[dict[SapLevel, tuple[str, ...]]]
SEQUENCE_RETURN: Final[str]
