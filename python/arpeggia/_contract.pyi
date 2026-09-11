"""Type-level Python exposure contract for arpeggia."""

from typing import Final, Literal, TypeAlias

SasaLevel: TypeAlias = Literal["atom", "residue", "chain"]
SapLevel: TypeAlias = Literal["atom", "residue"]
ProtonationMode: TypeAlias = Literal["all-charged", "heuristic", "explicit-only"]
NumberingScheme: TypeAlias = Literal["imgt", "martin", "chothia", "aho", "kabat"]
CdrDefinition: TypeAlias = Literal["auto", "imgt", "martin", "chothia", "aho", "kabat"]
GermlineSpecies: TypeAlias = Literal["human", "mouse", "alpaca", "rat", "rabbit"]
AlignmentMode: TypeAlias = Literal["global", "local", "semi-global"]
AtomSubset: TypeAlias = Literal["ca", "backbone", "heavy", "all"]
ClusteringMethod: TypeAlias = Literal["k-medoids"]
SequenceList: TypeAlias = list[tuple[str, str]]
DsasaComponents: TypeAlias = tuple[float, float, float, float]

EXPORTED_CLASSES: Final[tuple[str, ...]]
EXPORTED_FUNCTIONS: Final[tuple[str, ...]]
CONTACT_COLUMNS: Final[tuple[str, ...]]
SASA_COLUMNS: Final[dict[SasaLevel, tuple[str, ...]]]
PAIRWISE_RMSD_COLUMNS: Final[tuple[str, str, str]]
CLUSTER_COLUMNS: Final[tuple[str, str, str, str]]
