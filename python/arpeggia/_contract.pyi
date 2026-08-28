"""Type-level Python exposure contract for arpeggia."""

from typing import Final, Literal, TypeAlias

SasaLevel: TypeAlias = Literal["atom", "residue", "chain"]
SapLevel: TypeAlias = Literal["atom", "residue"]
ProtonationMode: TypeAlias = Literal["all-charged", "heuristic", "explicit-only"]
AtomSubset: TypeAlias = Literal["ca", "backbone", "heavy", "all"]
ClusteringMethod: TypeAlias = Literal["k-medoids"]
SequenceList: TypeAlias = list[tuple[str, str]]
DsasaComponents: TypeAlias = tuple[float, float, float, float]

EXPORTED_FUNCTIONS: Final[tuple[str, ...]]
SASA_LEVELS: Final[tuple[SasaLevel, ...]]
SAP_LEVELS: Final[tuple[SapLevel, ...]]
ATOM_SUBSETS: Final[tuple[AtomSubset, ...]]
CLUSTERING_METHODS: Final[tuple[ClusteringMethod, ...]]
DEFAULTS: Final[dict[str, dict[str, object]]]
CONTACT_COLUMNS: Final[tuple[str, ...]]
SASA_COLUMNS: Final[dict[SasaLevel, tuple[str, ...]]]
RELATIVE_SASA_COLUMNS: Final[tuple[str, ...]]
SAP_COLUMNS: Final[dict[SapLevel, tuple[str, ...]]]
SEQUENCE_RETURN: Final[str]
PAIRWISE_RMSD_COLUMNS: Final[tuple[str, str, str]]
CLUSTER_COLUMNS: Final[tuple[str, str, str, str]]
