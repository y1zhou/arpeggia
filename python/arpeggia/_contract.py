"""Shared Python exposure contract for arpeggia."""

EXPORTED_FUNCTIONS = (
    "rmsd",
    "pairwise_rmsd",
    "cluster_structs",
    "contacts",
    "sasa",
    "relative_sasa",
    "sap_score",
    "dsasa",
    "dsasa_components",
    "sc",
    "seq",
    "seqres",
)

CONTACT_COLUMNS = (
    "model",
    "interaction",
    "distance",
    "from_chain",
    "from_resn",
    "from_resi",
    "from_insertion",
    "from_altloc",
    "from_atomn",
    "from_atomi",
    "to_chain",
    "to_resn",
    "to_resi",
    "to_insertion",
    "to_altloc",
    "to_atomn",
    "to_atomi",
    "sc_centroid_dist",
    "sc_dihedral",
    "sc_centroid_angle",
)

SASA_COLUMNS = {
    "atom": (
        "atomi",
        "sasa",
        "chain",
        "resn",
        "resi",
        "insertion",
        "altloc",
        "atomn",
        "polarity",
    ),
    "residue": (
        "chain",
        "resn",
        "resi",
        "insertion",
        "sasa",
        "polar_sasa",
        "hydrophobic_sasa",
        "unclassified_sasa",
    ),
    "chain": (
        "chain",
        "sasa",
        "polar_sasa",
        "hydrophobic_sasa",
        "unclassified_sasa",
    ),
}

PAIRWISE_RMSD_COLUMNS = ("id_1", "id_2", "rmsd")
CLUSTER_COLUMNS = ("id", "cluster_id", "medoid_id", "rmsd_to_medoid")
