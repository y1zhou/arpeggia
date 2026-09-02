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

SASA_LEVELS = ("atom", "residue", "chain")
SAP_LEVELS = ("atom", "residue")
ATOM_SUBSETS = ("ca", "backbone", "heavy", "all")
CLUSTERING_METHODS = ("k-medoids",)

DEFAULTS = {
    "rmsd": {
        "model_num": 0,
        "superpose_residues": "",
        "rmsd_residues": "",
        "atoms": "ca",
    },
    "pairwise_rmsd": {
        "id_col": "id",
        "path_col": "path",
        "model_num": 0,
        "superpose_residues": "",
        "rmsd_residues": "",
        "atoms": "ca",
        "num_threads": 0,
        "bypass_mem_check": False,
    },
    "cluster_structs": {
        "method": "k-medoids",
        "max_iterations": 100,
        "model_num": 0,
        "superpose_residues": "",
        "rmsd_residues": "",
        "atoms": "ca",
        "num_threads": 0,
        "bypass_mem_check": False,
    },
    "contacts": {
        "groups": "/",
        "vdw_comp": 0.1,
        "dist_cutoff": 6.5,
        "ignore_zero_occupancy": False,
        "protonation": "all-charged",
        "ph": 7.4,
        "num_threads": 1,
    },
    "sasa": {
        "level": "atom",
        "probe_radius": 1.4,
        "n_points": 100,
        "model_num": 0,
        "chains": "",
        "num_threads": 1,
    },
    "dsasa": {
        "probe_radius": 1.4,
        "n_points": 100,
        "model_num": 0,
        "num_threads": 1,
    },
    "relative_sasa": {
        "probe_radius": 1.4,
        "n_points": 100,
        "model_num": 0,
        "chains": "",
        "num_threads": 1,
    },
    "sap_score": {
        "level": "residue",
        "probe_radius": 1.1,
        "n_points": 100,
        "model_num": 0,
        "sap_radius": 5.0,
        "chains": "",
        "num_threads": 1,
    },
    "sc": {
        "model_num": 0,
        "num_threads": 0,
    },
}

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

RELATIVE_SASA_COLUMNS = (
    "chain",
    "resn",
    "resi",
    "insertion",
    "sasa",
    "polar_sasa",
    "hydrophobic_sasa",
    "unclassified_sasa",
    "relative_sasa",
)

SAP_COLUMNS = {
    "atom": (
        "chain",
        "resn",
        "resi",
        "insertion",
        "atomn",
        "atomi",
        "sasa",
        "sap_score",
    ),
    "residue": (
        "chain",
        "resn",
        "resi",
        "insertion",
        "sc_sasa",
        "sap_score",
        "max_sc_asa",
        "relative_sc_sasa",
    ),
}

SEQUENCE_RETURN = "list[tuple[str, str]]"

PAIRWISE_RMSD_COLUMNS = ("id_1", "id_2", "rmsd")
CLUSTER_COLUMNS = ("id", "cluster_id", "medoid_id", "rmsd_to_medoid")
