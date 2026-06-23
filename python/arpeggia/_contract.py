"""Shared Python exposure contract for arpeggia."""

EXPORTED_FUNCTIONS = (
    "contacts",
    "sasa",
    "relative_sasa",
    "sap_score",
    "dsasa",
    "sc",
    "seq",
)

SASA_LEVELS = ("atom", "residue", "chain")
SAP_LEVELS = ("atom", "residue")

DEFAULTS = {
    "contacts": {
        "groups": "/",
        "vdw_comp": 0.1,
        "dist_cutoff": 6.5,
        "ignore_zero_occupancy": False,
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
        "probe_radius": 1.4,
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
    ),
    "residue": ("chain", "resn", "resi", "insertion", "sasa", "is_polar"),
    "chain": ("chain", "sasa"),
}

RELATIVE_SASA_COLUMNS = (
    "chain",
    "resn",
    "resi",
    "insertion",
    "sasa",
    "is_polar",
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
