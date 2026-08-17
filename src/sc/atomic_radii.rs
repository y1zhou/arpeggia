//! Atomic radii lookup for SC calculations.
//!
//! Uses embedded radii table from sc-rs for consistency (Lawrence & Colman).

use super::types::AtomRadius;

/// Embedded atomic radii table (from Lawrence & Colman).
/// Radii are in Angstroms.
pub static EMBEDDED_ATOMIC_RADII: &[AtomRadius] = &[
    // Specific residue radii
    AtomRadius {
        residue: "ALA",
        atom: "CB",
        radius: 1.95,
    },
    AtomRadius {
        residue: "ARG",
        atom: "NH*",
        radius: 1.70,
    },
    AtomRadius {
        residue: "ARG",
        atom: "CZ",
        radius: 1.80,
    },
    AtomRadius {
        residue: "ARG",
        atom: "NE",
        radius: 1.65,
    },
    AtomRadius {
        residue: "ARG",
        atom: "CD",
        radius: 1.90,
    },
    AtomRadius {
        residue: "ARG",
        atom: "CG",
        radius: 1.90,
    },
    AtomRadius {
        residue: "ASN",
        atom: "ND2",
        radius: 1.70,
    },
    AtomRadius {
        residue: "ASN",
        atom: "OD1",
        radius: 1.60,
    },
    AtomRadius {
        residue: "ASN",
        atom: "CG",
        radius: 1.80,
    },
    AtomRadius {
        residue: "ASP",
        atom: "OD*",
        radius: 1.60,
    },
    AtomRadius {
        residue: "ASP",
        atom: "CG",
        radius: 1.80,
    },
    AtomRadius {
        residue: "GLN",
        atom: "NE2",
        radius: 1.70,
    },
    AtomRadius {
        residue: "GLN",
        atom: "OE1",
        radius: 1.60,
    },
    AtomRadius {
        residue: "GLN",
        atom: "CD",
        radius: 1.80,
    },
    AtomRadius {
        residue: "GLN",
        atom: "CG",
        radius: 1.90,
    },
    AtomRadius {
        residue: "GLU",
        atom: "OE*",
        radius: 1.60,
    },
    AtomRadius {
        residue: "GLU",
        atom: "CD",
        radius: 1.80,
    },
    AtomRadius {
        residue: "GLU",
        atom: "CG",
        radius: 1.90,
    },
    AtomRadius {
        residue: "GLY",
        atom: "CA",
        radius: 1.90,
    },
    AtomRadius {
        residue: "HIS",
        atom: "CD2",
        radius: 1.90,
    },
    AtomRadius {
        residue: "HIS",
        atom: "NE2",
        radius: 1.65,
    },
    AtomRadius {
        residue: "HIS",
        atom: "CE1",
        radius: 1.90,
    },
    AtomRadius {
        residue: "HIS",
        atom: "ND1",
        radius: 1.65,
    },
    AtomRadius {
        residue: "HIS",
        atom: "CG",
        radius: 1.80,
    },
    AtomRadius {
        residue: "HOH",
        atom: "O**",
        radius: 1.70,
    },
    AtomRadius {
        residue: "ILE",
        atom: "CD1",
        radius: 1.95,
    },
    AtomRadius {
        residue: "ILE",
        atom: "CG1",
        radius: 1.90,
    },
    AtomRadius {
        residue: "ILE",
        atom: "CB",
        radius: 1.85,
    },
    AtomRadius {
        residue: "ILE",
        atom: "CG2",
        radius: 1.95,
    },
    AtomRadius {
        residue: "LEU",
        atom: "CD*",
        radius: 1.95,
    },
    AtomRadius {
        residue: "LEU",
        atom: "CG",
        radius: 1.85,
    },
    AtomRadius {
        residue: "LYS",
        atom: "NZ",
        radius: 1.75,
    },
    AtomRadius {
        residue: "LYS",
        atom: "CE",
        radius: 1.90,
    },
    AtomRadius {
        residue: "LYS",
        atom: "CD",
        radius: 1.90,
    },
    AtomRadius {
        residue: "LYS",
        atom: "CG",
        radius: 1.90,
    },
    AtomRadius {
        residue: "MET",
        atom: "CE",
        radius: 1.95,
    },
    AtomRadius {
        residue: "MET",
        atom: "CG",
        radius: 1.90,
    },
    AtomRadius {
        residue: "PHE",
        atom: "CD*",
        radius: 1.90,
    },
    AtomRadius {
        residue: "PHE",
        atom: "CE*",
        radius: 1.90,
    },
    AtomRadius {
        residue: "PHE",
        atom: "CZ",
        radius: 1.90,
    },
    AtomRadius {
        residue: "PHE",
        atom: "CG",
        radius: 1.80,
    },
    AtomRadius {
        residue: "PRO",
        atom: "CD",
        radius: 1.90,
    },
    AtomRadius {
        residue: "PRO",
        atom: "CG",
        radius: 1.90,
    },
    AtomRadius {
        residue: "SER",
        atom: "OG",
        radius: 1.70,
    },
    AtomRadius {
        residue: "SUL",
        atom: "S",
        radius: 1.90,
    },
    AtomRadius {
        residue: "SUL",
        atom: "O***",
        radius: 1.65,
    },
    AtomRadius {
        residue: "THR",
        atom: "CG2",
        radius: 1.95,
    },
    AtomRadius {
        residue: "THR",
        atom: "OG1",
        radius: 1.70,
    },
    AtomRadius {
        residue: "THR",
        atom: "CB",
        radius: 1.85,
    },
    AtomRadius {
        residue: "TRP",
        atom: "CE2",
        radius: 1.80,
    },
    AtomRadius {
        residue: "TRP",
        atom: "CE3",
        radius: 1.90,
    },
    AtomRadius {
        residue: "TRP",
        atom: "CD1",
        radius: 1.90,
    },
    AtomRadius {
        residue: "TRP",
        atom: "CD2",
        radius: 1.80,
    },
    AtomRadius {
        residue: "TRP",
        atom: "CZ*",
        radius: 1.90,
    },
    AtomRadius {
        residue: "TRP",
        atom: "CH2",
        radius: 1.90,
    },
    AtomRadius {
        residue: "TRP",
        atom: "NE1",
        radius: 1.65,
    },
    AtomRadius {
        residue: "TRP",
        atom: "CG",
        radius: 1.80,
    },
    AtomRadius {
        residue: "TYR",
        atom: "OH",
        radius: 1.70,
    },
    AtomRadius {
        residue: "TYR",
        atom: "CD*",
        radius: 1.90,
    },
    AtomRadius {
        residue: "TYR",
        atom: "CE*",
        radius: 1.90,
    },
    AtomRadius {
        residue: "TYR",
        atom: "CZ",
        radius: 1.80,
    },
    AtomRadius {
        residue: "TYR",
        atom: "CG",
        radius: 1.80,
    },
    AtomRadius {
        residue: "VAL",
        atom: "CG*",
        radius: 1.95,
    },
    AtomRadius {
        residue: "VAL",
        atom: "CB",
        radius: 1.85,
    },
    AtomRadius {
        residue: "WAT",
        atom: "O",
        radius: 1.70,
    },
    AtomRadius {
        residue: "WAT",
        atom: "O*",
        radius: 1.70,
    },
    // Generic patterns (matching any residue)
    AtomRadius {
        residue: "***",
        atom: "H",
        radius: 0.50,
    },
    AtomRadius {
        residue: "***",
        atom: "H*",
        radius: 0.50,
    },
    AtomRadius {
        residue: "***",
        atom: "H**",
        radius: 0.50,
    },
    AtomRadius {
        residue: "***",
        atom: "H***",
        radius: 0.50,
    },
    AtomRadius {
        residue: "***",
        atom: "CA",
        radius: 1.85,
    },
    AtomRadius {
        residue: "***",
        atom: "C",
        radius: 1.80,
    },
    AtomRadius {
        residue: "***",
        atom: "O",
        radius: 1.60,
    },
    AtomRadius {
        residue: "***",
        atom: "N",
        radius: 1.65,
    },
    AtomRadius {
        residue: "***",
        atom: "CB",
        radius: 1.90,
    },
    AtomRadius {
        residue: "***",
        atom: "OT*",
        radius: 1.60,
    },
    AtomRadius {
        residue: "***",
        atom: "OXT",
        radius: 1.60,
    },
    AtomRadius {
        residue: "***",
        atom: "S*",
        radius: 1.90,
    },
    AtomRadius {
        residue: "***",
        atom: "P",
        radius: 1.80,
    },
];

/// Wildcard match for residue/atom patterns.
/// `*` at start matches everything, `*` elsewhere matches suffix.
pub fn wildcard_match(query: &str, pattern: &str) -> bool {
    fn rtrim_spaces(s: &str) -> &str {
        let mut end = s.len();
        let b = s.as_bytes();
        while end > 0 && (b[end - 1] as char) == ' ' {
            end -= 1;
        }
        &s[..end]
    }

    let q = rtrim_spaces(query);
    let p = rtrim_spaces(pattern);

    if p.starts_with('*') {
        return true;
    }

    if let Some(star) = p.find('*') {
        let plen = star;
        if q.len() < plen {
            return false;
        }
        return q[..plen] == p[..plen];
    }

    // No '*' in pattern: exact match
    q == p
}
