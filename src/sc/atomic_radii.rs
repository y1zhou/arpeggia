//! Atomic radii lookup for SC calculations.
//!
//! Uses embedded radii table from sc-rs for consistency (Lawrence & Colman).

use super::types::AtomRadius;

/// Embedded atomic radii table (from Lawrence & Colman).
/// Radii are in Angstroms.
pub fn embedded_atomic_radii() -> Vec<AtomRadius> {
    vec![
        // Specific residue radii
        AtomRadius {
            residue: "ALA".to_string(),
            atom: "CB".to_string(),
            radius: 1.95,
        },
        AtomRadius {
            residue: "ARG".to_string(),
            atom: "NH*".to_string(),
            radius: 1.70,
        },
        AtomRadius {
            residue: "ARG".to_string(),
            atom: "CZ".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "ARG".to_string(),
            atom: "NE".to_string(),
            radius: 1.65,
        },
        AtomRadius {
            residue: "ARG".to_string(),
            atom: "CD".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "ARG".to_string(),
            atom: "CG".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "ASN".to_string(),
            atom: "ND2".to_string(),
            radius: 1.70,
        },
        AtomRadius {
            residue: "ASN".to_string(),
            atom: "OD1".to_string(),
            radius: 1.60,
        },
        AtomRadius {
            residue: "ASN".to_string(),
            atom: "CG".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "ASP".to_string(),
            atom: "OD*".to_string(),
            radius: 1.60,
        },
        AtomRadius {
            residue: "ASP".to_string(),
            atom: "CG".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "GLN".to_string(),
            atom: "NE2".to_string(),
            radius: 1.70,
        },
        AtomRadius {
            residue: "GLN".to_string(),
            atom: "OE1".to_string(),
            radius: 1.60,
        },
        AtomRadius {
            residue: "GLN".to_string(),
            atom: "CD".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "GLN".to_string(),
            atom: "CG".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "GLU".to_string(),
            atom: "OE*".to_string(),
            radius: 1.60,
        },
        AtomRadius {
            residue: "GLU".to_string(),
            atom: "CD".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "GLU".to_string(),
            atom: "CG".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "GLY".to_string(),
            atom: "CA".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "HIS".to_string(),
            atom: "CD2".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "HIS".to_string(),
            atom: "NE2".to_string(),
            radius: 1.65,
        },
        AtomRadius {
            residue: "HIS".to_string(),
            atom: "CE1".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "HIS".to_string(),
            atom: "ND1".to_string(),
            radius: 1.65,
        },
        AtomRadius {
            residue: "HIS".to_string(),
            atom: "CG".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "HOH".to_string(),
            atom: "O**".to_string(),
            radius: 1.70,
        },
        AtomRadius {
            residue: "ILE".to_string(),
            atom: "CD1".to_string(),
            radius: 1.95,
        },
        AtomRadius {
            residue: "ILE".to_string(),
            atom: "CG1".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "ILE".to_string(),
            atom: "CB".to_string(),
            radius: 1.85,
        },
        AtomRadius {
            residue: "ILE".to_string(),
            atom: "CG2".to_string(),
            radius: 1.95,
        },
        AtomRadius {
            residue: "LEU".to_string(),
            atom: "CD*".to_string(),
            radius: 1.95,
        },
        AtomRadius {
            residue: "LEU".to_string(),
            atom: "CG".to_string(),
            radius: 1.85,
        },
        AtomRadius {
            residue: "LYS".to_string(),
            atom: "NZ".to_string(),
            radius: 1.75,
        },
        AtomRadius {
            residue: "LYS".to_string(),
            atom: "CE".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "LYS".to_string(),
            atom: "CD".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "LYS".to_string(),
            atom: "CG".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "MET".to_string(),
            atom: "CE".to_string(),
            radius: 1.95,
        },
        AtomRadius {
            residue: "MET".to_string(),
            atom: "CG".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "PHE".to_string(),
            atom: "CD*".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "PHE".to_string(),
            atom: "CE*".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "PHE".to_string(),
            atom: "CZ".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "PHE".to_string(),
            atom: "CG".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "PRO".to_string(),
            atom: "CD".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "PRO".to_string(),
            atom: "CG".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "SER".to_string(),
            atom: "OG".to_string(),
            radius: 1.70,
        },
        AtomRadius {
            residue: "SUL".to_string(),
            atom: "S".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "SUL".to_string(),
            atom: "O***".to_string(),
            radius: 1.65,
        },
        AtomRadius {
            residue: "THR".to_string(),
            atom: "CG2".to_string(),
            radius: 1.95,
        },
        AtomRadius {
            residue: "THR".to_string(),
            atom: "OG1".to_string(),
            radius: 1.70,
        },
        AtomRadius {
            residue: "THR".to_string(),
            atom: "CB".to_string(),
            radius: 1.85,
        },
        AtomRadius {
            residue: "TRP".to_string(),
            atom: "CE2".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "TRP".to_string(),
            atom: "CE3".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "TRP".to_string(),
            atom: "CD1".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "TRP".to_string(),
            atom: "CD2".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "TRP".to_string(),
            atom: "CZ*".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "TRP".to_string(),
            atom: "CH2".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "TRP".to_string(),
            atom: "NE1".to_string(),
            radius: 1.65,
        },
        AtomRadius {
            residue: "TRP".to_string(),
            atom: "CG".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "TYR".to_string(),
            atom: "OH".to_string(),
            radius: 1.70,
        },
        AtomRadius {
            residue: "TYR".to_string(),
            atom: "CD*".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "TYR".to_string(),
            atom: "CE*".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "TYR".to_string(),
            atom: "CZ".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "TYR".to_string(),
            atom: "CG".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "VAL".to_string(),
            atom: "CG*".to_string(),
            radius: 1.95,
        },
        AtomRadius {
            residue: "VAL".to_string(),
            atom: "CB".to_string(),
            radius: 1.85,
        },
        AtomRadius {
            residue: "WAT".to_string(),
            atom: "O".to_string(),
            radius: 1.70,
        },
        AtomRadius {
            residue: "WAT".to_string(),
            atom: "O*".to_string(),
            radius: 1.70,
        },
        // Generic patterns (matching any residue)
        AtomRadius {
            residue: "***".to_string(),
            atom: "H".to_string(),
            radius: 0.50,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "H*".to_string(),
            radius: 0.50,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "H**".to_string(),
            radius: 0.50,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "H***".to_string(),
            radius: 0.50,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "CA".to_string(),
            radius: 1.85,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "C".to_string(),
            radius: 1.80,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "O".to_string(),
            radius: 1.60,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "N".to_string(),
            radius: 1.65,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "CB".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "OT*".to_string(),
            radius: 1.60,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "OXT".to_string(),
            radius: 1.60,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "S*".to_string(),
            radius: 1.90,
        },
        AtomRadius {
            residue: "***".to_string(),
            atom: "P".to_string(),
            radius: 1.80,
        },
    ]
}

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
