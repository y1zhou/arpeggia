//! Atomic radii lookup for SC calculations.
//!
//! Uses embedded radii table from sc-rs for consistency.

use super::types::AtomRadius;
use super::vector3::ScValue;

/// Embedded atomic radii table (from Lawrence & Colman).
const ATOMIC_RADII_JSON: &str = include_str!("atomic_radii.json");

pub fn embedded_atomic_radii() -> Vec<AtomRadius> {
    read_atomic_radii_from_str(ATOMIC_RADII_JSON).unwrap_or_default()
}

#[derive(serde::Deserialize)]
struct RadiusRecord {
    residue: String,
    atom: String,
    radius: ScValue,
}

pub fn read_atomic_radii_from_str(data: &str) -> std::io::Result<Vec<AtomRadius>> {
    let recs: Vec<RadiusRecord> = serde_json::from_str(data)
        .map_err(|e| std::io::Error::other(format!("invalid radii json: {e}")))?;
    Ok(recs
        .into_iter()
        .filter(|r| r.radius > 0.0)
        .map(|r| AtomRadius {
            residue: r.residue,
            atom: r.atom,
            radius: r.radius,
        })
        .collect())
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
