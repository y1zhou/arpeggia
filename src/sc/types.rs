//! SC calculation types.

use super::vector3::{ScValue, Vec3};

/// Atom attention/visibility state.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
#[allow(clippy::enum_variant_names)]
pub enum Attention {
    /// Far from the interface; not considered for surface emission
    Far,
    /// Intermediate state for geometric constructions (from sc-rs algorithm).
    /// Currently unused as our attention assignment only uses Far/Buried states.
    /// The checks against this variant exist for algorithmic completeness.
    #[allow(dead_code)]
    Consider,
    /// Buried and flagged for interface processing
    #[default]
    Buried,
}

#[derive(Clone, Debug, Default)]
pub struct Atom {
    pub natom: i32,
    pub molecule: usize,
    pub radius: ScValue,
    /// Per-atom sampling density (~15 dots/Å²)
    pub density: ScValue,
    pub attention: Attention,
    /// Is atom accessible to solvent/contact surface
    pub accessible: bool,
    pub atom: String,
    pub residue: String,
    pub coor: Vec3,
    /// Neighbor indices on same molecule
    pub neighbor_indices: Vec<usize>,
    /// Neighbor indices on opposite molecule that bury this atom
    pub buried_by_indices: Vec<usize>,
}

impl Atom {
    pub fn distance_squared(&self, other: &Atom) -> ScValue {
        self.coor.distance_squared(other.coor)
    }
}

#[derive(Clone, Debug)]
pub struct Probe {
    /// Indices of the three atoms defining the probe center
    pub atom_indices: [usize; 3],
    pub height: ScValue,
    pub point: Vec3,
    pub alt: Vec3,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum DotKind {
    Contact,
    Reentrant,
    Cavity,
}

#[derive(Clone, Debug)]
pub struct Dot {
    /// Discretized surface point
    pub coor: Vec3,
    /// Outward unit normal at the point
    pub outnml: Vec3,
    pub area: ScValue,
    pub buried: bool,
    // The following fields are set during surface generation but currently
    // only used internally. They're kept for potential debugging and future use.
    #[allow(dead_code)]
    /// Type of surface (contact, reentrant, cavity)
    pub kind: DotKind,
    #[allow(dead_code)]
    /// Index of the atom this dot belongs to
    pub atom_index: usize,
}

#[derive(Clone, Debug, Default)]
pub struct DotStats {
    pub convex: usize,
    pub toroidal: usize,
    pub concave: usize,
}

#[derive(Clone, Debug, Default)]
pub struct SurfaceStats {
    pub n_atoms: usize,
    pub n_buried_atoms: usize,
    pub n_blocked_atoms: usize,
    pub d_mean: ScValue,
    pub d_median: ScValue,
    pub s_mean: ScValue,
    pub s_median: ScValue,
    pub n_all_dots: usize,
    pub n_trimmed_dots: usize,
    pub trimmed_area: ScValue,
}

#[derive(Clone, Debug, Default)]
pub struct Results {
    pub valid: i32,
    pub n_atoms: usize,
    pub surfaces: [SurfaceStats; 2],
    pub combined: SurfaceStats,
    pub dots: DotStats,
    pub sc: ScValue,
    pub distance: ScValue,
    pub area: ScValue,
}

#[derive(Clone, Debug, Default)]
pub struct AtomRadius {
    pub residue: String,
    pub atom: String,
    pub radius: ScValue,
}
