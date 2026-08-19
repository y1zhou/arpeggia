//! SC calculation types.

use std::collections::HashMap;

use super::vector3::Vec3;

/// Atom attention/visibility state.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum Attention {
    /// Far from the interface; not considered for surface emission
    Far,
    /// Buried and flagged for interface processing
    #[default]
    Buried,
}

#[derive(Clone, Debug, Default)]
pub struct ScAtom {
    pub atomi: usize,
    pub molecule: usize,
    pub radius: f64,
    /// Per-atom sampling density (~15 dots/Å²)
    pub density: f64,
    pub attention: Attention,
    /// Is atom accessible to solvent/contact surface
    pub accessible: bool,
    pub atomn: String,
    pub resn: String,
    pub coor: Vec3,
    /// Indices of all neighbors and their distance^2 to limit search space
    pub neighbors_atomi_dist2: HashMap<usize, f64>,
    /// Neighbor indices on same molecule
    pub neighbor_indices: Vec<usize>,
    /// Neighbor indices on opposite molecule that bury this atom
    pub buried_by_indices: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct Probe {
    /// Indices of the three atoms defining the probe center
    pub atom_indices: [usize; 3],
    pub height: f64,
    pub point: Vec3,
    pub alt: Vec3,
}

#[derive(Clone, Debug)]
pub struct Dot {
    /// Discretized surface point
    pub coor: Vec3,
    /// Outward unit normal at the point
    pub outnml: Vec3,
    pub buried: bool,
}

#[derive(Clone, Debug, Default)]
pub struct SurfaceStats {
    pub s_median: f64,
}

#[derive(Clone, Debug, Default)]
pub struct Results {
    pub surfaces: [SurfaceStats; 2],
    pub sc: f64,
}

#[derive(Clone, Debug, Default)]
pub struct AtomRadius {
    pub residue: &'static str,
    pub atom: &'static str,
    pub radius: f64,
}
