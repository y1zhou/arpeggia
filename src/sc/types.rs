//! SC calculation types.

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
    pub attention: Attention,
    /// Is atom accessible to solvent/contact surface
    pub accessible: bool,
    pub atomn: String,
    pub resn: String,
    pub coor: Vec3,
    /// Same-molecule neighbors sorted by distance.
    pub neighbor_indices: Vec<usize>,
    /// Same-molecule neighbor distances sorted by compact atom index.
    pub neighbor_distances: Vec<(usize, f64)>,
}

impl ScAtom {
    pub fn neighbor_distance_squared(&self, index: usize) -> Option<f64> {
        self.neighbor_distances
            .binary_search_by_key(&index, |&(neighbor, _)| neighbor)
            .ok()
            .map(|position| self.neighbor_distances[position].1)
    }
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
