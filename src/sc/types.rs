//! SC calculation types.

use super::vector3::Vec3;

/// Atom attention/visibility state.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(super) enum Attention {
    /// Far from the interface; not considered for surface emission
    Far,
    /// Buried and flagged for interface processing
    #[default]
    Buried,
}

#[derive(Clone, Debug, Default)]
pub(super) struct ScAtom {
    pub(super) atomi: usize,
    pub(super) molecule: u8,
    pub(super) radius: f64,
    pub(super) attention: Attention,
    /// Is atom accessible to solvent/contact surface
    pub(super) accessible: bool,
    pub(super) atomn: String,
    pub(super) resn: String,
    pub(super) coor: Vec3,
    /// Same-molecule neighbors sorted by distance.
    pub(super) neighbor_indices: Vec<usize>,
    /// Same-molecule neighbor distances sorted by compact atom index.
    pub(super) neighbor_distances: Vec<(usize, f64)>,
}

impl ScAtom {
    pub(super) fn neighbor_distance_squared(&self, index: usize) -> Option<f64> {
        self.neighbor_distances
            .binary_search_by_key(&index, |&(neighbor, _)| neighbor)
            .ok()
            .map(|position| self.neighbor_distances[position].1)
    }
}

#[derive(Clone, Debug)]
pub(super) struct Probe {
    /// Indices of the three atoms defining the probe center
    pub(super) atom_indices: [usize; 3],
    pub(super) height: f64,
    pub(super) point: Vec3,
    pub(super) alt: Vec3,
}

#[derive(Clone, Debug)]
pub(super) struct Dot {
    /// Discretized surface point
    pub(super) coor: Vec3,
    /// Outward unit normal at the point
    pub(super) outnml: Vec3,
    pub(super) buried: bool,
}

#[derive(Clone, Debug, Default)]
pub(super) struct SurfaceStats {
    pub(super) s_median: f64,
}

#[derive(Clone, Debug, Default)]
pub(super) struct Results {
    pub(super) surfaces: [SurfaceStats; 2],
    pub(super) sc: f64,
}

#[derive(Clone, Debug, Default)]
pub(super) struct AtomRadius {
    pub(super) residue: &'static str,
    pub(super) atom: &'static str,
    pub(super) radius: f64,
}
