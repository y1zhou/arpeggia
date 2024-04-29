use crate::residues::ResidueExt;
use pdbtbx::*;
use rayon::prelude::*;

pub trait ChainExt {
    fn pdb_seq(&self) -> Vec<&str>;
}

impl ChainExt for Chain {
    fn pdb_seq(&self) -> Vec<&str> {
        // Load the amino acid sequence for each chain
        self.par_residues().map(|res| res.resn().unwrap()).collect()
    }
}
