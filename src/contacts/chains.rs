use super::residues::ResidueExt;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::load_model;

    #[test]
    fn test_pdb_seq() {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/1ubq.pdb");

        let (pdb, _) = load_model(&path);
        let chain = pdb.model(0).unwrap().chain(0).unwrap();
        let seq = chain.pdb_seq().join("");
        let ubiquitin =
            "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
        assert!(
            seq.starts_with(ubiquitin),
            "Expected: {ubiquitin}\nFound: {seq}"
        );
        let all_waters = seq.matches('O').count();
        assert_eq!(
            all_waters, 58,
            "Found {all_waters} waters in the model instead of 58"
        );
    }
}
