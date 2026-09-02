//! Sequence extraction from protein structures.
//!
//! This module provides functions for extracting amino acid sequences
//! from PDB structures.

use crate::contacts::one_letter_code;
use pdbtbx::*;

/// Get observed chain sequences from one model in a PDB structure.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `model_num` - Model number to select (`0` selects the first model)
///
/// # Returns
///
/// A `Vec` of `(chain_id, sequence)` pairs.
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_sequences};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let pdb = load_model(&input_file).unwrap().value;
/// let sequences = get_sequences(&pdb, 0).unwrap();
/// for (chain_id, seq) in sequences {
///     println!("Chain {}: {}", chain_id, seq);
/// }
/// ```
pub fn get_sequences(pdb: &PDB, model_num: usize) -> crate::ArpeggiaResult<Vec<(String, String)>> {
    let model = crate::structure::selected_model(pdb, model_num)?;
    Ok(model
        .chains()
        .map(|chain| {
            let sequence = chain
                .residues()
                .filter_map(|residue| residue.name().and_then(one_letter_code))
                .collect();
            (chain.id().to_string(), sequence)
        })
        .filter(|(_, sequence): &(String, String)| !sequence.is_empty())
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::load_model;

    #[test]
    fn observed_sequence_is_polymer_coordinates_only() {
        let path = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
        let pdb = load_model(&path).unwrap().value;
        assert_eq!(
            get_sequences(&pdb, 0).unwrap(),
            [(
                "A".into(),
                "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
                    .into()
            )]
        );
    }

    #[test]
    fn observed_sequence_keeps_supported_modified_monomers_but_not_water() {
        let input =
            b"HETATM    1  CA  MSE A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
HETATM    2  O   HOH A   2       1.000   0.000   0.000  1.00 20.00           O  \n\
ATOM      3  CA  GLY A   3       2.000   0.000   0.000  1.00 20.00           C  \n\
END                                                                             \n";
        let pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(std::io::BufReader::new(input.as_slice()))
            .unwrap()
            .0;
        assert_eq!(get_sequences(&pdb, 0).unwrap(), [("A".into(), "MG".into())]);
    }

    #[test]
    fn observed_sequence_omits_nonpolymer_only_chains() {
        let input =
            b"ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
HETATM    2  O   HOH W   1       2.000   0.000   0.000  1.00 20.00           O  \n\
END                                                                             \n";
        let pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(std::io::BufReader::new(input.as_slice()))
            .unwrap()
            .0;

        assert_eq!(get_sequences(&pdb, 0).unwrap(), [("A".into(), "A".into())]);
    }

    #[test]
    fn observed_sequence_selects_one_model() {
        let input = b"MODEL        7\n\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ENDMDL\n\
MODEL        9\n\
ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ENDMDL\nEND\n";
        let pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(std::io::BufReader::new(input.as_slice()))
            .unwrap()
            .0;

        assert_eq!(get_sequences(&pdb, 0).unwrap(), [("A".into(), "A".into())]);
        assert_eq!(get_sequences(&pdb, 9).unwrap(), [("A".into(), "G".into())]);
        assert!(matches!(
            get_sequences(&pdb, 3),
            Err(crate::ArpeggiaError::InvalidArgument(_))
        ));
    }
}
