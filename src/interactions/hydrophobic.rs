use super::structs::Interaction;

use pdbtbx::*;

const HYDROPHOBIC_CONTACT_DIST: f64 = 4.5;

/// Search for hydrophobic contacts.
///
/// Check if the distance between two hydrophobics is within [`HYDROPHOBIC_CONTACT_DIST`].
pub fn find_hydrophobic_contact(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
) -> Option<Interaction> {
    let e1_atom = entity1.atom();
    let e2_atom = entity2.atom();
    match is_hydrophobic(entity1.residue().name().unwrap(), e1_atom.name())
        & is_hydrophobic(entity2.residue().name().unwrap(), e2_atom.name())
        & (e1_atom.distance(e2_atom) <= HYDROPHOBIC_CONTACT_DIST)
    {
        true => Some(Interaction::HydrophobicContact),
        false => None,
    }
}

/// Check if the entity has an atom that belongs to the hydrophobic category.
fn is_hydrophobic(res_name: &str, atom_name: &str) -> bool {
    // Carbon beta of all other than Glycine/Serine
    if (atom_name == "CB") & (res_name != "SER") {
        return true;
    }
    matches!(
        (res_name, atom_name),
        ("ARG", "CG") | ("GLN", "CG")
            | ("GLU", "CG")
            | ("ILE", "CG1")
            | ("ILE", "CD1")
            | ("ILE", "CG2")
            | ("LEU", "CG")
            | ("LEU", "CD1")
            | ("LEU", "CD2")
            | ("LYS", "CG")
            | ("LYS", "CD")
            | ("MET", "CG")
            | ("MET", "SD") // sulfur in CYS has a hydrogen and is polarized
            | ("MET", "CE")
            | ("PHE", "CG")
            | ("PHE", "CD1")
            | ("PHE", "CD2")
            | ("PHE", "CE1")
            | ("PHE", "CE2")
            | ("PHE", "CZ")
            | ("PRO", "CG")
            | ("THR", "CG2")
            | ("TRP", "CG")
            | ("TRP", "CD2")
            | ("TRP", "CE3")
            | ("TRP", "CZ3")
            | ("TRP", "CH2")
            | ("TRP", "CZ2")
            | ("TYR", "CG")
            | ("TYR", "CD1")
            | ("TYR", "CD2")
            | ("TYR", "CE1")
            | ("TYR", "CE2")
            | ("VAL", "CG1")
            | ("VAL", "CG2")
    )
}
