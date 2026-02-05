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
    if is_hydrophobic(entity1.residue().name().unwrap(), e1_atom.name())
        && is_hydrophobic(entity2.residue().name().unwrap(), e2_atom.name())
        && (e1_atom.distance(e2_atom) <= HYDROPHOBIC_CONTACT_DIST)
    {
        Some(Interaction::HydrophobicContact)
    } else {
        None
    }
}

/// Check if the entity has an atom that belongs to the hydrophobic category.
fn is_hydrophobic(res_name: &str, atom_name: &str) -> bool {
    // Carbon beta of all other than Glycine/Serine
    if (atom_name == "CB") && (res_name != "SER") {
        return true;
    }
    matches!(
        (res_name, atom_name),
        ("ARG" | "GLN" | "GLU" | "PRO", "CG")
            | ("ILE", "CG1" | "CD1" | "CG2")
            | ("LEU", "CG" | "CD1" | "CD2")
            | ("LYS", "CG" | "CD")
            | ("MET", "CG" | "CE" | "SD") // sulfur in CYS has a hydrogen and is polarized
            | ("PHE", "CG" | "CD1" | "CD2" | "CE1" | "CE2" | "CZ")
            | ("THR", "CG2")
            | ("TRP", "CG" | "CD2" | "CE3" | "CZ3" | "CH2" | "CZ2")
            | ("TYR", "CG" | "CD1" | "CD2" | "CE1" | "CE2" )
            | ("VAL", "CG1" | "CG2")
    )
}
