use crate::residues::ResidueExt;
use pdbtbx::*;
use polars::prelude::*;
use rust_sasa::calculate_sasa_internal;
use rust_sasa::Atom as SASAAtom;

use crate::interactions::InteractingEntity;

pub fn get_atom_sasa(pdb: &PDB) -> DataFrame {
    // Calculate the SASA for each atom
    let atoms = pdb
        .atoms_with_hierarchy()
        .filter(|x| {
            let resn = x.residue().resn().unwrap();
            resn != "O" && resn != "X"
        })
        .map(|x| SASAAtom {
            position: nalgebra::Point3::new(
                x.atom().pos().0 as f32,
                x.atom().pos().1 as f32,
                x.atom().pos().2 as f32,
            ),
            radius: x
                .atom()
                .element()
                .unwrap()
                .atomic_radius()
                .van_der_waals
                .unwrap() as f32,
            id: x.atom().serial_number(),
            parent_id: None,
        })
        .collect::<Vec<_>>();
    let atom_sasa = calculate_sasa_internal(&atoms, None, Some(960));

    // Create a DataFrame with the results
    let atom_annotations = pdb
        .atoms_with_hierarchy()
        .map(|x| {
            let (resi, insertion) = x.residue().id();
            let altloc = match insertion {
                Some(insertion) => insertion.to_string(),
                None => "".to_string(),
            };
            InteractingEntity {
                chain: x.chain().id().to_string(),
                resn: x.residue().name().unwrap().to_string(),
                resi,
                altloc,
                atomn: x.atom().name().to_string(),
                atomi: x.atom().serial_number(),
            }
        })
        .collect::<Vec<InteractingEntity>>();
    let atom_annot_df = df!(
        "chain" => atom_annotations.iter().map(|x| x.chain.to_owned()).collect::<Vec<String>>(),
        "resn" => atom_annotations.iter().map(|x| x.resn.to_owned()).collect::<Vec<String>>(),
        "resi" => atom_annotations.iter().map(|x| x.resi as i64).collect::<Vec<i64>>(),
        "altloc" => atom_annotations.iter().map(|x| x.altloc.to_owned()).collect::<Vec<String>>(),
        "atomn" => atom_annotations.iter().map(|x| x.atomn.to_owned()).collect::<Vec<String>>(),
        "atomi" => atom_annotations.iter().map(|x| x.atomi as i64).collect::<Vec<i64>>(),
    )
    .unwrap();

    df!(
        "atomi" => atoms.iter().map(|x| x.id as i64).collect::<Vec<i64>>(),
        "sasa" => atom_sasa
    )
    .unwrap()
    .join(
        &atom_annot_df,
        ["atomi"],
        ["atomi"],
        JoinArgs::new(JoinType::Inner),
    )
    .unwrap()
    .sort(["chain", "resi", "altloc", "atomi"], Default::default())
    .unwrap()
}
