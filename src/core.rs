use crate::interactions::{InteractionComplex, Interactions, ResultEntry};
use crate::utils::load_model;

use pdbtbx::*;
use polars::prelude::*;
use std::path::{Path, PathBuf};

/// Core function of the library. Argument documentation can be found in [`main::Args`].
pub fn core(
    input_file: PathBuf,
    groups: &str,
    vdw_comp: f64,
    dist_cutoff: f64,
) -> Result<(DataFrame, DataFrame, InteractionComplex, Vec<PDBError>), Vec<PDBError>> {
    // Make sure `input` exists
    let input_path = Path::new(&input_file).canonicalize().unwrap();
    let input_file: String = input_path.to_str().unwrap().parse().unwrap();

    // Load file as complex structure
    let (pdb, errors) = load_model(&input_file)?;

    let i_complex = InteractionComplex::new(pdb, groups, vdw_comp, dist_cutoff);

    // Find interactions
    let atomic_contacts = i_complex.get_atomic_contacts();
    let df_atomic = results_to_df(&atomic_contacts);

    let mut ring_contacts: Vec<ResultEntry> = Vec::new();
    ring_contacts.extend(i_complex.get_ring_atom_contacts());
    ring_contacts.extend(i_complex.get_ring_ring_contacts());
    let df_ring = results_to_df(&ring_contacts)
        // .drop_many(&["from_atomn", "from_atomi", "to_atomn", "to_atomi"])
        .sort(
            [
                "from_chain",
                "from_resi",
                "from_altloc",
                "to_chain",
                "to_resi",
                "to_altloc",
            ],
            Default::default(),
        )
        .unwrap();

    Ok((df_atomic, df_ring, i_complex, errors))
}

fn results_to_df(res: &[ResultEntry]) -> DataFrame {
    df!(
        "interaction" => res.iter().map(|x| x.interaction.to_string()).collect::<Vec<String>>(),
        "distance" => res.iter().map(|x| x.distance).collect::<Vec<f64>>(),
        "from_chain" => res.iter().map(|x| x.ligand.chain.to_owned()).collect::<Vec<String>>(),
        "from_resn" => res.iter().map(|x| x.ligand.resn.to_owned()).collect::<Vec<String>>(),
        "from_resi" => res.iter().map(|x| x.ligand.resi as i64).collect::<Vec<i64>>(),
        "from_altloc" => res.iter().map(|x| x.ligand.altloc.to_owned()).collect::<Vec<String>>(),
        "from_atomn" => res.iter().map(|x| x.ligand.atomn.to_owned()).collect::<Vec<String>>(),
        "from_atomi" => res.iter().map(|x| x.ligand.atomi as i64).collect::<Vec<i64>>(),
        "to_chain" => res.iter().map(|x| x.receptor.chain.to_owned()).collect::<Vec<String>>(),
        "to_resn" => res.iter().map(|x| x.receptor.resn.to_owned()).collect::<Vec<String>>(),
        "to_resi" => res.iter().map(|x| x.receptor.resi as i64).collect::<Vec<i64>>(),
        "to_altloc" => res.iter().map(|x| x.receptor.altloc.to_owned()).collect::<Vec<String>>(),
        "to_atomn" => res.iter().map(|x| x.receptor.atomn.to_owned()).collect::<Vec<String>>(),
        "to_atomi" => res.iter().map(|x| x.receptor.atomi as i64).collect::<Vec<i64>>(),
    )
    .unwrap()
}
