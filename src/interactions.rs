use pdbtbx::*;
use std::collections::HashSet;

pub struct InteractionComplex {
    pub model: PDB,
    pub ligand: HashSet<String>,
    pub receptor: HashSet<String>,
    pub vdw_comp_factor: f64,
    pub interacting_threshold: f64,
}

impl InteractionComplex {
    pub fn new(model: PDB, groups: &str, vdw_comp_factor: f64, interacting_threshold: f64) -> Self {
        // Parse the first two fields in groups
        let sel_vec: Vec<&str> = groups.split("/").collect();
        let ligand_chains = sel_vec.get(0).unwrap();
        let receptor_chains = sel_vec.get(1).unwrap();

        // Create a HashSet of chains for the ligand and receptor
        let mut ligand: HashSet<String> = ligand_chains
            .split(",")
            .map(|c| c.to_string())
            .filter(|c| !c.is_empty())
            .collect();
        let mut receptor: HashSet<String> = receptor_chains
            .split(",")
            .map(|c| c.to_string())
            .filter(|c| !c.is_empty())
            .collect();

        // If there are no ligand or receptor chains, use all remaining chains
        if ligand.is_empty() ^ receptor.is_empty() {
            if ligand.is_empty() {
                ligand = model
                    .chains()
                    .map(|c| c.id().to_string())
                    .filter(|c| !receptor.contains(c))
                    .collect();
            } else {
                receptor = model
                    .chains()
                    .map(|c| c.id().to_string())
                    .filter(|c| !ligand.contains(c))
                    .collect();
            }
        }

        println!("Parsed ligands {ligand:?}; receptors {receptor:?}");
        Self {
            model,
            ligand,
            receptor,
            vdw_comp_factor,
            interacting_threshold,
        }
    }
}

pub trait Interactions {
    fn initialize(&self);
    fn run_arpeggio(&self);
    fn get_contacts(&self);
}

impl Interactions for InteractionComplex {
    fn initialize(&self) {
        todo!()
    }

    fn run_arpeggio(&self) {
        todo!()
    }

    fn get_contacts(&self) {
        todo!()
    }
}
