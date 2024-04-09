use pdbtbx::*;

pub struct InteractionComplex {
    pub model: PDB,
    pub vdw_comp_factor: f64,
    pub interacting_threshold: f64,
}

pub trait Interactions {
    fn add_hydrogens(&self);
    fn initialize(&self);
    fn run_arpeggio(&self);
    fn get_contacts(&self);
}

impl Interactions for InteractionComplex {
    fn add_hydrogens(&self) {
        todo!()
    }

    fn get_contacts(&self) {
        todo!()
    }

    fn initialize(&self) {
        todo!()
    }

    fn run_arpeggio(&self) {
        todo!()
    }
}
