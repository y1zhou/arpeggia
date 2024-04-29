use core::fmt;

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum Interaction {
    StericClash,
    Covalent,
    VanDerWaals,

    // Electrostatic,
    Ionic,
    HydrogenBond,
    WeakHydrogenBond, // C-H...O hydrogen bond
    // HalogenBond
    Polar, // hydrogen bonding without angle terms

    // Aromatic
    PiDisplaced, // staggered stacking, parallel displaced
    PiT,         // perpendicular T-shaped
    PiSandwich,  // direct stacking, repulsive
    CationPi,

    Hydrophobic,
}

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct InteractingEntity {
    pub chain: String,
    pub resn: String,
    pub resi: isize,
    pub altloc: String,
    pub atomn: String,
    pub atomi: usize,
}
#[derive(Debug, Clone)]
pub struct ResultEntry {
    pub interaction: Interaction,
    pub ligand: InteractingEntity,
    pub receptor: InteractingEntity,
}

impl fmt::Display for InteractingEntity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Chain {chain}, Residue {resn} {resi}{altloc}, Atom {atom_name} {atom_idx}",
            chain = self.chain,
            resn = self.resn,
            resi = self.resi,
            altloc = self.altloc,
            atom_name = self.atomn,
            atom_idx = self.atomi
        )
    }
}

impl fmt::Display for ResultEntry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Ligand [{ligand}] has {intxn:?} with Receptor [{receptor}]",
            ligand = self.ligand,
            intxn = self.interaction,
            receptor = self.receptor
        )
    }
}
