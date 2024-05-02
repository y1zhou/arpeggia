use core::fmt;

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum Interaction {
    StericClash,
    CovalentBond,
    VanDerWaalsContact,

    // Electrostatic,
    IonicBond,
    HydrogenBond,
    WeakHydrogenBond, // C-H...O hydrogen bond
    // HalogenBond
    PolarContact,     // hydrogen bonding without angle terms
    WeakPolarContact, // C-H...O bonding without angle terms

    // Aromatic
    PiDisplacedStacking, // staggered stacking, parallel displaced
    PiTStacking,         // perpendicular T-shaped
    PiSandwichStacking,  // direct stacking, repulsive
    CationPi,

    HydrophobicContact,
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
