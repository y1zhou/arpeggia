use core::fmt;

/// Interaction types.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Interaction {
    /// Within covalent radii
    StericClash,
    /// Covalent bonded
    CovalentBond,
    // CB-S-S-CB bond
    Disulfide,
    /// Within van der Waals radii
    VanDerWaalsContact,

    // Electrostatic
    /// Negative and positive charged ions
    IonicBond,
    /// Donor-H-acceptor hydrogen bonds
    HydrogenBond,
    /// C-H...O hydrogen bond
    WeakHydrogenBond,
    // HalogenBond,
    /// Hydrogen bonding without angle terms
    PolarContact,
    /// C-H...O bonding without angle terms
    WeakPolarContact,

    // Aromatic
    /// Pi-pi staggered stacking, parallel displaced (`of`)
    PiDisplacedStacking,
    /// Pi-pi perpendicular T-shaped stacking (`fe` and `ef`)
    PiTStacking,
    /// Pi-pi direct stacking, repulsive (`ff`)
    PiSandwichStacking,
    /// Pi-pi parallel in-plane interactions (`ee`)
    PiParallelInPlaneStacking,
    /// Pi-pi tilted interactions (`ft`, `ot`, and `et`)
    PiTiltedStacking,
    /// Pi-pi cogwheel or L-shaped interactions (`oe`)
    PiLStacking,
    /// Cation-pi interaction
    CationPi,

    /// Hydrophobic interaction
    HydrophobicContact,
}

/// Entity interacting with another.
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct InteractingEntity {
    /// Chain identifier
    pub chain: String,
    /// Residue name
    pub resn: String,
    /// Residue index
    pub resi: isize,
    /// Alternate location identifier
    pub altloc: String,
    /// Atom name
    pub atomn: String,
    /// Atom index
    pub atomi: usize,
}

/// Entry passed to the results.
#[derive(Debug, Clone)]
pub struct ResultEntry {
    /// Interaction type
    pub interaction: Interaction,
    /// Ligand residue and atom
    pub ligand: InteractingEntity,
    /// Receptor residue and atom
    pub receptor: InteractingEntity,
    /// Distance between ligand and receptor atoms
    pub distance: f64,
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
            "Ligand [{ligand}] has {intxn:?} with Receptor [{receptor}], d={dist:.3}Å",
            ligand = self.ligand,
            intxn = self.interaction,
            receptor = self.receptor,
            dist = self.distance
        )
    }
}

impl fmt::Display for Interaction {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
        // or, alternatively:
        // fmt::Debug::fmt(self, f)
    }
}
