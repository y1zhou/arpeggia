use core::fmt;
use pdbtbx::*;

/// Interaction types.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Interaction {
    /// A non-disulfide bond explicitly declared by the input structure.
    Covalent,
    /// A disulfide explicitly declared by the input structure.
    Disulfide,
    /// A nonbonded overlap shorter than the potential-covalent band.
    StericClash,
    /// Within the covalent-distance band without input bond evidence
    PotentialCovalent,
    /// Cysteine SG geometry consistent with a disulfide, without input bond evidence
    PotentialDisulfide,
    /// Nonbonded overlap inside the van der Waals envelope
    VanDerWaalsClash,
    /// Within the compensated outer van der Waals envelope
    VanDerWaalsContact,

    // Electrostatic
    /// Negative and positive charged ions
    IonicBond,
    /// An ionic bond that depends on unresolved histidine protonation.
    PotentialIonicBond,
    /// Donor-H-acceptor hydrogen bonds
    HydrogenBond,
    /// C-H...O hydrogen bond
    WeakHydrogenBond,
    // HalogenBond,
    /// Hydrogen bonding without angle terms
    PolarContact,
    /// C-H...O bonding without angle terms
    WeakPolarContact,
    /// Like charges repelling each other
    IonicRepulsion,
    /// Like-charge repulsion that depends on unresolved histidine protonation.
    PotentialIonicRepulsion,
    /// Presence of both ionic and hydrogen bonding
    SaltBridge,

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
    /// Cation-pi interaction that depends on unresolved histidine protonation.
    PotentialCationPi,

    /// Hydrophobic interaction
    HydrophobicContact,
}

impl Interaction {
    pub(crate) const fn as_str(&self) -> &'static str {
        match self {
            Self::Covalent => "Covalent",
            Self::Disulfide => "Disulfide",
            Self::StericClash => "StericClash",
            Self::PotentialCovalent => "PotentialCovalent",
            Self::PotentialDisulfide => "PotentialDisulfide",
            Self::VanDerWaalsClash => "VanDerWaalsClash",
            Self::VanDerWaalsContact => "VanDerWaalsContact",
            Self::IonicBond => "IonicBond",
            Self::PotentialIonicBond => "PotentialIonicBond",
            Self::HydrogenBond => "HydrogenBond",
            Self::WeakHydrogenBond => "WeakHydrogenBond",
            Self::PolarContact => "PolarContact",
            Self::WeakPolarContact => "WeakPolarContact",
            Self::IonicRepulsion => "IonicRepulsion",
            Self::PotentialIonicRepulsion => "PotentialIonicRepulsion",
            Self::SaltBridge => "SaltBridge",
            Self::PiDisplacedStacking => "PiDisplacedStacking",
            Self::PiTStacking => "PiTStacking",
            Self::PiSandwichStacking => "PiSandwichStacking",
            Self::PiParallelInPlaneStacking => "PiParallelInPlaneStacking",
            Self::PiTiltedStacking => "PiTiltedStacking",
            Self::PiLStacking => "PiLStacking",
            Self::CationPi => "CationPi",
            Self::PotentialCationPi => "PotentialCationPi",
            Self::HydrophobicContact => "HydrophobicContact",
        }
    }
}

/// Policy for interpreting histidine protonation in ionic contacts.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, clap::ValueEnum)]
pub enum ProtonationMode {
    /// Treat every histidine as potentially charged, matching legacy Arpeggio typing.
    #[default]
    AllCharged,
    /// Resolve explicit evidence first, then use an intrinsic-pKa population heuristic.
    Heuristic,
    /// Treat histidine as charged only when the input contains explicit evidence.
    ExplicitOnly,
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
    /// residue insertion code
    pub insertion: String,
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
    /// Model identifier
    pub model: usize,
    /// Interaction type
    pub interaction: Interaction,
    /// Ligand residue and atom
    pub ligand: InteractingEntity,
    /// Receptor residue and atom
    pub receptor: InteractingEntity,
    /// Distance between ligand and receptor atoms
    pub distance: f64,
}

impl InteractingEntity {
    /// Create a new interacting entity
    pub fn new(
        chain: &str,
        resi: isize,
        insertion: &str,
        altloc: &str,
        resn: &str,
        atomn: &str,
        atomi: usize,
    ) -> Self {
        Self {
            chain: chain.to_string(),
            resi,
            insertion: insertion.to_string(),
            altloc: altloc.to_string(),
            resn: resn.to_string(),
            atomn: atomn.to_string(),
            atomi,
        }
    }
    /// Helper function to convert an [`pdbtbx::AtomConformerResidueChainModel`] to a human-readable format
    pub fn from_hier(hierarchy: &AtomConformerResidueChainModel<'_>) -> InteractingEntity {
        Self::new(
            hierarchy.chain().id(),
            hierarchy.residue().serial_number(),
            hierarchy.residue().insertion_code().unwrap_or(""),
            hierarchy.conformer().alternative_location().unwrap_or(""),
            hierarchy.residue().name().unwrap_or(""),
            hierarchy.atom().name(),
            hierarchy.atom().serial_number(),
        )
    }
}

impl fmt::Display for InteractingEntity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Chain {chain}, Residue {resn} {resi}{insertion} {altloc}, Atom {atom_name} {atom_idx}",
            chain = self.chain,
            resn = self.resn,
            resi = self.resi,
            insertion = self.insertion,
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
        f.write_str(self.as_str())
    }
}
