//! Antibody variable-domain numbering with explicit region conventions.

mod alignment;
mod core;
mod germline;
mod impute;
pub use alignment::{AntibodyAlignment, align_antibodies};
pub use germline::{GermlineHit, GermlineMatch, GermlineReference, GermlineSpecies};

use crate::{ArpeggiaError, ArpeggiaResult};
use immunum::{Chain, Scheme};
use serde::Serialize;
use std::fmt::{Display, Formatter};

/// Supported numbering convention. Chothia resolves to Martin.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, clap::ValueEnum, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum NumberingScheme {
    /// IMGT unique numbering.
    #[default]
    Imgt,
    /// Martin/enhanced Chothia numbering.
    #[value(alias = "chothia")]
    Martin,
    /// AHo structural numbering.
    Aho,
    /// Kabat sequence numbering.
    Kabat,
}

impl NumberingScheme {
    pub(crate) fn backend(self) -> Scheme {
        match self {
            Self::Imgt => Scheme::IMGT,
            Self::Martin => Scheme::Martin,
            Self::Aho => Scheme::Aho,
            Self::Kabat => Scheme::Kabat,
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Imgt => "imgt",
            Self::Martin => "martin",
            Self::Aho => "aho",
            Self::Kabat => "kabat",
        }
    }
}

/// CDR region convention, independent of numbering labels.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, clap::ValueEnum)]
pub enum CdrDefinition {
    /// Follow the selected numbering scheme.
    #[default]
    Auto,
    /// IMGT CDR boundaries.
    Imgt,
    /// AbM boundaries in Martin numbering.
    Martin,
    /// AHo structural-loop boundaries.
    Aho,
    /// Kabat CDR boundaries.
    Kabat,
    /// Distinct 2021 consensus Chothia CDR boundaries.
    Chothia,
}

impl CdrDefinition {
    fn backend(self, scheme: NumberingScheme) -> Scheme {
        match self {
            Self::Auto => scheme.backend(),
            Self::Imgt => Scheme::IMGT,
            Self::Martin => Scheme::Martin,
            Self::Aho => Scheme::Aho,
            Self::Kabat => Scheme::Kabat,
            Self::Chothia => Scheme::Chothia,
        }
    }
}

/// Options for numbering one unaligned variable-domain sequence.
#[derive(Clone, Debug, Default)]
pub struct NumberingOptions {
    /// Display name; an empty name uses `Seq001`.
    pub name: String,
    /// Omission selects IMGT. An explicit CDR override requires an explicit scheme.
    pub scheme: Option<NumberingScheme>,
    /// Region convention; automatic follows the numbering scheme.
    pub cdr_definition: CdrDefinition,
    /// Restrict germline matching; empty searches all bundled species.
    pub species: Vec<GermlineSpecies>,
}

/// A numbered position with an optional single-letter insertion code.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Serialize)]
#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object, module = "arpeggia")
)]
pub struct NumberedPosition {
    /// Position number in the selected scheme.
    pub number: u8,
    /// Insertion letter, when present.
    pub insertion: Option<char>,
}

impl Display for NumberedPosition {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.number)?;
        if let Some(letter) = self.insertion {
            write!(f, "{letter}")?;
        }
        Ok(())
    }
}

impl From<immunum::Position> for NumberedPosition {
    fn from(position: immunum::Position) -> Self {
        Self {
            number: position.number,
            insertion: position.insertion,
        }
    }
}

impl NumberedPosition {
    pub(crate) fn order(self, scheme: Scheme, chain: Chain) -> (u8, i16) {
        let insertion = self.insertion.map_or(0, |c| i16::from(c as u8 - b'A' + 1));
        let reversed = core::rules(scheme, chain).iter().any(|rule| {
            matches!(rule.insertion, immunum::Insertion::Symmetric { right, .. } if right == self.number)
        });
        (self.number, if reversed { -insertion } else { insertion })
    }
}

/// A residue and its correspondence to the supplied antibody sequence.
#[derive(Clone, Debug, Serialize)]
#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object, module = "arpeggia")
)]
pub struct NumberedResidue {
    /// Numbered position.
    pub position: NumberedPosition,
    /// Supplied amino-acid symbol.
    pub amino_acid: char,
    /// Zero-based index in the original input; absent for imputed residues.
    pub input_index: Option<usize>,
    /// FR1, CDR1, FR2, CDR2, FR3, CDR3 or FR4.
    pub region: String,
    /// Source reference IDs for an imputed residue; empty for supplied residues.
    pub imputed_from: Vec<String>,
}

/// One variable domain, retaining the complete input and original domain span.
#[derive(Clone, Debug, Serialize)]
#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object, module = "arpeggia")
)]
pub struct NumberedAntibody {
    /// Display name.
    pub name: String,
    /// Complete normalized input, including unnumbered tails.
    pub input_sequence: String,
    /// Zero-based, half-open numbered span in the original input.
    pub domain_span: (usize, usize),
    /// Detected H, K or L chain class.
    pub chain: String,
    /// Resolved numbering convention.
    pub scheme: String,
    /// Resolved CDR convention.
    pub cdr_definition: String,
    /// Ordered residue correspondence.
    pub residues: Vec<NumberedResidue>,
    /// Profile confidence heuristic; not a probability of correct numbering.
    pub confidence: f32,
    /// Distinct profile positions with matched input residues, excluding insertions.
    pub matched_profile_positions: usize,
    /// Recoverable limitations of this annotation.
    pub diagnostics: Vec<String>,
    /// Best qualifying V similarities, or none when reference evidence is insufficient.
    pub v_match: Option<GermlineMatch>,
    /// Best qualifying J similarities, or none when reference evidence is insufficient.
    pub j_match: Option<GermlineMatch>,
}

impl NumberedAntibody {
    /// Numbered sequence, excluding unnumbered input tails.
    pub fn sequence(&self) -> String {
        self.residues.iter().map(|r| r.amino_acid).collect()
    }

    /// Residues assigned to an FR/CDR region, in input order.
    pub fn region_sequence(&self, region: &str) -> String {
        self.residues
            .iter()
            .filter(|r| r.region.eq_ignore_ascii_case(region))
            .map(|r| r.amino_acid)
            .collect()
    }

    /// First framework sequence.
    pub fn fr1(&self) -> String {
        self.region_sequence("FR1")
    }
    /// First complementarity-determining region sequence.
    pub fn cdr1(&self) -> String {
        self.region_sequence("CDR1")
    }
    /// Second framework sequence.
    pub fn fr2(&self) -> String {
        self.region_sequence("FR2")
    }
    /// Second complementarity-determining region sequence.
    pub fn cdr2(&self) -> String {
        self.region_sequence("CDR2")
    }
    /// Third framework sequence.
    pub fn fr3(&self) -> String {
        self.region_sequence("FR3")
    }
    /// Third complementarity-determining region sequence.
    pub fn cdr3(&self) -> String {
        self.region_sequence("CDR3")
    }
    /// Final framework sequence.
    pub fn fr4(&self) -> String {
        self.region_sequence("FR4")
    }
}

/// Number one antibody variable domain with explicit numbering/CDR conventions.
///
/// IMGT is the default; `chothia` is a CLI/Python alias for Martin numbering.
/// Inputs may have tags or constant tails, but a second recognized domain fails.
/// Only terminal FR1/FR4 truncations retaining the intervening core are supported.
/// Numbering does not impute missing residues.
///
/// Conventions: <https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html>,
/// <https://pmc.ncbi.nlm.nih.gov/articles/PMC10939163/>, and
/// <https://pubs.rsc.org/en/content/articlehtml/2019/me/c9me00021f>.
pub fn number_antibody(
    sequence: &str,
    options: &NumberingOptions,
) -> ArpeggiaResult<NumberedAntibody> {
    if options.scheme.is_none() && options.cdr_definition != CdrDefinition::Auto {
        return Err(ArpeggiaError::InvalidArgument(
            "an explicit cdr_definition requires an explicit scheme".into(),
        ));
    }
    core::number(sequence, options)
}
