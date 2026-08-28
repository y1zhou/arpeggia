//! Stable analysis warnings and public errors.

use std::fmt::{Display, Formatter};

/// Machine-searchable warning categories emitted by analyses.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum WarningCode {
    /// The input parser recovered from a non-fatal problem.
    Parser,
    /// One alternate conformer was selected from several.
    ConformerSelected,
    /// The first coordinate model was selected from a multi-model input.
    ModelSelected,
    /// A donor had no directly associated explicit hydrogen.
    MissingDonorHydrogen,
    /// A full-atom reference calculation received no hydrogen atoms.
    HydrogenFreeInput,
    /// Histidine charge could not be resolved from the selected policy.
    UnresolvedHistidine,
    /// Histidine naming and explicit ring hydrogens disagree.
    InconsistentHistidine,
    /// A declared polymer monomer could not be mapped to a sequence letter.
    UnsupportedMonomer,
    /// An atom could not be classified by the selected SASA polarity scheme.
    UnsupportedPolarity,
    /// An optional geometric feature could not be built from the available atoms.
    IncompleteGeometry,
}

impl Display for WarningCode {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let code = match self {
            Self::Parser => "PARSER",
            Self::ConformerSelected => "CONFORMER_SELECTED",
            Self::ModelSelected => "MODEL_SELECTED",
            Self::MissingDonorHydrogen => "MISSING_DONOR_HYDROGEN",
            Self::HydrogenFreeInput => "HYDROGEN_FREE_INPUT",
            Self::UnresolvedHistidine => "UNRESOLVED_HISTIDINE",
            Self::InconsistentHistidine => "INCONSISTENT_HISTIDINE",
            Self::UnsupportedMonomer => "UNSUPPORTED_MONOMER",
            Self::UnsupportedPolarity => "UNSUPPORTED_POLARITY",
            Self::IncompleteGeometry => "INCOMPLETE_GEOMETRY",
        };
        f.write_str(code)
    }
}

/// A recoverable, stable diagnostic accompanying a successful result.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AnalysisWarning {
    /// Stable category intended for programmatic handling.
    pub code: WarningCode,
    /// Human-readable context.
    pub message: String,
}

impl AnalysisWarning {
    /// Construct a warning.
    pub fn new(code: WarningCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
        }
    }
}

impl Display for AnalysisWarning {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}] {}", self.code, self.message)
    }
}

/// A successful value together with recoverable warnings.
#[derive(Clone, Debug)]
pub struct Analysis<T> {
    /// Successful result.
    pub value: T,
    /// Recoverable warnings emitted while producing the result.
    pub warnings: Vec<AnalysisWarning>,
}

impl<T> Analysis<T> {
    /// Construct a successful analysis.
    pub fn new(value: T, warnings: Vec<AnalysisWarning>) -> Self {
        Self { value, warnings }
    }

    /// Transform the value while preserving warnings.
    pub fn map<U>(self, map: impl FnOnce(T) -> U) -> Analysis<U> {
        Analysis::new(map(self.value), self.warnings)
    }
}

/// Errors that prevent an analysis from producing a valid value.
#[derive(Debug)]
#[non_exhaustive]
pub enum ArpeggiaError {
    /// The input or output could not be read or written.
    Io(std::io::Error),
    /// The structure parser reported an invalidating error.
    Parse(String),
    /// A public argument was invalid.
    InvalidArgument(String),
    /// A valid request could not produce a complete scientific value.
    Calculation(String),
}

impl Display for ArpeggiaError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => write!(f, "I/O error: {error}"),
            Self::Parse(message) => write!(f, "structure parse error: {message}"),
            Self::InvalidArgument(message) => write!(f, "invalid argument: {message}"),
            Self::Calculation(message) => write!(f, "calculation failed: {message}"),
        }
    }
}

impl std::error::Error for ArpeggiaError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::Parse(_) | Self::InvalidArgument(_) | Self::Calculation(_) => None,
        }
    }
}

impl From<std::io::Error> for ArpeggiaError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

/// Public result type for invalidating failures.
pub type ArpeggiaResult<T> = Result<T, ArpeggiaError>;
