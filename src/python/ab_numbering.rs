use super::{display_options, python_error, value_enum};
use crate::{AntibodyAlignment, NumberedAntibody, NumberedPosition, NumberedResidue};
use pyo3::prelude::*;
use std::ffi::CString;

fn warnings(py: Python<'_>, diagnostics: &[String]) -> PyResult<()> {
    for diagnostic in diagnostics {
        let message = CString::new(diagnostic.as_str()).map_err(|_| {
            pyo3::exceptions::PyValueError::new_err("diagnostic contains a NUL byte")
        })?;
        PyErr::warn(
            py,
            &py.get_type::<pyo3::exceptions::PyUserWarning>(),
            &message,
            1,
        )?;
    }
    Ok(())
}

/// Number one antibody variable domain and find its closest bundled V/J references.
///
/// Args:
///     sequence (str): Unaligned amino-acid string, 30–10000 residues. Lowercase
///         is normalized; B/Z/X/U/O are accepted. Gaps, whitespace and stops fail.
///         Tags and constant tails are retained, but a second detected domain fails.
///         Partial inputs must retain the core from FR1 through FR4.
///     name (str): Display name, default "Seq001".
///     scheme (str | None): "imgt" (default), "martin", "aho", or "kabat".
///         "chothia" is an alias for Martin/enhanced Chothia numbering.
///     cdr_definition (str): "auto" (default) follows the numbering scheme.
///         Explicit "imgt", "martin", "aho", "kabat", or "chothia" requires
///         an explicit scheme. Martin uses AbM; explicit Chothia uses the distinct
///         2021 consensus boundaries. Mixed conventions retain the requested labels.
///     species (str | Sequence[str] | None): Restrict references to "human",
///         "mouse", or "alpaca", or a sequence of these names. None searches all.
///         Alpaca references cover heavy chains; species describe matched references.
///
/// Returns:
///     NumberedAntibody: Read-only residues, zero-based input indices and half-open
///         domain_span, chain/conventions, FR/CDR sequence properties, confidence,
///         diagnostics, and separate v_match/j_match similarities. Each match's
///         hits retain exact score ties and their source references. Known-residue
///         evidence excludes ambiguous symbols. Missing evidence produces None and
///         a warning. Numbering does not fill missing residues; use result.impute().
///
/// Raises:
///     ValueError: Invalid input or options.
///     RuntimeError: No supported domain, multiple domains, severe truncation,
///         or numbering beyond the supported single-letter insertion limits.
///
/// Examples:
///     >>> antibody = arpeggia.number_antibody(sequence, name="VHH", species="alpaca")
///     >>> antibody.cdr3
///     >>> print(antibody.format(width=80, color="never"))
///
/// Conventions and sources:
/// https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md
/// IMGT: https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
/// Martin/AbM: https://pmc.ncbi.nlm.nih.gov/articles/PMC10939163/
/// AHo loops: https://pubs.rsc.org/en/content/articlehtml/2019/me/c9me00021f
#[pyfunction]
#[pyo3(signature=(sequence, *, name="Seq001", scheme=None, cdr_definition="auto", species=None))]
fn number_antibody(
    py: Python<'_>,
    sequence: &str,
    name: &str,
    scheme: Option<&str>,
    cdr_definition: &str,
    species: Option<&Bound<'_, PyAny>>,
) -> PyResult<NumberedAntibody> {
    let species = species
        .map(|s| {
            s.extract::<String>()
                .map(|s| vec![s])
                .or_else(|_| s.extract::<Vec<String>>())
        })
        .transpose()?
        .unwrap_or_default();
    let options = crate::NumberingOptions {
        name: name.into(),
        scheme: scheme
            .map(|s| value_enum(s, "scheme must be imgt, martin (or chothia), aho, or kabat"))
            .transpose()?,
        cdr_definition: value_enum(
            cdr_definition,
            "cdr_definition must be auto, imgt, martin, aho, kabat, or chothia",
        )?,
        species: species
            .iter()
            .map(|s| value_enum(s, "species must be human, mouse, or alpaca"))
            .collect::<PyResult<Vec<_>>>()?,
    };
    let result = py
        .detach(|| crate::number_antibody(sequence, &options))
        .map_err(python_error)?;
    warnings(py, &result.diagnostics)?;
    Ok(result)
}

/// Align numbered antibodies through their shared position labels.
///
/// Args:
///     antibodies (Sequence[NumberedAntibody]): One or more numbered antibodies,
///         all in the same scheme and all heavy or all light chains. K/L mixtures
///         and different per-row CDR definitions are allowed.
///     reference_index (int): Zero-based input row used for comparisons, default 0.
///
/// Returns:
///     AntibodyAlignment: Read-only antibodies, positions and aligned_sequences
///         in input order. Displays place the selected reference first and hide
///         germline rows. The grid has no multiple-sequence-alignment score.
///
/// Raises:
///     ValueError: Empty/incompatible inputs or an invalid reference index.
///
/// Examples:
///     >>> comparison = arpeggia.align_antibodies([wild_type, mutant], reference_index=0)
///     >>> print(comparison.format(width=80, rulers=False))
#[pyfunction]
#[pyo3(signature=(antibodies, *, reference_index=0))]
fn align_antibodies(
    py: Python<'_>,
    antibodies: Vec<Py<NumberedAntibody>>,
    reference_index: usize,
) -> PyResult<AntibodyAlignment> {
    let antibodies = antibodies.iter().map(|a| a.borrow(py).clone()).collect();
    py.detach(|| crate::align_antibodies(antibodies, reference_index))
        .map_err(python_error)
}

#[pymethods]
impl NumberedAntibody {
    /// Numbered sequence, including explicit imputation and excluding input tails.
    #[getter(sequence)]
    fn sequence_python(&self) -> String {
        self.sequence()
    }
    /// First framework sequence under the selected CDR definition.
    #[getter(fr1)]
    fn fr1_python(&self) -> String {
        self.fr1()
    }
    /// First CDR sequence under the selected CDR definition.
    #[getter(cdr1)]
    fn cdr1_python(&self) -> String {
        self.cdr1()
    }
    /// Second framework sequence under the selected CDR definition.
    #[getter(fr2)]
    fn fr2_python(&self) -> String {
        self.fr2()
    }
    /// Second CDR sequence under the selected CDR definition.
    #[getter(cdr2)]
    fn cdr2_python(&self) -> String {
        self.cdr2()
    }
    /// Third framework sequence under the selected CDR definition.
    #[getter(fr3)]
    fn fr3_python(&self) -> String {
        self.fr3()
    }
    /// Third CDR sequence under the selected CDR definition.
    #[getter(cdr3)]
    fn cdr3_python(&self) -> String {
        self.cdr3()
    }
    /// Final framework sequence under the selected CDR definition.
    #[getter(fr4)]
    fn fr4_python(&self) -> String {
        self.fr4()
    }

    /// Return a new antibody with supported missing FR1/FR4 terminal residues filled.
    ///
    /// Args:
    ///     v_reference (str | None): Exact GermlineReference.id from a tied V hit.
    ///         None requires every tied reference to agree on presence and residue.
    ///     j_reference (str | None): Equivalent selector for tied J hits.
    ///
    /// Returns:
    ///     NumberedAntibody: Original input and domain_span are preserved. Added
    ///         residues have input_index=None and imputed_from reference IDs.
    ///         Internal gaps, CDRs and unknown input residues remain unchanged.
    ///         Unsupported or conflicting reference evidence stays unresolved.
    ///
    /// Raises:
    ///     ValueError: A selector is absent from the tied reference set.
    #[pyo3(name="impute",signature=(*, v_reference=None,j_reference=None))]
    fn impute_python(
        &self,
        py: Python<'_>,
        v_reference: Option<&str>,
        j_reference: Option<&str>,
    ) -> PyResult<Self> {
        let result = py
            .detach(|| self.impute(v_reference, j_reference))
            .map_err(python_error)?;
        warnings(py, &result.diagnostics[self.diagnostics.len()..])?;
        Ok(result)
    }

    /// Render the input above its combined V/J germlines, with CDR highlighting.
    ///
    /// Each block contains one CDR marker row, input ruler and sequence, germline
    /// ruler and sequence, then operations relative to the input. CDR1/2/3 bands
    /// are gray/pink/cyan across every row except operations. Yellow backgrounds mark
    /// imputed residues and their summary count. All tied reference names appear
    /// in the summary; only representative V/J sequences are shown.
    ///
    /// Args:
    ///     width (int | None): Total visible columns including both gene labels.
    ///         None detects terminal width or falls back to 80.
    ///     color (str): "auto", "always", or "never". Auto uses sys.stdout and
    ///         respects NO_COLOR. Stored data remain plain.
    ///     rulers (bool): Show one-based residue positions at every tenth residue.
    ///         Supplied-input coordinates survive imputation; imputed positions
    ///         are blank. Germline counts continue from V through J, ignoring gaps. False
    ///         hides rulers but retains endpoint numbers and CDR markers.
    ///
    /// Returns:
    ///     str: Wrapped sequence rows and separate V/J similarities for the supplied
    ///         input. Outer germline padding is blank; unknown junction hyphens
    ///         are gray with blank operations. A '+' marks a germline insertion,
    ///         '-' a deletion, ':' a positive-BLOSUM62 substitution, and 'x' any
    ///         other mismatch. Matches have blank operations.
    ///
    /// Raises:
    ///     ValueError: Invalid color policy or insufficient width.
    #[pyo3(name="format",signature=(width=None,color="auto",rulers=true))]
    fn format_python(
        &self,
        py: Python<'_>,
        width: Option<usize>,
        color: &str,
        rulers: bool,
    ) -> PyResult<String> {
        let (width, color) = display_options(py, width, color)?;
        self.render(width, color, rulers).map_err(python_error)
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.format_python(py, None, "auto", true)
    }
    fn __str__(&self, py: Python<'_>) -> PyResult<String> {
        self.__repr__(py)
    }
}

#[pymethods]
impl AntibodyAlignment {
    /// Render antibody comparisons with the selected reference first and germlines hidden.
    ///
    /// The selected reference defines the shared CDR bands and the convention
    /// summary. Each row retains its own original input coordinates. Imputed
    /// residues have yellow backgrounds; the summary count totals all antibodies.
    ///
    /// Args:
    ///     width (int | None): Total visible width, detected from the terminal
    ///         when omitted, with an 80-column fallback.
    ///     color (str): "auto" (default), "always", or "never".
    ///     rulers (bool): Show one-based original input positions at every tenth
    ///         residue. Imputed positions are blank. False hides rulers while
    ///         retaining endpoint numbers and the single CDR marker row per block.
    ///     reference_index (int | None): Override the zero-based comparison row
    ///         for this display. None uses the object's stored reference_index.
    ///
    /// Returns:
    ///     str: Wrapped comparisons. Stored antibodies, positions, aligned_sequences
    ///         and reference_index remain unchanged.
    ///
    /// Raises:
    ///     ValueError: Invalid reference index, color policy, or width.
    #[pyo3(name="format",signature=(width=None,color="auto",rulers=true,*,reference_index=None))]
    fn format_python(
        &self,
        py: Python<'_>,
        width: Option<usize>,
        color: &str,
        rulers: bool,
        reference_index: Option<usize>,
    ) -> PyResult<String> {
        let (width, color) = display_options(py, width, color)?;
        self.render(width, color, rulers, reference_index)
            .map_err(python_error)
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.format_python(py, None, "auto", true, None)
    }
    fn __str__(&self, py: Python<'_>) -> PyResult<String> {
        self.__repr__(py)
    }
}

#[pymethods]
impl NumberedPosition {
    fn __str__(&self) -> String {
        self.to_string()
    }
    fn __repr__(&self) -> String {
        format!("NumberedPosition({:?})", self.to_string())
    }
}

#[pymethods]
impl NumberedResidue {
    fn __repr__(&self) -> String {
        format!("{}:{} ({})", self.position, self.amino_acid, self.region)
    }
}

pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<NumberedAntibody>()?;
    m.add_class::<AntibodyAlignment>()?;
    m.add_class::<NumberedPosition>()?;
    m.add_class::<NumberedResidue>()?;
    m.add_class::<crate::GermlineReference>()?;
    m.add_class::<crate::GermlineHit>()?;
    m.add_class::<crate::GermlineMatch>()?;
    m.add_function(wrap_pyfunction!(number_antibody, m)?)?;
    m.add_function(wrap_pyfunction!(align_antibodies, m)?)?;
    Ok(())
}
