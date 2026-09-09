# Arpeggia Scientific Model

Arpeggia provides reference-compatible structural biology calculations while
making intentional deviations explicit and reviewable.

## Language

**Reference-Compatible**:
A calculation whose scientific semantics follow its named upstream
implementation or published method unless an Accepted Deviation says otherwise.
_Avoid_: Inspired by, approximately compatible

**Accepted Deviation**:
An intentional departure from a reference method supported by a recorded
maintainer rationale or an unmistakably intentional code comment.
Tests, historical persistence, and accidental compatibility do not establish one.
_Avoid_: Legacy behavior, existing behavior

**Scientific Correction**:
A change that restores reference semantics or removes behavior unsupported by an
Accepted Deviation. It may change numerical or tabular output in a minor release
when the change and its expected effects are documented.
_Avoid_: Breaking change, cleanup

**Feature Reference**:
The named authoritative method for one public calculation. Published methods
govern scientific semantics; their official implementations resolve details the
publication leaves unspecified.
_Avoid_: Global reference, approximate inspiration

**Scientific Warning**:
An unresolved doubt about the scientific validity of retained behavior. It marks
an open question requiring confirmation and never establishes an Accepted
Deviation.
_Avoid_: Accepted limitation, intentional behavior

**Unsupported Classification**:
A scientific label whose required evidence is absent from the input. Arpeggia
does not present such a label as a definitive result.
_Avoid_: Best guess, heuristic result

**Covalent Contact**:
A non-disulfide atomic contact supported by bonding information extracted from
the input structure. Its public interaction label is `Covalent`.
_Avoid_: Distance-inferred bond

**Disulfide Contact**:
A CYS SG--SG bond supported by a PDB `SSBOND` or mmCIF disulfide declaration.
Its public interaction label is `Disulfide`.
_Avoid_: Geometry-inferred disulfide

**Potential Covalent Contact**:
An atomic contact that satisfies Arpeggia's covalent-distance criterion but lacks
bonding evidence in the input structure. Its public interaction label is
`PotentialCovalent`.
_Avoid_: Covalent bond, inferred covalent bond

**Potential Disulfide Contact**:
A CYS SG--SG contact in Arpeggia's covalent-distance band whose absolute
CB--SG--SG--CB dihedral is 60--120 degrees, but whose residue pair lacks an
`SSBOND` or mmCIF disulfide declaration. Its public interaction label is
`PotentialDisulfide`.
_Avoid_: Disulfide bond, inferred disulfide bond

**Resolved Explicit Bond**:
An input bond declaration that matches two atoms in the selected model and
Selected Conformer. Contact calculation represents it with compact selected-atom
identities while retaining the qualified declaration as its evidence.
_Avoid_: Distance-inferred bond, unqualified bond

**Peptide-Adjacent Residues**:
Two residues in one selected chain whose peptide C--N geometry is continuous and
not separated by an explicit chain break. Coordinate order or consecutive
residue numbering alone does not establish adjacency.
_Avoid_: Consecutive coordinate residues, neighboring residue numbers

**Valid Explicit Chain Break**:
A PDB `TER` record placed after a complete residue and identifying that same
residue by name, chain, sequence number, and insertion code. A break inside a
residue or identifying a different residue is invalid input, not a boundary.
_Avoid_: Any TER record, implicit coordinate gap

**Van der Waals Clash**:
A nonbonded atomic contact inside the combined van der Waals envelope but outside
the severe steric-clash and potential-covalent regions.
_Avoid_: Van der Waals contact

**Van der Waals Contact**:
A nonbonded atomic contact inside the compensated outer van der Waals envelope
without overlapping the unmodified envelope.
_Avoid_: Van der Waals clash

**Observed Sequence**:
The ordered polymer residues present in atomic coordinates for a selected model
and chain. Model `0` selects the first model; solvent, ligands, and polymer
residues without coordinates are omitted.
_Avoid_: Complete sequence, declared sequence, all-model sequence

**Declared Sequence**:
The polymer sequence declared by PDB `SEQRES` or mmCIF entity metadata, including
polymer residues without coordinates.
_Avoid_: Observed sequence, coordinate sequence

**Selected Conformer**:
The single alternate conformer Arpeggia chooses for a residue: the unambiguous
blank/default conformer, otherwise highest occupancy with deterministic `A`
tie-breaking. Automatic selection produces an alternate-conformer warning.
_Avoid_: First conformer, all conformers

**Protonation Mode**:
The evidence policy used to interpret histidine charge: `AllCharged` preserves
Arpeggio's unconditional positive-ionisable typing, `Heuristic` combines explicit
evidence with a documented pH-dependent estimate, and `ExplicitOnly` never
guesses when input evidence is absent.
_Avoid_: Legacy mode, strict mode

**Histidine Tautomer Evidence**:
An explicit HID/HSD, HIE/HSE, or HIP/HSP identity, or the corresponding HD1 and
HE2 hydrogen pattern, that determines histidine H-bond donor and acceptor atoms.
Hydrogen-free plain HIS remains ambiguous and cannot prove a geometric hydrogen
bond.
_Avoid_: Histidine charge mode, all-ring donor/acceptor typing

**Rosetta-Partitioned SASA**:
Arpeggia's Shrake–Rupley atom SASA partitioned by Rosetta `SasaFilter` atom
polarity. It reproduces the polarity scheme, not Rosetta's LeGrand SASA method.
_Avoid_: Rosetta SASA, LeGrand-compatible SASA

**Prepared Full-Atom Structure**:
A caller-supplied structure containing the hydrogens, heavy atoms, termini, and
residue variants required by the selected reference calculation. Arpeggia does
not reconstruct missing chemistry for numerical-compatibility calculations.
_Avoid_: Raw structure, automatically completed structure

**Rosetta-Numeric Compatibility**:
Agreement with a pinned Rosetta calculation within the accepted error bounds
when both programs receive the same Prepared Full-Atom Structure. It does not
require algorithm or implementation parity.
_Avoid_: Rosetta method parity, approximate correlation

**Canonical Compatibility Set**:
A fixed set of Prepared Full-Atom Structures containing supported canonical
proteins, used to evaluate Rosetta-Numeric Compatibility claims.
_Avoid_: Training set, heterogeneous corpus, holdout set

**Prepared-Input Diagnostic**:
A conservative warning that detects obvious evidence of absent or inconsistent
full-atom preparation without claiming to prove chemical completeness.
_Avoid_: Preparation validation, structure completion

**Calculation Failure**:
A validly requested analysis for which Arpeggia cannot produce a scientifically
meaningful complete value. Rust returns a typed error, Python raises an
exception, and the CLI exits unsuccessfully; no null or partial scalar is
reported.
_Avoid_: Missing result, nullable score, partial success

**Sampled SC Interface**:
The paired buried surface-dot populations that remain after SC peripheral-band
trimming and can produce both directional scores. If either population is empty,
there is no complete SC value.
_Avoid_: Raw surface proximity, zero SC

**Definition-Derived Calculation**:
A deterministic calculation whose empirical parameters come from a named
scientific or reference definition rather than being fitted to benchmark output.
_Avoid_: First-principles calculation, trained model

**Polar SASA**:
The sum of SASA assigned to Rosetta donor, acceptor, or polar-hydrogen atoms.
_Avoid_: SASA of polar residues

**Hydrophobic SASA**:
The sum of SASA assigned to atoms outside Rosetta's polar atom classes when their
Rosetta-compatible classification is known.
_Avoid_: SASA of hydrophobic residues

**Unclassified SASA**:
The sum of SASA for atoms whose Rosetta polarity cannot be reproduced from the
available chemical identity and topology.
_Avoid_: Hydrophobic SASA, unknown total

**Structure Observation**:
One selected coordinate model from one input structure that contributes one
object to an RMSD comparison or structure clustering calculation.
_Avoid_: Input file, ensemble

**Structure Selection**:
The recognized amino-acid atoms retained from a Structure Observation by its
chain-residue ranges and atom subset.
_Avoid_: Atom filter, alignment

**Superposition Selection**:
The paired atoms eligible to determine the rigid-body transform between two
Structure Observations. Refinement may exclude outlying pairs from fitting.
_Avoid_: Fit selection, alignment selection, sequence alignment

**RMSD Selection**:
The paired atoms whose residual distances are evaluated under the transform
determined by the Superposition Selection and any refinement. Rejection from
fitting does not remove a pair from this evaluation set.
_Avoid_: Fit selection, superposition selection

**Exact Atom Correspondence**:
A one-to-one pairing in which two Structure Selections contain the same atom
identities, including chain, author residue number, insertion code, residue
name, and atom name.
_Avoid_: Common atoms, atom intersection

**Sequence Alignment**:
A scored correspondence between two amino-acid sequences, including paired
residues, gaps, and any unaligned terminal segments.
_Avoid_: Superposition, antibody numbering

**Semi-Global Sequence Alignment**:
An alignment consuming the entire second sequence while allowing unaligned
terminal segments of the first sequence without penalty.
_Avoid_: Symmetric overlap, local alignment

**Alignment Operation**:
A match, substitution, insertion, or deletion along a chosen alignment, directed
from the reference sequence to the query sequence. Unaligned terminal segments
are outside these operations.
_Avoid_: Minimum edit script, sequence edit distance

**Sequence Identity**:
The identical-residue pair count divided by either alignment-column count
(including gaps) or shorter full input length, with the denominator named.
_Avoid_: Substitution score, sequence coverage

**Paired-Residue Coverage**:
The nongap residue-pair count, including substitutions, divided by either
alignment-column count or shorter full input length, with the denominator named.
_Avoid_: Sequence identity, aligned span

**Sequence Edit Distance**:
The minimum number of single-residue substitutions, insertions, and deletions
needed to convert one complete input sequence into the other.
_Avoid_: Protein alignment score, edits along a chosen alignment

**Chain Correspondence**:
The pairing of chains between two Structure Observations within which residue
correspondence is established.
_Avoid_: Chain order, equal chain identifiers

**Sequence-Derived Residue Correspondence**:
Nongap residue pairs from an alignment of observed chain sequences, including
substitutions. Author numbering can differ between paired residues.
_Avoid_: Identical residues, equal residue numbers, antibody numbering

**Core RMSD**:
The RMSD of retained fitting pairs under their final fitted transform.
_Avoid_: Full-selection RMSD, final evaluation RMSD

**Evaluation RMSD**:
The RMSD of the RMSD Selection under the final fitting transform, including
selected pairs rejected from fitting.
_Avoid_: Core RMSD, refitted evaluation region

**Refinement Cycle**:
One inspection of current fitting-pair residuals, rejection of outliers, and
refitting of survivors. The initial fit precedes all refinement cycles.
_Avoid_: Initial fit, sequence realignment, numerical solver iteration

**Medoid Structure**:
An observed structure selected as a cluster representative because it minimizes
the cluster method's dissimilarity objective. It is an input structure, not an
average coordinate model.
_Avoid_: Centroid, average structure

**Pairwise RMSD Matrix**:
The symmetric dissimilarity values produced by comparing every pair of Structure
Observations with shared Superposition Selection and RMSD Selection policies.
_Avoid_: RMSD features, coordinate matrix

**Pairwise RMSD Table**:
The long-form interchange representation of a Pairwise RMSD Matrix, with one
unordered structure pair per row and columns `id_1`, `id_2`, and `rmsd`. It
excludes diagonal and reverse-duplicate rows.
_Avoid_: Wide RMSD matrix, square RMSD table

**Pairwise RMSD Cache**:
A Pairwise RMSD Table at the exact CLI output path that may replace RMSD
recalculation when its complete ID set matches the current inputs. Cache reuse
does not verify source coordinates, Superposition Selection settings, or RMSD
Selection settings.
_Avoid_: Verified RMSD cache, provenance-matched matrix

**Fixed Cluster Count**:
A cluster count specified by the caller before clustering begins.
_Avoid_: Selected cluster count, automatic clustering

**Bounded Automatic Cluster Count**:
A cluster count selected by the clustering method within caller-supplied lower
and upper bounds.
_Avoid_: Parameter-free clustering, unrestricted cluster count
