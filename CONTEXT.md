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
An atomic contact supported by bonding information extracted from the input
structure. Its public interaction label is `Covalent`.
_Avoid_: Distance-inferred bond

**Potential Covalent Contact**:
An atomic contact that satisfies Arpeggia's covalent-distance criterion but lacks
bonding evidence in the input structure. Its public interaction label is
`PotentialCovalent`.
_Avoid_: Covalent bond, inferred covalent bond

**Van der Waals Clash**:
A nonbonded atomic contact inside the combined van der Waals envelope but outside
the severe steric-clash and potential-covalent regions.
_Avoid_: Van der Waals contact

**Van der Waals Contact**:
A nonbonded atomic contact inside the compensated outer van der Waals envelope
without overlapping the unmodified envelope.
_Avoid_: Van der Waals clash

**Observed Sequence**:
The ordered polymer residues present in atomic coordinates for a chain. It omits
solvent, ligands, and polymer residues without coordinates.
_Avoid_: Complete sequence, declared sequence

**Declared Sequence**:
The polymer sequence declared by PDB `SEQRES` or mmCIF entity metadata, including
polymer residues without coordinates.
_Avoid_: Observed sequence, coordinate sequence

**Selected Conformer**:
The single alternate conformer Arpeggia chooses for a residue: the unambiguous
blank/default conformer, otherwise highest occupancy with deterministic `A`
tie-breaking. Automatic selection produces a Scientific Warning.
_Avoid_: First conformer, all conformers

**Protonation Mode**:
The evidence policy used to interpret histidine charge: `AllCharged` preserves
Arpeggio's unconditional positive-ionisable typing, `Heuristic` combines explicit
evidence with a documented pH-dependent estimate, and `ExplicitOnly` never
guesses when input evidence is absent.
_Avoid_: Legacy mode, strict mode

**Rosetta-Partitioned SASA**:
Arpeggia's Shrake–Rupley atom SASA partitioned by Rosetta `SasaFilter` atom
polarity. It reproduces the polarity scheme, not Rosetta's LeGrand SASA method.
_Avoid_: Rosetta SASA, LeGrand-compatible SASA

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
