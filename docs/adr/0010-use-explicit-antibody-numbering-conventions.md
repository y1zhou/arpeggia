# Use explicit antibody numbering conventions

Arpeggia will accept `chothia` as an alias for Martin/enhanced Chothia numbering,
as requested for the antibody-numbering API. This deliberately favors the
structurally corrected convention over compatibility with historical Chothia
outputs. The alias is an Arpeggia API choice: the original conventions remain
distinct in the [scheme authors' numbering service](https://www.bioinf.org.uk/abs/abnum/).

The initial CLI and Python APIs accept named amino-acid strings. A
`NumberedAntibody` represents one variable domain, retaining its input sequence
and domain span. Terminal tags and constant-region tails are allowed;
additional detected variable domains produce an error. Structure-file
integration and constant-region numbering are deferred. Recognizable partial
domains can return results with coverage and diagnostics; numbering itself
does not fill missing sequence.

Germline comparison reports separate V and J reference similarities, coverage
and tied references. It does not reconstruct a unique ancestral antibody or
infer a D segment. References ship offline as a versioned, attributed IMGT
subset covering human and mouse H/K/L and alpaca heavy-chain/VHH references.
Search covers all bundled species by default and accepts an explicit species
restriction. Report matched-reference species rather than presumed input origin.

`cdr_definition="auto"` selects the CDR convention associated with the chosen
numbering scheme. A caller choosing a different CDR definition must explicitly
provide both arguments. This keeps the usual scheme/region pairing convenient
while making a mixed convention intentional. IMGT and Kabat use their matching
definitions. Martin uses AbM/Martin boundaries, as in the
[Martin group's study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10939163/).
AHo uses CDR1 25–40, CDR2 58–77 and CDR3 109–137 from the published
[structural-loop convention](https://pubs.rsc.org/en/content/articlehtml/2019/me/c9me00021f).
Public docstrings must cite these definitions. The meaning of an explicit
`chothia` CDR-definition argument remains open.

An explicit imputation method returns a new `NumberedAntibody` and preserves
the original object and supplied sequence. It fills only missing chunks at the
beginning of FR1 or the end of FR4, using the closest aligned germline references.
It does not fill internal gaps or replace unknown input residues. A position
is filled only when tied references agree on its presence and amino acid;
otherwise it remains unresolved unless the caller selects a reference explicitly.
Record each inferred residue and its germline source so later alignments can
distinguish supplied and inferred residues.

Engine qualification compares RIOT, Immunum and AntPack on the
[AntPack test set](https://github.com/jlparkI/AntPack/tree/main/tests/test_data).
Pin the dataset and tool versions, retain failures, and distinguish agreement
with fixture labels from independent numbering accuracy. Report residue
numbering, domain coverage and failures separately.
The [comparison report](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md)
records the completed 26,365-input run and unresolved qualification issues.

Use Immunum's native Rust core with default features disabled. Arpeggia owns
the result types, Python bindings and CLI/display formatting. Its four required
schemes and fixture agreement support this choice without establishing
independent accuracy. Integration must correct the conversion/span defect,
require meaningful aligned-domain coverage, and qualify long insertions and
multiple-domain detection. Pin the qualified core version and record any
upstream corrections.

The remaining result contract is open. Supporting evidence is in the
[numbering research](https://github.com/y1zhou/arpeggia/blob/master/docs/research/antibody-numbering-schemes-and-tools.md).
