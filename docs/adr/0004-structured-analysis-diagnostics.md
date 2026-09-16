# Return structured scientific diagnostics

Structure and pairwise-sequence calculations return an analysis value containing the result
and stable, machine-readable warning codes. The CLI renders those warnings to
stderr, while Python emits standard warnings and continues returning its existing
convenient value types. This makes conformer selection, missing donor
hydrogens, unsupported topology, and recoverable parser problems visible without
turning scientifically usable partial results into failures.

Antibody results retain diagnostics on their result objects instead; their
display and imputation behavior is specified in
[ADR 0010](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md).

Contact evidence categories follow
[ADR 0003](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0003-evidence-based-structure-preparation.md).

DataFrame identifier columns retain the compact schema adopted for memory and
serialized-file efficiency: model and atom identifiers are unsigned 32-bit
integers, while residue identifiers are signed 32-bit integers because residue
numbers can be negative. These ranges are sufficient for regular structure
files. An identifier outside them produces a Calculation Failure instead of
panicking or silently truncating the value.
