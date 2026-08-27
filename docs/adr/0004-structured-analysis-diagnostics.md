# Return structured scientific diagnostics

Fallible Rust calculations return an analysis value containing both the result
and stable, machine-readable warning codes. The CLI renders those warnings to
stderr, while Python emits standard warnings and continues returning its existing
convenient value types. This makes automatic conformer selection, missing donor
hydrogens, unsupported topology, and recoverable parser problems visible without
turning scientifically usable partial results into failures.

Contact output also distinguishes evidence levels and physical regions:
file-backed `Disulfide` and `Covalent`, geometry-backed `PotentialDisulfide`,
distance-backed `PotentialCovalent`, `StericClash`, `VanDerWaalsClash`, and
`VanDerWaalsContact`. Bond-related labels follow that evidence order: resolved
disulfide declarations, other resolved connectivity, qualifying undeclared CYS
geometry, then generic covalent-distance inference. These category changes are
treated as a breaking API change rather than hidden behind misleading aliases.

DataFrame identifier columns retain the compact schema adopted for memory and
serialized-file efficiency: model and atom identifiers are unsigned 32-bit
integers, while residue identifiers are signed 32-bit integers because residue
numbers can be negative. These ranges are sufficient for regular structure
files. An identifier outside them produces a Calculation Failure instead of
panicking or silently truncating the value.
