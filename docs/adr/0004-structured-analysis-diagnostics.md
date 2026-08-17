# Return structured scientific diagnostics

Fallible Rust calculations return an analysis value containing both the result
and stable, machine-readable warning codes. The CLI renders those warnings to
stderr, while Python emits standard warnings and continues returning its existing
convenient value types. This makes automatic conformer selection, missing donor
hydrogens, unsupported topology, and recoverable parser problems visible without
turning scientifically usable partial results into failures.

Contact output also distinguishes evidence levels and physical regions:
file-backed `Covalent`, distance-backed `PotentialCovalent`, `StericClash`,
`VanDerWaalsClash`, and `VanDerWaalsContact`. These category changes are treated
as a breaking API change rather than hidden behind misleading aliases.
