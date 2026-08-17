# Unify standard SASA records

Standard atom, residue, and chain SASA use one canonical selected atom population,
ProtOr radii with documented elemental fallback, and one rust-sasa Shrake–Rupley
kernel result per molecular context. Residue and chain values aggregate those
same atom records with a small floating-point tolerance. Selected non-solvent,
non-ion heavy atoms participate in occlusion and atom output, while polymer-level
aggregations include polymer atoms unless callers explicitly request non-polymer
groups.

Atom SASA is partitioned with Rosetta `SasaFilter` polarity while retaining the
rust-sasa Shrake–Rupley numerical method. Atom output carries `polarity` as
`polar`, `hydrophobic`, or `unknown`; residue and chain output carry `sasa`,
`polar_sasa`, `hydrophobic_sasa`, and `unclassified_sasa`. Unsupported chemical
identities contribute to `unclassified_sasa` with a structured warning. The same
partition produces `polar_dsasa`, `hydrophobic_dsasa`, and
`unclassified_dsasa` under Arpeggia's two-sided dSASA convention. The misleading
residue-level `is_polar` field is removed.

Partition sums and residue/chain aggregation must reproduce their underlying atom
SASA within a narrow floating-point tolerance. Numerical similarity to Rosetta's
LeGrand implementation is measured on a small canonical benchmark and reported
as an error distribution; it is not guaranteed, inferred from correlation, or
enforced by broad golden-test tolerances.

SAP deliberately retains its separately calibrated elemental-radius exposure
model until direct Rosetta benchmarks support changing it. The residue SAP fix
preserves positive-only score accumulation while aggregating complete side-chain
SASA and retaining eligible zero/nonpositive-score residues.
