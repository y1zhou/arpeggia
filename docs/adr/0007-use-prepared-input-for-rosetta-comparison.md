# Require prepared input for Rosetta numerical compatibility

Arpeggia targets numerical agreement with pinned Rosetta calculations only when
both programs receive the same caller-prepared full-atom structure. Adding
hydrogens, reconstructing missing atoms, assigning terminal variants, and
completing molecular chemistry remain outside the project scope. Numerical
agreement does not require method parity.

## Compatibility policy

Claims apply to a fixed canonical compatibility set and the selected default
parameters. Modified residues, ligands, solvent, and unsupported chemistry are
reported separately. Unprepared input remains calculable from observed atoms,
with conservative diagnostics and no numerical-compatibility claim. An absence
of preparation warnings does not prove chemical completeness.

Each claimed metric must independently meet overall Spearman correlation of
at least 0.99 and metric-specific absolute-error limits set before production
changes. Reports include MAE, median absolute error, RMSE, 95th-percentile and
maximum absolute error, and signed bias. Relative error, percentage change,
linear slope, and post-hoc fitting do not establish compatibility. Whole-structure
metrics order SAP candidates; per-residue absolute errors are also reported,
but per-residue Spearman correlation is not an acceptance requirement.

Candidates are deterministic definitions, not trained models, so the fixed
corpus is evaluated as one set rather than split into training and holdout data.
SASA partitions must reconcile with one surface, and SAP's atom population
must match its maximum-area calibration with recorded upstream provenance.
Empirical scaling, per-residue fitted corrections, and generated replacement
calibrations are prohibited. Among passing definitions, choose the simplest,
then fastest, then lowest-median-error candidate. Runtime and peak memory are
reported for review without an automatic acceptance ceiling.

SASA and SAP choose their atom populations independently because their
observables and calibrations differ. If the lean allowed definitions fail the
gate, report the residual error and narrow the claim rather than port Rosetta's
legacy method or relax the gate. Experimental combinations are not public API
without a separate demonstrated scientific need.

## Evidence and selected definitions

The [v0.9.0 validation report](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/v0.9.0-validation.html) records the
reference comparisons, file hashes, atom counts, parameters, warnings, timings,
and reproducibility details. Rosetta revision
`597b55d6600c3939574ffee30a4469b26c3337bd` is the documented reference, not a
production dependency or a version encoded into calculation behavior.

The experiment compared heavy/all-atom populations, ProtOr/Reduce radii, and
100/162/500 sampling points for SASA at a 1.4 Å probe. SAP compared heavy/all-atom
populations, elemental/Reduce radii, and 1.4/1.1 Å probes with matched upstream
maximum-area calibrations. Reference preparation imported and wrote a full-atom
pose without relaxation or coordinate optimization; both programs reloaded the
same saved structure.

Selected SAP uses all supplied atoms, Reduce SASA radii, a 1.1 Å probe,
Rosetta's precise hydrophobicity constants, and maximum side-chain areas from
`SapDatabase::generate_max_sasa()`. It excludes `OXT` from side-chain membership
and normalizes MSE to methionine. No benchmark fitting or Rosetta database
parser is shipped.

Selected standard SASA retains heavy-atom Shrake–Rupley at 100 points and a
1.4 Å probe with ProtOr radii and elemental fallback. The validation report
quantifies its lower total, polar, and hydrophobic SASA errors on 85 identically
prepared structures; it does not establish numerical equivalence to Rosetta's
LeGrand method. [ADR 0006](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0006-unify-standard-sasa-records.md) records the shared
atom-record and polarity-partition contracts.

One final human-readable report per benchmark retains its provenance and
acceptance evidence. Generated structures, machine-result dumps, intermediate
reports, and experimental implementations remain local or in feature-branch
history, except deliberately selected regression fixtures. Final reports are
documentation and excluded from Rust and Python packages.
