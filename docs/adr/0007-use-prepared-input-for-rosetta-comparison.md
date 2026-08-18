# Require prepared input for Rosetta numerical compatibility

Arpeggia targets numerical agreement with pinned Rosetta calculations only when
both programs receive the same caller-prepared full-atom structure. Adding
hydrogens, reconstructing missing atoms, assigning terminal variants, and
otherwise completing molecular chemistry remain outside the project scope.
Numerical agreement does not require method parity.

The existing v0.9 implementation is preserved as a commit boundary before its
SASA and SAP calculations are refined in place. A result may be described as
Rosetta-numerically compatible only on a canonical compatibility set.
Compatibility is evaluated with absolute-error distributions and overall
Spearman rank correlation, not relative error, percentage change, linear slope,
or post-hoc fitting. Heterogeneous and unsupported structures are reported
separately.

Prepared input does not predetermine the calculation's atom population.
Arpeggia will benchmark all-atom and heavy-only surfaces independently because
its Shrake--Rupley method may not benefit numerically from copying Rosetta's
all-atom population in isolation. The implementation selected in place is the
scientifically defined, non-fitted variant that passes the compatibility gate;
its atom population is documented as part of that definition.

Unprepared or obviously hydrogen-free input remains calculable from its
observed atoms but emits a structured warning and carries no numerical
compatibility claim. Empirical output scaling, per-residue fitted corrections,
and other post-hoc corpus fitting are prohibited. Rosetta commit
`597b55d6600c3939574ffee30a4469b26c3337bd` is the documented v0.9 reference,
but the revision is not encoded in production behavior.

SASA partitioning and SAP select their heavy-only or all-atom populations
independently because their observables and calibrations differ. If the allowed
lean variants fail the gate, Arpeggia reports the residual error and narrows its
claim; it does not port Rosetta's legacy method or weaken the acceptance gate.
The compatibility claim covers supported canonical amino-acid proteins only.
Concrete OXT and MSE bookkeeping defects are corrected, but modified residues,
ligands, solvent, and other unsupported chemistry remain outside the claim.

Prepared-input diagnostics are intentionally conservative rather than chemical
validation. They scan for an entirely hydrogen-free selected model and for
histidine name/atom patterns that leave protonation unresolved or internally
inconsistent. These checks emit warnings and never add atoms, infer missing
chemistry, or claim that an input without warnings is complete.

The preparation benchmark may use the pinned Rosetta executable to import and
write a full-atom pose, after which both programs reload the identical saved
file. It records file hashes and atom counts and performs no relaxation,
minimization, or other expensive coordinate optimization. Rosetta remains a
benchmark-only tool rather than a build or runtime dependency. Histidine
preparation diagnostics reuse `ExplicitOnly` evidence without changing the
caller's selected protonation mode.

Before production calculations change, the benchmark tests a bounded factorial.
SASA crosses heavy-only and all-atom populations, ProtOr and atom-typed Reduce
radii, and 100, 162, and 500 points at a 1.4 Angstrom probe. SAP crosses both
atom populations, elemental and atom-typed Reduce radii, and 1.4 and 1.1
Angstrom probes. SAP maximum areas are taken from pinned upstream calibration
data, transformed without adjusting values, and matched to the candidate's
heavy-only or all-atom population. Existing mismatched combinations remain a
diagnostic baseline only; Arpeggia does not generate or fit a new empirical
maximum-area dataset.

Every claimed metric must pass the compatibility gate independently. Overall
Spearman correlation must be at least 0.99. The report records MAE, median
absolute error, RMSE, 95th-percentile absolute error, maximum error, and signed
bias; the maintainer sets metric-specific absolute-error limits at the review
checkpoint before production changes. Among passing definitions, Arpeggia
selects the simplest, then the faster, then the one with lower median error. The
complete matrix remains available so the selected implementation is auditable.

Runtime and peak memory are measured for every candidate but have no automatic
acceptance ceiling; the maintainer reviews their trade-off after seeing the
complete report. Total SASA may change with the selected additive surface, but
its absolute-error distribution and Spearman correlation remain first-class
acceptance evidence alongside the polar and hydrophobic partitions.

Atom population and radius schemes remain internal during the experiment. If
the complete results demonstrate that multiple schemes have distinct,
scientifically useful strengths, the v0.9 release design may expose them through
a public API; experimental combinations are not published speculatively.
Compatibility applies only to the selected default parameters.

The factorial ends at a review checkpoint before production calculations
change. Its durable evidence consists of a human-readable HTML report and a
machine-readable result file containing source and prepared hashes, atom counts,
parameters, Rosetta revision, per-structure metrics, warnings, failures,
timings, and summary statistics. Generated prepared structures stay outside Git
except for deliberately selected regression fixtures.

The canonical compatibility set is evaluated as one fixed corpus rather than
split into development and holdout subsets. Candidate formulas are deterministic
and derived from their stated scientific definitions; no benchmark values train
or fit them. Whole-structure metrics order SAP candidates. Per-residue SAP
absolute errors are reported because residue output is public, but per-residue
Spearman correlation is not an acceptance requirement.

Scientific consistency is a hard gate independent of numerical agreement. A
candidate is rejected if its SASA partitions do not reconcile with one surface,
or if its SAP atom population and maximum-area population are mismatched or the
upstream calibration provenance is missing.

Upstream Rosetta calibration data may be fetched, deterministically cleaned,
and evaluated locally during the benchmark. Neither the upstream nor cleaned
dataset is committed during this task; only the cleaner, source provenance, and
result artifacts are candidates for version control. Arpeggia itself is
GPL-3.0, and its README and Python package metadata point to the repository
`LICENSE` file.
