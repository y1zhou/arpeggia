# Cleanup and package-size audit

Audited against commit `7e06c82549c1eb94a1e75b2e1d5c2d38f426bc99`
(v0.9.1 code), comparing release metadata with v0.9.0.

## Package-size findings

Published compressed artifact sizes, in decimal MB:

| Artifact | v0.9.0 | v0.9.1 |
| --- | ---: | ---: |
| CPython 3.13 Linux x86-64 wheel | 6.42 | 21.65 |
| Linux x86-64 GNU CLI archive | 9.37 | 24.65 |

Sources: [v0.9.0 PyPI metadata](https://pypi.org/pypi/arpeggia/0.9.0/json),
[v0.9.1 PyPI metadata](https://pypi.org/pypi/arpeggia/0.9.1/json),
[v0.9.0 GitHub release](https://github.com/y1zhou/arpeggia/releases/tag/v0.9.0),
and [v0.9.1 GitHub release](https://github.com/y1zhou/arpeggia/releases/tag/v0.9.1).

The major avoidable growth comes from enabling Polars' `lazy` feature for a
single NDJSON reader. Combined with JSON/Parquet support, it activates the
query planner, execution engines, streaming engine, and cloud/network support.
Removing `lazy` reduces the normal dependency graph from 311 to 196 package
nodes. Parquet support remains enabled because both reading and writing it
are supported features. RMSD/clustering also legitimately add `kmedoids` and
`sysinfo`; those remain in use.

The replacement uses Polars' eager JSON reader, infers and projects the schema
before materializing columns, and checks expected cache row counts in a
preliminary buffered pass. It stops at the first excess nonblank row, rejects
wrong heights before allocating a DataFrame, and preserves input validation.
Valid NDJSON caches incur an extra sequential read; this avoids carrying a
query engine solely for its row-limit option. A subsequent
[5B8C contacts benchmark](benchmarks/ndjson-5b8c.md) found essentially unchanged
contact-generation timings and unchanged or faster NDJSON loading on a
2,574-row table.

Same-machine release wheel builds with the same CPython 3.13 interpreter,
Rust toolchain, lockfile, and maturin settings measured:

| Content | Before cleanup | After cleanup |
| --- | ---: | ---: |
| Compressed wheel | 22,376,505 bytes | 9,480,021 bytes |
| Uncompressed extension | 72,859,664 bytes | 34,657,952 bytes |

The wheel is 57.6% smaller. These local wheels target manylinux 2.35; their
absolute sizes are not directly comparable to the published manylinux 2.17
builds. Remaining growth over v0.9.0 is not isolated by this comparison.

`docs/` was already absent from binary release archives, whose workflow adds
only the binary, README, and license. The wheel also includes no docs and
exactly one extension, despite stale development extensions in the checkout.
The Cargo manifest previously excluded only benchmark HTML; it now excludes
all `docs/**`, as does maturin's wheel/sdist configuration. Archive inspection
confirmed no `docs/` entries in the crate, sdist, or wheel. Test structures
remain in source packages because packaged tests and doctests use them.

## Code audit coverage and decisions

The delegated audit read all production Rust and Python code, Rust/Python
tests, and release/test workflows. The scientific-kernel and RMSD/clustering
subtasks covered 13,152 initial lines; the coordinating reviewer covered the
remaining interfaces, helpers, and workflows. Retained code was evaluated
against current callers and supported behavior, without adding explanatory
comments to every line.

| Area | Removed or simplified | Why the remaining code is needed |
| --- | --- | --- |
| `src/contacts/` | Repeated hydrogen filtering, temporary contact vectors, duplicate match arms, eager fallback atom extraction, unused SVD output | Chemical evidence, protonation policy, chain/model boundaries, and geometry determine supported contact classifications. |
| `src/sc/` | Custom RTree point wrapper, redundant vector clearing, manual trailing-space trimming and wildcard prefix matching | Surface construction, trimming, directional scoring, and calibrated radii define shape complementarity. Native RTree points and `GeomWithData` provide the removed wrapper's behavior. |
| `src/sap.rs`, `src/sasa.rs` | Unnecessary chain collection and repeated SAP sorting; corrected SAP aggregation documentation | SASA populations, polarity, exposure definitions, calibration tables, and preparation warnings implement distinct scientific contracts. |
| `src/rmsd.rs`, `src/pairwise_rmsd.rs`, `src/clustering.rs` | Recomputed residual centroid; replaced lazy NDJSON reader | Exact correspondence, numerical conditioning, degeneracy checks, matrix validation, memory limits, tie handling, and convergence checks are exercised supported behavior. |
| `src/metadata.rs`, `src/structure.rs`, `src/sequences.rs` | No speculative rewrite | PDB/mmCIF metadata, explicit bonds and chain breaks, conformer selection, and observed/declared sequence distinctions require their separate paths. |
| `src/python.rs`, `python/arpeggia/` | Duplicate scalar dSASA calculation and unused contract defaults/schema constants; stale contacts argument documentation | Public signatures, type aliases, consumed schema contracts, exception mapping, thread handling, and warning propagation serve Python callers. |
| `src/main.rs`, `src/cli/`, `src/lib.rs`, `src/utils.rs`, `src/diagnostics.rs` | Repeated writer error conversion and commented-out export | Commands, public exports, typed errors, diagnostics, output validation, and no-clobber/alias checks protect supported interfaces. |
| Tests and `.github/workflows/` | Duplicate contacts test; repaired signed angle tolerances and removed an incorrect angle assertion | Scientific regressions, failure paths, wheel smoke tests, and the platform/interpreter matrix cover real behavior. |

No scientific classification policy, calibration table, numerical safeguard,
or public API was removed. The geometry test had accepted negative errors of
arbitrary magnitude; it now uses absolute tolerances. Its removed assertion
incorrectly expected 90 degrees for geometry already checked as 45 degrees.

## Documentation

README features now map Python APIs to CLI commands in a compact table.
Chain Groups Specification is under Usage; Scientific conventions follows
Development. The table distinguishes scalar/sequence results from DataFrames.

CONTEXT.md no longer defines the already-completed v0.9 release gate.
Selected Conformer now refers to the concrete alternate-conformer warning,
avoiding confusion with the glossary's unresolved scientific-review warning.
Related terms such as declared/observed sequences, explicit/potential bonds,
and numerical/method compatibility remain separate because they describe
different behavior rather than duplicated statements.

## Validation

- 168 Rust library tests, 3 CLI unit tests, 13 integration tests, and 8 doctests passed.
- 16 Python tests passed against the newly built, separately installed wheel.
- Installed-wheel sequence smoke test passed.
- `ty check python` and all applicable `prek` hooks passed, including all-target/all-feature Clippy.
- Inspected actual wheel, sdist, and crate contents; docs are excluded and the wheel contains one extension.
- Added focused regressions for bounded/projected NDJSON reading, dSASA scalar/component agreement, and radius wildcard matching.

The pre-existing `uv.lock` working-tree edit was preserved.
