# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `number_antibody()` / `number-antibody` provide antibody variable-domain numbering with IMGT, Martin, AHo and Kabat,
  explicit CDR conventions, residue/input correspondence, and guards for weak
  matches, severe truncations, multiple domains and unsupported insertions.
  Partial domains must span the IMGT 23–118 framework-anchor interval to retain
  context for alignment, numbering and CDR conversion.
  Inputs outside 30–10,000 residues fail before sequence encoding.
  Python position labels compare and hash by number and insertion code for
  residue/column lookup; CLI diagnostics escape control characters in input names.
- Offline human, mouse, alpaca, rat and rabbit V/J germline similarities, with explicit
  species restrictions, known-residue coverage, tied gene/allele references,
  and attributed IMGT release 202636-7 data. `match_germlines=False` /
  `--no-germlines` skip matching without changing numbering or CDRs; matching
  remains enabled by default.
  `germlines_searched` distinguishes skipped searches from absent matches;
  skipped results display only input annotations and require matching before
  imputation. Binary release archives include
  the project README and a distinct `IMGT-GERMLINES.md` attribution file.
- Explicit terminal FR1/FR4 imputation returns a new numbered antibody and
  preserves original input indices and reference provenance. Tied references
  must agree, with optional explicit reference selection; internal gaps, CDRs
  and unknown input residues remain unchanged. Uncovered reference endpoints
  produce explicit coverage diagnostics.
- `align_antibodies()` / `align-antibodies` return `AntibodyAlignment` objects
  for one or more numbered antibodies through the ordered
  union of positions, preserving input row order and per-row CDR definitions.
  Heavy/light mixing and incompatible schemes fail; K/L mixtures are supported.
- Antibody displays place the input above combined V/J references, with separate
  gene labels, gray uncovered junctions and blank outer padding. Per-sequence
  rulers preserve original input coordinates through imputation; reference-defined
  CDR bands span the marker, ruler and sequence rows, and yellow backgrounds identify
  imputed residues and their count. Germline rulers count continuously through
  V then J, with endpoint and gene labels beside the visible sequence.
  Summaries provide colored CDR legends and tied V/J names.
  Multi-antibody views show the selected reference first and hide
  germlines, with wrapping and color/ruler controls.
- Pairwise protein sequence alignment through Rust/Python `align_seqs` and CLI
  `align-seqs`: global, local, and query-full semi-global modes with BLOSUM62
  and configurable affine gap costs with hundredth precision and a 0.01 minimum.
  `SeqAlignment` includes gapped strings,
  operations, input spans, identity/coverage ratios, gap statistics, and full-input
  edit distance.
- Shared CLI/Python alignment displays with blue `:` markers for positive-score
  substitutions, colored edits and clipped tails,
  custom sequence names, terminal-width wrapping, position rulers, and
  width/color/ruler controls.
- Optional observed-sequence correspondence before two-structure RMSD, with
  reference-based selections, explicit chain maps or unique maximum-score
  inference, and diagnostics for omitted atoms.
  Indistinguishable reference chains require explicit mapping before pairwise scoring.
- Optional atom-wise rejection and refitting. Final evaluation retains all mapped
  selected pairs, including rejected fitting pairs. Unchanged fitting sets retain
  their core RMSD; degenerate surviving fits fail.

### Changed

- Grouped RMSD selection, fitting, correspondence, and pairwise calculations
  under one module, and sequence alignment with its matrix data and display.
  Existing crate-root exports and inline unit-test organization are preserved.
- Python `rmsd()` returns a read-only `RmsdResult` instead of a scalar; use `.rmsd`
  for full evaluation RMSD and `.core_rmsd` for retained fitting pairs. Rust
  `get_rmsd` accepts `RmsdOptions` and returns `Analysis<RmsdResult>`. CLI RMSD
  output reports details by default and supports structured `--json` output.
- Python sequence and antibody displays share `__repr__` through Python's
  standard `str()` fallback, with `.format()` for explicit display controls.
- Python structural and antibody diagnostics share warning emission, preserving
  `UserWarning` filters and reporting only newly added imputation diagnostics.
- Expanded CLI help and Google-style Python/IDE docstrings with reference
  selection syntax, independent defaults, units, result semantics, and examples.
  Probe-radius guidance explains crevice access and how SASA changes depend on
  the structure. Consolidated structure-comparison usage in a dedicated guide,
  linked repository files through GitHub for installed-package users, and
  clarified Python/CLI names in README.

### Validation

- Qualified Hyalite against 2,745 Biopython reference cases and checked structural
  correspondence, refinement, and terminal displays. See the
  [alignment validation and benchmarks](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/sequence-alignment.md).

## [0.9.2] - 2026-09-08

### Changed

- Simplified chain-group splitting and duplicate observation-ID checks while
  preserving accepted grammar, path uniqueness checks, and validation order.
  Documented the already-sorted input expected by SAP DataFrame construction.
- Shared strong/weak hydrogen-bond geometry and aromatic plane fitting while
  retaining donor rules, classification thresholds, and missing-geometry warnings.
- Removed the ring-ring candidate-index vector by classifying borrowed pairs
  within the parallel ring loop. Atom and ring preparation remain outside it;
  atom-atom and ring-atom candidate handling is unchanged.
- Replaced the Polars lazy NDJSON scanner with eager reading and schema
  projection, retaining bounded cache row-count validation. Removing the lazy
  engine eliminates 115 normal dependency nodes and reduced a same-machine
  release wheel from 22.38 MB to 9.48 MB (57.6%). These are local build results,
  not measurements of published release artifacts.
- Excluded `docs/` from Rust crates and Python source distributions and wheels.
- Removed redundant contact geometry work, reused native spatial-index and
  string-matching operations in SC, and avoided repeated RMSD and SAP work.
- Shared Python dSASA and CLI output handling, and removed unused internal
  Python contract declarations while retaining supported APIs and schemas.
- Consolidated development documentation into ADRs, retained additional research
  and benchmark evidence, and reorganized README features and usage guidance.
  Removed machine-specific paths, unavailable benchmark commands, and repeated prose.
  Scientific conventions now include a contact-identification decision diagram
  in a dedicated [document](https://github.com/y1zhou/arpeggia/blob/master/docs/scientific-conventions.md).

### Fixed

- Replaced the ineffective zero-occupancy test with a fixture that verifies
  exactly which contacts disappear and preserves unaffected rows. Replaced
  stub source-text assertions with typed public-API usage checked by `ty`.
- Propagated SC surface-sampling errors instead of silently omitting failed
  patches. Legitimately empty geometry remains valid. Removed redundant
  first-atom guards with comments documenting their caller-established
  preconditions, and simplified an infallible probe helper.
- Corrected geometry-test tolerances to check absolute errors and removed an
  inconsistent angle assertion.
- Made Python tests explicitly assert expected missing-hydrogen, unresolved
  histidine, and incomplete-geometry warnings from their fixtures.

### Validation

- Isolated ring-ring benchmarks on 5B8C and 6BFT were faster with one and eight
  threads and retained identical outputs. A synthetic 5,632-ring case reduced
  peak process RSS from about 249 MiB to 7 MiB. These are ring-classification
  measurements, not end-to-end contact-analysis speedups; details are in the
  [cleanup audit](https://github.com/y1zhou/arpeggia/blob/master/docs/research/v0.9.2-cleanup-audit.md#ring-ring-enumeration-benchmark).
- The [5B8C benchmark and cleanup audit](https://github.com/y1zhou/arpeggia/blob/master/docs/research/v0.9.2-cleanup-audit.md)
  records essentially unchanged contact generation (+0.8–1.2%) and unchanged
  or faster NDJSON loading on a 2,574-row contacts table after removing the
  lazy engine. These measurements do not establish performance for all inputs.

## [0.9.1] - 2026-09-02

### Added

- Added uniform-weight Kabsch RMSD for exact protein-atom correspondence through
  Rust, Python, and the `rmsd` CLI command, with chain/residue and atom-subset
  selection. Superposition and RMSD residue selections can differ, allowing a
  structure to be aligned on one region while motion is measured in another.
- Added deterministic fixed-count FasterPAM and bounded automatic DynMSC
  structure clustering, packed parallel pairwise RMSD, heuristic memory guards,
  Polars outputs, and ID-only CLI pair-table reuse.

### Changed

- Simplified tabular IO to CSV, Parquet, and NDJSON. NDJSON reads use Polars'
  bounded lazy scanner with projection pushdown; ordinary JSON tables are no
  longer accepted or produced.
- Reduced the Rust public surface to supported scientific APIs. Contact and SC
  implementation details, chain-group parsing, and CLI-only table writers are
  now private; Python functions and CLI behavior are unchanged.

### Fixed

- Hardened Kabsch superposition and RMSD rescaling for degenerate, extreme, and
  subnormal coordinate cases.
- Strengthened pairwise RMSD cache, matrix, selection, and output validation.

## [0.9.0] - 2026-08-27

### Scientific corrections

- Restored both directional SC surfaces after fixing the inverted atom-2
  toroidal-surface branch and model-local plane traversal. Ringless structures
  are valid; absolute parity with the pinned `sc-rs` CLI remains explicitly
  documented in the regression test.
- Split resolved `Disulfide` and `Covalent` evidence from geometry-inferred
  `PotentialDisulfide` and distance-inferred `PotentialCovalent`, and named
  steric clash, van der Waals clash, and van der Waals contact regions
  separately.
- Associated hydrogen geometry with the specific donor rather than every
  hydrogen in its residue, restricted weak donors to carbon sites bearing
  hydrogen, restored explicitly protonated terminal-Pro donation, and added
  missing-hydrogen diagnostics.
- Added `AllCharged`, `Heuristic`, and `ExplicitOnly` histidine policies with
  potential ionic and cation-pi categories for inferred charge, including
  consistent treatment of histidine aliases and contradictory-input warnings.
- Preserved model, alternate-location, residue-name, and mmCIF label identity
  when applying explicit connectivity to selected atoms.
- Applied deterministic highest-occupancy alternate-conformer selection with an
  `A` tie-break and visible warnings across analyses.
- Unified standard atom, residue, and chain SASA over one selected atom
  population with ProtOr/fallback radii. Added Rosetta `SasaFilter` atom
  polarity and additive polar, hydrophobic, and unclassified areas.
- Aligned SAP with Rosetta's full-atom Reduce-radius exposure definition,
  1.1 Å probe, precise hydrophobicities, and runtime maximum side-chain areas.
  Exposure ratios remain unclamped; residue output retains complete side-chain
  SASA and eligible zero/nonpositive-score residues.
- Kept dSASA as the established two-sided buried area, corrected its docs,
  required disjoint groups, and added polarity components.
- Fixed coordinate-observed `seq` output to omit solvent and added `seqres` for
  PDB `SEQRES` and mmCIF entity-polymer declarations.

### API and schema changes

- Rust parsing and analysis entry points now return typed errors and successful
  `Analysis<T>` values with stable warning codes instead of panicking.
- Python emits `UserWarning`, `ValueError`, `RuntimeError`, or `OSError` while
  retaining direct DataFrame/scalar/list results. CPU-bound parsing and
  calculations release the GIL.
- Contact interaction labels and SASA/RSA columns changed as described above.
  Model/atom identifiers in contact and surface DataFrames are unsigned 32-bit
  integers and residue identifiers are signed 32-bit integers. Python adds
  `dsasa_components()` and `seqres()`.
- Removed redundant Rust convenience wrappers in favor of the canonical
  `analyze_contacts`, `get_dsasa_components`, and `get_sc_details` entry points.
  Python function names, return values, and CLI command names are unchanged.

### Performance and security

- Prepare dSASA structures once, avoid debug-only DataFrame work unless debug
  logging is enabled, and remove an arbitrary SC sampling ceiling.
- Restrict output filenames to a single normal path component.
- Pin release-critical GitHub Actions to immutable commits, track Cargo and uv
  lockfiles, use locked/frozen CI and release commands, and smoke-test every
  built wheel before upload.
- Lock Python artifacts to official PyPI and publish only the CPython/platform
  combinations supported by the release dependencies.
- Updated locked `h2` metadata to 0.4.16 for RUSTSEC-2026-0258. The optional
  dependency is not compiled by Arpeggia, but the release lock remains clean.
- Stream the narrow mmCIF metadata parser and discard unrelated loop values
  without allocating them.

### Simplified

- Replaced dynamic SC radii with static data, deleted duplicate geometry and
  sequence abstractions, removed redundant SAP maps/wrappers and unused Cargo
  features, removed write-only SC surface state, centralized CLI input
  diagnostics, and replaced stale build notes.

## [0.8.1] - 2026-06-23

### Changed

- Changed `pdb2seq` output from `dict` to `list[tuple]` to preserve chain order
- Updated dependency versions (clap, nalgebra, polars, pyo3, pyo3-polars, rayon, tracing, and tracing-subscriber)
- Fixed issue where `model_num=0` did not default to the first model in the structure
- Refactored Rust library for improved performance and robustness

## [0.8.0] - 2026-02-05

### Added

- **Shape Complementarity (SC) score calculation**: New `sc` module for computing geometric complementarity at protein-protein interfaces following the Lawrence & Colman (1993) algorithm
  - Full Connolly molecular surface generation with contact, toroidal, and concave surface dots
  - RTree-based spatial indexing for efficient peripheral band trimming and neighbor search
  - Attention classification system for dot filtering
  - New `sc` CLI command: `arpeggia sc -i structure.pdb -g "H,L/A"`
  - Python binding: `arpeggia.sc()` function returning SC score as float
  - Multi-threading support with `--threads` parameter
  - Model selection with `--model` parameter

### Changed

- Performance optimizations in SC module using RTree spatial indexing (6x speedup)
- SC calculation uses pdbtbx's VdW radii as fallback when embedded radii table has no match

## [0.7.0] - 2026-01-30

### Added

- **Spatial Aggregation Propensity (SAP) score calculation**: New functionality to predict aggregation-prone regions in proteins based on the Chennamsetty et al. paper "Developability Index: A Rapid In Silico Tool for the Screening of Antibody Aggregation Propensity" (DOI: 10.1002/jps.22758)
  - New `sap.rs` module with Black & Mould (1991) hydrophobicity scale, glycine-normalized
  - `get_per_atom_sap_score()` for atom-level SAP calculation with R-tree spatial indexing
  - `get_per_residue_sap_score()` for residue-level SAP aggregation
  - New `sap` CLI command with `--level` flag for atom/residue level calculation
  - Python binding: `sap_score()` function with `level`, `sap_radius`, and `chains` parameters

- **Chain filtering for SASA and SAP functions**: New `chains` parameter to filter structures to specific chains before calculation
  - Added `--chains` / `-c` flag to `sasa`, `relative-sasa`, and `sap` CLI commands
  - Python bindings: `chains` parameter added to `sasa()`, `relative_sasa()`, and `sap_score()`
  - Empty string (default) keeps all chains; comma-separated chain IDs filter to specified chains (e.g., "H,L")
  - Example: `arpeggia sap -i antibody.pdb -o output/ -c "H,L"` to analyze only heavy and light chains

### Changed

- `prepare_pdb_for_sasa()` now accepts a `chains` parameter for chain filtering
- All SASA functions (`get_atom_sasa`, `get_residue_sasa`, `get_chain_sasa`, `get_relative_sasa`) updated with `chains` parameter
- Updated documentation with SAP score examples and chain filtering usage

## [0.6.0] - 2026-01-26

### Added

- **rust-sasa v0.9.2 upgrade**: Updated from v0.3.2 with API changes, performance improvements, and insertion code support
- New `--level` option for `sasa` CLI command to calculate SASA at different granularities:
  - `atom` (default): Per-atom SASA values
  - `residue`: Aggregated SASA by residue with `is_polar` classification
  - `chain`: Aggregated SASA by chain
- New `dsasa` CLI command to calculate buried surface area at the interface between chain groups
- New `relative-sasa` CLI command to calculate relative solvent accessible surface area (RSA) normalized by Tien et al. (2013) MaxASA values
- New library functions: `get_residue_sasa`, `get_chain_sasa`, `get_dsasa`, `get_relative_sasa`, `get_max_asa`
- Preprocessing to remove solvent, ions, and hydrogens before SASA calculations (`prepare_pdb_for_sasa`)
- Utility functions: `get_num_threads` for thread management, `sum_sasa` for SASA aggregation
- Python bindings:
  - `sasa()` now accepts a `level` parameter ("atom", "residue", "chain")
  - New `dsasa()` function for buried surface area calculation
  - New `relative_sasa()` function for RSA calculation
- Added `num_threads` parameter to all SASA functions for parallel processing control
- 16 new tests for SASA functionality

### Changed

- Refactored `lib.rs` into separate modules: `contacts.rs`, `sasa.rs`, `sequences.rs`
- SASA functions now use rust-sasa's new `SASAOptions<T>` builder API
- Renamed `--name` to `--filename` (short: `-f`) in `sasa` and `contacts` CLI commands for consistency
- Consolidated duplicate code: Python `dsasa()` now uses `get_dsasa()` library function directly

### Fixed

- CLI short flag conflict between `--name` and `--num-points` (both used `-n`)

## [0.5.1] - 2026-01-20

### Added

- New `--ignore-zero-occupancy` flag for the `contacts` CLI command to filter out atoms with zero occupancy
- New `ignore_zero_occupancy` parameter for the Python `contacts()` function

## [0.5.0] - 2025-12-24

### Added

- PyO3 bindings for interaction detection, SASA calculation, and sequence extraction
- New `arpeggia` Python package for easy installation and usage
- GitHub Actions workflow for building and testing the Python package

## [0.4.2] - 2025-03-17

### Added

- The new `Plane` struct for better abstraction of sidechain centroids and normals
- `get_contacts` now returns a single DataFrame, with additional sidechain centroid distance and dihedral columns
- Salt bridges (when a hydrogen bond and an ionic bond are both present) are now correctly identified
- Tests for some IO and interaction detection functions

### Fixed

- Duplicated rows when chains appear on both sides of the `group` CLI argument
- Ignore planes when there are fewer than three atoms in the sidechain

### Changed

- Better logging messages and documentation of methods

## [0.4.1] - 2025-02-19

### Added

- Better error messages and flag documentations

### Fixed

- Create parent directories if the output directory does not exist
- Only interactions within the same model of the input file is considered
- Rows in the output of `contacts` are now sorted more naturally
- Use one thread by default, as using more rarely gives any gains in performance
- Use a distance cutoff of 6.5Å for searching neighbor atoms by default, as the previous 4.5Å could miss certain Pi interactions

## [0.4.0] - 2025-02-17

### Added

- Detection of repulsion between like charges
- Support of parquet, json, and ndjson output formats
- Added a `-name` flag in the CLI to rename the output file of `contacts` and `sasa` commands

### Fixed

- Checks for polar contacts were skipped when hydrogen bond criteria are not satisfied
- Better error messages when rings have missing atoms for finding the center and normal vector
- Nomenclature mix of residue insertion codes and alternative locations; the two are now stored under separate columns (`*_insertion` and `*_altloc`) in the output files

### Changed

- Distance cutoff for T-shaped Pi-stacking lowered from 6Å to 5Å
- Added hydrogen bond distance check to better differentiate Hbonds and polar contacts
- Performance/memory footprint improvement by switching from 64-bit numbers to 32-bit
- Logging is now less verbose

## [0.3.1]

### Fixed

- Wrong theta angle specification for pi-pi interactions
- Reduced unnecessary cloning of objects and strings

### Changed

- Defaults to searching for all intra- and inter-chain interactions when `-g '/'` is passed to `contacts`

## [0.3.0] - 2024-08-14

### Added

- `sasa` command to calculate the atom level SASA
- `seq` command to extract protein sequences from PDB files

### Fixed

- Only report chains that are part of the ligand or receptor for `contacts`

### Changed

- Moved previous top-level command to `contacts` sub-command
- Better path parsing

## [0.2.0] - 2024-08-08

### Added

- Separate CLI and core methods to prepare for future Python tooling
- Dump results to CSV file

### Fixed

- Use `pdbtbx` version that can deal with non-standard PDB rows

### Changed

- As a consequence of the `pdbtbx` update, only atomic coordinates are now parsed

## [0.1.0] - 2024-05-09

### Added

- Initial release
- Detection of common protein-protein interactions in a PDB or mmCIF file

[Unreleased]: https://github.com/y1zhou/arpeggia/compare/v0.9.2...HEAD
[0.9.2]: https://github.com/y1zhou/arpeggia/compare/v0.9.1...v0.9.2
[0.9.1]: https://github.com/y1zhou/arpeggia/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/y1zhou/arpeggia/compare/v0.8.1...v0.9.0
[0.8.1]: https://github.com/y1zhou/arpeggia/releases/tag/v0.8.1
[0.8.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.8.0
[0.7.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.7.0
[0.6.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.6.0
[0.5.1]: https://github.com/y1zhou/arpeggia/releases/tag/v0.5.1
[0.5.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.5.0
[0.4.2]: https://github.com/y1zhou/arpeggia/releases/tag/v0.4.2
[0.4.1]: https://github.com/y1zhou/arpeggia/releases/tag/v0.4.1
[0.4.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.4.0
[0.3.1]: https://github.com/y1zhou/arpeggia/releases/tag/v0.3.1
[0.3.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.3.0
[0.2.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.2.0
[0.1.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.1.0
