# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/y1zhou/arpeggia/compare/v0.4.2...HEAD
[0.4.2]: https://github.com/y1zhou/arpeggia/releases/tag/v0.4.2
[0.4.1]: https://github.com/y1zhou/arpeggia/releases/tag/v0.4.1
[0.4.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.4.0
[0.3.1]: https://github.com/y1zhou/arpeggia/releases/tag/v0.3.1
[0.3.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.3.0
[0.2.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.2.0
[0.1.0]: https://github.com/y1zhou/arpeggia/releases/tag/v0.1.0
