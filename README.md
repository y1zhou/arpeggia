# Arpeggia

This is a port of the [Arpeggio](https://github.com/PDBeurope/arpeggio/) library to Rust, with a focus on identifying certain protein-protein interactions in PDB and mmCIF files.

# Features

- [x] Parse PDB and mmCIF files
- [x] Parse user selection of chain groups
- [x] Extract protein chains and residues
- [x] Calculate distances between residues
- Identify protein-protein interactions
  - [x] Steric clashes
  - [x] VdW interactions
  - [x] Hydrophobic interactions
  - [x] Aromatic interactions
  - [x] Cation-pi interactions
  - [x] Ionic interactions
  - [x] Hydrogen bonds
  - [x] Weak hydrogen bonds
  - [x] Disulfide bonds
  - [x] Covalent bonds
- [x] Output results in various formats (e.g., JSON, CSV)
- [ ] Bundle into a `PyO3` extension module
- [ ] ~~Add hydrogens to the structure model~~

# Installation

## Cargo install

> TODO: publish to crates.io

## Pre-built binaries

Binaries for tagged releases can be found at the [GitHub Releases page](https://github.com/y1zhou/arpeggia/releases/).

## Build from source

1. Clone the repository:

    ```bash
    git clone https://github.com/y1zhou/arpeggia.git
    cd arpeggia
    ```

2. Build and install using Cargo:

    ```bash
    cargo install --path .
    ```

This will install the `arpeggia` binary to your Cargo binary directory (usually `~/.cargo/bin`).

# Usage

To see all available commands and options:

```bash
arpeggia help
```

> TODO: document how the `contacts`, `sasa`, and `seq` commands can be used.
