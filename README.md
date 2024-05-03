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
  - [ ] Aromatic interactions
  - [ ] Cation-pi interactions
  - [x] Ionic interactions
  - [x] Hydrogen bonds
  - [x] Weak hydrogen bonds
  - [ ] Disulfide bonds
  - [ ] Covalent bonds
- [ ] Output results in various formats (e.g., JSON, CSV)
- [ ] Bundle into a `PyO3` extension module
- [ ] Add hydrogens to the structure model
