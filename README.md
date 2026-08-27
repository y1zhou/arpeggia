# Arpeggia

This is a port of the [Arpeggio](https://github.com/PDBeurope/arpeggio/) library to Rust, with a focus on identifying certain protein-protein interactions in PDB and mmCIF files.

[![PyPI version](https://img.shields.io/pypi/v/arpeggia)](https://pypi.org/project/arpeggia/)
![License](https://img.shields.io/badge/license-GPL--3.0-blue.svg)
![Python versions](https://img.shields.io/pypi/pyversions/arpeggia)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/y1zhou/arpeggia)

## Features

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
  - [x] Resolved disulfide/covalent bonds and inferred potential contacts
- [x] Calculate SASA (Solvent Accessible Surface Area) at atom, residue, and chain levels
- [x] Calculate relative SASA (RSA) normalized by MaxASA values
- [x] Calculate SAP (Spatial Aggregation Propensity) scores for aggregation prediction
- [x] Calculate Shape Complementarity (SC) scores at protein-protein interfaces
- [x] Filter calculations to specific chains
- [x] Output results in various formats (e.g., JSON, CSV, Parquet)
- [x] Python bindings via PyO3
- [x] Returns Polars DataFrames for efficient data manipulation

## Scientific conventions

- Contact rows use `Disulfide` for resolved PDB `SSBOND` or mmCIF disulfide
  declarations and `Covalent` for other resolved `LINK`, `CONECT`, or
  `_struct_conn` bonds. Undeclared CYS pairs matching the original distance and
  CB--SG--SG--CB dihedral rule produce `PotentialDisulfide`; other contacts in
  the covalent-distance band produce `PotentialCovalent`. Clash and van der
  Waals regions are separately named.
- Explicit hydrogen-bond geometry uses only hydrogens associated with the donor
  atom. Missing donor hydrogens produce warnings; Arpeggia does not protonate
  input structures.
- Histidines use `AllCharged` by default for Arpeggio-compatible
  positive-ionisable typing. `Heuristic` applies explicit evidence followed by
  a pH-dependent intrinsic-pKa prior, while `ExplicitOnly` never guesses.
  Inferred histidine charge produces potential ionic, repulsion, and cation-pi
  labels rather than definitive ones.
- All analyses deterministically choose the highest-occupancy alternate
  conformer, with `A` as the tie-breaker, and warn when selection occurs.
- Standard atom, residue, and chain SASA use one atom population and ProtOr
  radii with elemental fallback. Polar/hydrophobic columns follow Rosetta's
  legacy `SasaFilter` atom partition; numerical areas remain Shrake–Rupley.
- dSASA is the two-sided buried area
  `SASA(group 1) + SASA(group 2) - SASA(complex)`. Divide by two only when a
  one-sided interface-area convention is required.
- SAP uses the Rosetta-compatible full-atom Reduce-radius exposure definition
  with a 1.1 Å default probe and sums positive score contributions while
  reporting complete side-chain SASA. Arpeggia does not add missing atoms, so
  direct Rosetta comparison requires the same caller-prepared full-atom input.
  Monomers without a Rosetta calibration are omitted with a warning.

## Installation

### Python Package (Recommended)

Install using pip:

```bash
pip install arpeggia
```

Published wheels support CPython 3.10–3.14 on x86-64 Linux and Windows, and
Arm64 Linux and macOS.

Or install from source using maturin:

```bash
git clone https://github.com/y1zhou/arpeggia.git
cd arpeggia
uv sync --frozen --all-extras
uv run maturin develop --uv --release --features python --locked
```

### Rust Binary

For the command-line tool, you can install pre-built binaries from the [GitHub Releases page](https://github.com/y1zhou/arpeggia/releases/), or build from source:

```bash
git clone https://github.com/y1zhou/arpeggia.git
cd arpeggia
cargo install --path .
```

This will install the `arpeggia` binary to your Cargo binary directory (usually `~/.cargo/bin`).

## Usage

### Python API

```python
import arpeggia

# Analyze protein contacts
contacts_df = arpeggia.contacts(
    "structure.pdb",
    groups="/",                    # All-to-all chain interactions
    vdw_comp=0.1,                 # VdW radii compensation
    dist_cutoff=6.5,              # Distance cutoff in Ångströms
    ignore_zero_occupancy=False,  # Set True to ignore zero occupancy atoms
    protonation="all-charged",    # Or "heuristic" / "explicit-only"
    ph=7.4,
)
print(f"Found {len(contacts_df)} contacts")
print(contacts_df.head())

# Calculate solvent accessible surface area
# Atom-level (default)
sasa_df = arpeggia.sasa("structure.pdb", level="atom", probe_radius=1.4, n_points=100, model_num=0)
print(f"Calculated SASA for {len(sasa_df)} atoms")

# Residue-level SASA
residue_sasa = arpeggia.sasa("structure.pdb", level="residue")
print(f"Calculated SASA for {len(residue_sasa)} residues")

# Chain-level SASA for specific chains only
chain_sasa = arpeggia.sasa("structure.pdb", level="chain", chains="A,B")
print(f"Calculated SASA for chains A and B")

# Calculate relative SASA (RSA) normalized by Tien et al. (2013) MaxASA values
rsa_df = arpeggia.relative_sasa("structure.pdb")
print(f"Calculated RSA for {len(rsa_df)} residues")

# Calculate Spatial Aggregation Propensity (SAP) scores for aggregation prediction
sap_df = arpeggia.sap_score("antibody.pdb", level="residue")
print(f"Calculated SAP for {len(sap_df)} residues")

# SAP for specific chains (e.g., antibody heavy and light chains)
sap_hl = arpeggia.sap_score("antibody.pdb", chains="H,L", sap_radius=5.0)
print(f"Calculated SAP for H and L chains")

# Calculate buried surface area at the interface
bsa = arpeggia.dsasa("structure.pdb", groups="A,B/C,D")
print(f"Buried surface area: {bsa:.2f} Å²")

# Additive two-sided polarity components
total, polar, hydrophobic, unknown = arpeggia.dsasa_components(
    "structure.pdb", groups="A,B/C,D"
)

# Calculate Shape Complementarity at an interface
sc_score = arpeggia.sc("antibody_antigen.pdb", groups="H,L/A")
print(f"Shape Complementarity: {sc_score:.3f}")  # Typical values: 0.5-0.7

# Extract protein sequences
sequences = arpeggia.seq("structure.pdb")
for chain_id, seq in sequences:
    print(f"Chain {chain_id}: {seq}")

# Extract declared SEQRES/entity-polymer sequences, including missing coordinates
declared_sequences = arpeggia.seqres("structure.pdb")
```

The functions return [Polars](https://pola.rs/) DataFrames for efficient data manipulation. You can easily convert to pandas if needed:

```python
import polars as pl

# Convert to pandas
contacts_pd = contacts_df.to_pandas()

# Or save directly to various formats
contacts_df.write_csv("contacts.csv")
contacts_df.write_parquet("contacts.parquet")
```

### Command-Line Interface

The CLI provides the same functionality:

```bash
# Analyze contacts
arpeggia contacts -i structure.pdb -o output_dir -g "A,B/C,D" -t csv

# Analyze contacts, ignoring atoms with zero occupancy
arpeggia contacts -i structure.pdb -o output_dir --ignore-zero-occupancy

# Calculate SASA at different levels (atom, residue, chain)
arpeggia sasa -i structure.pdb -o output_dir --level atom
arpeggia sasa -i structure.pdb -o output_dir --level residue
arpeggia sasa -i structure.pdb -o output_dir --level chain

# Calculate SASA for specific chains only
arpeggia sasa -i structure.pdb -o output_dir --level residue --chains "A,B"

# Calculate relative SASA (RSA) for each residue
arpeggia relative-sasa -i structure.pdb -o output_dir

# Calculate SAP scores for aggregation prediction
arpeggia sap -i antibody.pdb -o output_dir --level residue

# Calculate SAP for specific chains (e.g., antibody H and L chains)
arpeggia sap -i antibody.pdb -o output_dir --chains "H,L"

# Calculate buried surface area at the interface
arpeggia dsasa -i structure.pdb -g "A,B/C,D"

# Calculate Shape Complementarity at an interface
arpeggia sc -i antibody_antigen.pdb -g "H,L/A"

# Extract sequences
arpeggia seq structure.pdb

# Extract declared SEQRES/entity-polymer sequences
arpeggia seqres structure.pdb
```

To see all available options:

```bash
arpeggia help
arpeggia contacts --help
```

## Chain Groups Specification

The `groups` parameter allows you to specify which chains interact with each other:

- `"/"` - All chains interact with all chains (including self) for contacts
- `"A,B/C,D"` - Chains A,B interact with chains C,D
- `"A/"` - Chain A interacts with all other chains
- `"A,B/"` - Chains A,B interact with all remaining chains

dSASA and SC require two disjoint, non-empty groups; `"/"` is therefore
invalid for those calculations.

## Development

To build the Python package in development mode:

```bash
uv sync --frozen --all-extras
uv run maturin develop --uv --features python --locked
uv run pytest
```

To run Rust tests:

```bash
cargo test --locked
```

## License

GNU General Public License v3.0 - see [LICENSE](LICENSE) for details.

## Credit

This project would not be possible without the following resources:

- [Arpeggio](https://github.com/PDBeurope/arpeggio/): Original Python library for protein-protein interaction analysis.
- [pdbtbx](https://github.com/douweschulte/pdbtbx/): The structural file parser doing all the heavy lifting.
- [RustSASA](https://github.com/maxall41/RustSASA): Library for calculating solvent accessible surface area.
- [sc-rs](https://github.com/cytokineking/sc-rs/): Library for calculating the Shape Complementarity by Lawrence & Colman (1993).
- [Rosetta](https://github.com/RosettaCommons/rosetta): Where the Spatial Aggregation Propensity (SAP) score calculations are inspired from.
