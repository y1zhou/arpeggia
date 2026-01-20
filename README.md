# Arpeggia

This is a port of the [Arpeggio](https://github.com/PDBeurope/arpeggio/) library to Rust, with a focus on identifying certain protein-protein interactions in PDB and mmCIF files.

![PyPI - Version](https://img.shields.io/pypi/v/arpeggia)

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
  - [x] Disulfide bonds
  - [x] Covalent bonds
- [x] Output results in various formats (e.g., JSON, CSV, Parquet)
- [x] Python bindings via PyO3
- [x] Returns Polars DataFrames for efficient data manipulation

## Installation

### Python Package (Recommended)

Install using pip:

```bash
pip install arpeggia
```

Or install from source using maturin:

```bash
git clone https://github.com/y1zhou/arpeggia.git
cd arpeggia
pip install maturin
maturin develop -v --release --features python
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
    ignore_zero_occupancy=False   # Filter out atoms with zero occupancy
)
print(f"Found {len(contacts_df)} contacts")
print(contacts_df.head())

# Calculate solvent accessible surface area
sasa_df = arpeggia.sasa(
    "structure.pdb",
    probe_radius=1.4,    # Probe radius in Ångströms
    n_points=100,        # Number of surface points
    model_num=0          # Model number (0 = first model)
)
print(f"Calculated SASA for {len(sasa_df)} atoms")

# Extract protein sequences
sequences = arpeggia.pdb2seq("structure.pdb")
for chain_id, seq in sequences.items():
    print(f"Chain {chain_id}: {seq}")
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

# Calculate SASA
arpeggia sasa -i structure.pdb -o output_dir -r 1.4 -n 100

# Extract sequences
arpeggia seq structure.pdb
```

To see all available options:

```bash
arpeggia help
arpeggia contacts --help
```

## Chain Groups Specification

The `groups` parameter allows you to specify which chains interact with each other:

- `"/"` - All chains interact with all chains (including self)
- `"A,B/C,D"` - Chains A,B interact with chains C,D
- `"A/"` - Chain A interacts with all other chains
- `"A,B/"` - Chains A,B interact with all remaining chains

## Development

To build the Python package in development mode:

```bash
pip install maturin polars
maturin develop -v --release --features python
python python/test_arpeggia.py
```

To run Rust tests:

```bash
cargo test
```

## License

MIT License - see LICENSE file for details.
