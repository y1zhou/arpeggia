# Quick Start Guide

## For Python Users

### Installation

```bash
# Clone the repository
git clone https://github.com/y1zhou/arpeggia.git
cd arpeggia

# Install the Python package
pip install maturin polars
maturin develop -v --release --features python
```

### Basic Usage

```python
import arpeggia
import polars as pl

# Analyze protein contacts
contacts_df = arpeggia.contacts(
    "test-data/1ubq.pdb",
    groups="/",           # All-to-all interactions
    vdw_comp=0.1,        # VdW compensation
    dist_cutoff=6.5      # Distance cutoff (Å)
)

print(f"Found {len(contacts_df)} contacts")
print(contacts_df.head())

# Calculate SASA
sasa_df = arpeggia.sasa(
    "test-data/1ubq.pdb",
    probe_radius=1.4,
    n_points=100
)

print(f"Calculated SASA for {len(sasa_df)} atoms")

# Extract sequences
sequences = arpeggia.pdb2seq("test-data/1ubq.pdb")
for chain_id, seq in sequences.items():
    print(f"Chain {chain_id}: {seq[:50]}...")  # First 50 residues
```

### Working with Results

```python
# Filter contacts by interaction type
hydrogen_bonds = contacts_df.filter(
    pl.col("interaction") == "HydrogenBond"
)

# Calculate statistics
interaction_counts = contacts_df.group_by("interaction").count()
print(interaction_counts)

# Convert to pandas if needed
import pandas as pd
contacts_pd = contacts_df.to_pandas()

# Save to file
contacts_df.write_csv("contacts.csv")
contacts_df.write_parquet("contacts.parquet")
```

### Chain Groups Specification

- `"/"` - All chains with all chains (self-interactions included)
- `"A,B/C,D"` - Chains A,B interact with chains C,D
- `"A/"` - Chain A interacts with all other chains  
- `"/C,D"` - All chains interact with chains C,D

## For CLI Users

### Installation

```bash
# Build and install
cargo install --path .
```

### Basic Usage

```bash
# Analyze contacts
arpeggia contacts \
    -i structure.pdb \
    -o output_dir \
    -g "A,B/C,D" \
    -t csv

# Calculate SASA
arpeggia sasa \
    -i structure.pdb \
    -o output_dir \
    -r 1.4 \
    -n 100

# Extract sequences
arpeggia seq structure.pdb
```

### Output Formats

Supported formats: `csv`, `parquet`, `json`, `ndjson`

```bash
arpeggia contacts -i input.pdb -o output/ -t parquet
```

## Examples

### Example 1: Find Hydrogen Bonds

```python
import arpeggia
import polars as pl

df = arpeggia.contacts("structure.pdb")

# Filter for hydrogen bonds
h_bonds = df.filter(
    pl.col("interaction").is_in([
        "HydrogenBond",
        "WeakHydrogenBond"
    ])
)

# Group by residue pairs
pairs = h_bonds.group_by([
    "from_resn", "from_resi",
    "to_resn", "to_resi"
]).agg(pl.count())

print(pairs)
```

### Example 2: Analyze Interface

```python
import arpeggia
import polars as pl

# Get contacts between chains A and B
df = arpeggia.contacts("structure.pdb", groups="A/B")

# Calculate interface residues
interface_residues = df.select([
    pl.col("from_resn"),
    pl.col("from_resi"),
    pl.col("from_chain")
]).unique()

print(f"Interface has {len(interface_residues)} residues")
```

### Example 3: SASA Analysis

```python
import arpeggia
import polars as pl

sasa = arpeggia.sasa("structure.pdb")

# Find buried residues (low SASA)
buried = sasa.filter(pl.col("sasa") < 10.0)
print(f"Found {len(buried)} buried atoms")

# Calculate average SASA per residue
avg_sasa = sasa.group_by([
    "chain", "resi", "resn"
]).agg(pl.col("sasa").mean())

print(avg_sasa)
```

## Troubleshooting

### Import Error

```python
ImportError: cannot import name 'arpeggia'
```

**Solution:** Make sure you've built the package with maturin:
```bash
maturin develop --release --features python
```

### Build Errors

**Missing Python.h:**
```bash
# Ubuntu/Debian
sudo apt-get install python3-dev

# macOS
brew install python
```

**Version conflicts:**
```bash
cargo clean
maturin develop --release --features python
```

## Next Steps

- Read the [full README](README.md) for more details
- Check [BUILD.md](BUILD.md) for advanced build options
- See the [test script](python/test_arpeggia.py) for more examples
- Report issues on [GitHub](https://github.com/y1zhou/arpeggia/issues)
