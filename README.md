# Arpeggia

This is a port of the [Arpeggio](https://github.com/PDBeurope/arpeggio/) library to Rust, with a focus on identifying certain protein-protein interactions in PDB and mmCIF files.

[![PyPI version](https://img.shields.io/pypi/v/arpeggia)](https://pypi.org/project/arpeggia/)
![License](https://img.shields.io/badge/license-GPL--3.0-blue.svg)
![Python versions](https://img.shields.io/pypi/pyversions/arpeggia)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/y1zhou/arpeggia)

## Features

| Python API | CLI command | Purpose |
| --- | --- | --- |
| `contacts()` | `contacts` | Atomic and aromatic contacts. |
| `sasa()` | `sasa` | Solvent accessible area per atom, residue, or chain. |
| `relative_sasa()` | `relative-sasa` | Residue SASA normalized by reference maximum areas. |
| `sap_score()` | `sap` | Spatial Aggregation Propensity per atom or residue. |
| `dsasa()`, `dsasa_components()` | `dsasa` | Two-sided buried interface area and polarity components. |
| `sc()` | `sc` | Shape complementarity between chain groups. |
| `seq()` | `seq` | Coordinate-observed protein sequences. |
| `seqres()` | `seqres` | Declared sequences, including residues without coordinates. |
| `align_seqs()` | `align-seqs` | Gapped pairwise alignment, identity, score, and edit distance. |
| `number_antibody()` | `number-antibody` | Antibody numbering, CDRs, germline similarities, and terminal imputation. |
| `align_antibodies()` | `align-antibodies` | Compare antibodies at shared numbered positions. |
| `rmsd()` | `rmsd` | Structural fit/evaluation with optional sequence alignment and rejection. |
| `pairwise_rmsd()` | `cluster-structs --pairwise-rmsd`¹ | Pairwise RMSD table for exactly corresponding structures. |
| `cluster_structs()` | `cluster-structs` | K-medoids clustering and representative structures. |

¹ The CLI writes the pair table as part of clustering; Python can calculate it independently.

See [antibody numbering](https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md),
[sequence alignment](https://github.com/y1zhou/arpeggia/blob/master/docs/sequence-alignment.md),
[structure comparison](https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md), and
[scientific conventions](https://github.com/y1zhou/arpeggia/blob/master/docs/scientific-conventions.md) for usage and assumptions.

Analyses accept PDB and mmCIF files with chain selections. Tabular Python
results are Polars DataFrames; CLI tables support CSV, Parquet, and NDJSON.

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
uv run --extra dev maturin develop --uv --release --features python --locked
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

# Superpose on chain A, measure chains B/C, and cluster a pairwise matrix
result = arpeggia.rmsd(
    "reference.cif",
    "query.cif",
    superpose_residues="A",
    rmsd_residues="B,C",
    atoms="ca",
)
print(result.rmsd, result.core_rmsd)
alignment = arpeggia.align_seqs("ACDEFGHIK", "ACDEYGHIK")
pairs = arpeggia.pairwise_rmsd("structures/", num_threads=8)
clusters = arpeggia.cluster_structs(pairwise_rmsd=pairs, max_clusters=10)
```

Tabular results are [Polars](https://pola.rs/) DataFrames for efficient data manipulation. You can easily convert to pandas if needed:

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

# RMSD and fixed-count clustering; save the reusable long pair table
arpeggia rmsd reference.cif query.cif --atoms backbone
arpeggia cluster-structs -i structures/ -o results/ \
  --num-clusters 5 --pairwise-rmsd --num-threads 8
```

To see all available options:

```bash
arpeggia help
arpeggia contacts --help
```

### Chain Groups Specification

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
uv run --extra dev maturin develop --uv --features python --locked
uv run --extra dev pytest
```

Rerun `maturin develop` after changing Rust code or switching branches with
native API changes. Editable installs expose Python edits immediately but retain
the compiled Rust extension until rebuilt.

To run Rust tests:

```bash
cargo test --locked
```

## License

GNU General Public License v3.0 - see [LICENSE](https://github.com/y1zhou/arpeggia/blob/master/LICENSE) for details.

## Credit

This project would not be possible without the following resources:

- [Arpeggio](https://github.com/PDBeurope/arpeggio/): Original Python library for protein-protein interaction analysis.
- [pdbtbx](https://github.com/douweschulte/pdbtbx/): The structural file parser doing all the heavy lifting.
- [RustSASA](https://github.com/maxall41/RustSASA): Library for calculating solvent accessible surface area.
- [sc-rs](https://github.com/cytokineking/sc-rs/): Library for calculating the Shape Complementarity by Lawrence & Colman (1993).
- [Rosetta](https://github.com/RosettaCommons/rosetta): Where the Spatial Aggregation Propensity (SAP) score calculations are inspired from.
