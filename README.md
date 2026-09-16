# Arpeggia

Arpeggia provides protein contact analysis, surface measurements, structural
comparison, sequence alignment and antibody numbering through Rust, Python and a
CLI. Contact analysis is based on [Arpeggio](https://github.com/PDBeurope/arpeggio/).

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

Structure analyses accept PDB and mmCIF files. Tabular Python results are
[Polars](https://pola.rs/) DataFrames; CLI tables support CSV, Parquet and NDJSON.
Sequence commands accept unaligned amino-acid strings.

## Installation

```bash
pip install arpeggia
```

Published wheels support CPython 3.10–3.14 on x86-64 Linux and Windows, and
Arm64 Linux and macOS. Download the CLI from
[GitHub Releases](https://github.com/y1zhou/arpeggia/releases/), or see
[BUILD.md](https://github.com/y1zhou/arpeggia/blob/master/BUILD.md) to build either interface from source.

## Usage

See [Python and CLI examples](https://github.com/y1zhou/arpeggia/blob/master/docs/examples.md)
for chain-group syntax, surface options, table output and reusable RMSD calculations.

### Python API

```python
import arpeggia

contacts = arpeggia.contacts("structure.pdb", groups="A/B")
contacts.write_parquet("contacts.parquet")
residue_sasa = arpeggia.sasa("structure.pdb", level="residue")
print(residue_sasa.sort("sasa").head(10))
relative_sasa = arpeggia.relative_sasa("structure.pdb")
sap = arpeggia.sap_score("structure.pdb", level="residue")
print(sap.sort("sap_score", descending=True).head(10))
total, polar, hydrophobic, unknown = arpeggia.dsasa_components(
    "structure.pdb", groups="A/B"
)
sc_score = arpeggia.sc("antibody_antigen.pdb", groups="H,L/A")
print(f"Shape Complementarity: {sc_score:.3f}")

print(arpeggia.align_seqs("ACDEFGHIK", "ACDEYGHIK"))
result = arpeggia.rmsd(
    "reference.cif", "query.cif",
    superpose_residues="A", rmsd_residues="B,C", atoms="ca",
)
print(result.rmsd, result.core_rmsd)
```

Use `help(arpeggia.contacts)` for arguments and defaults. The
[contact-table examples](https://github.com/y1zhou/arpeggia/blob/master/docs/scientific-conventions.md#contact-table-examples)
cover hydrogen-bond counts and interface residues; the feature guides above
cover alignment and numbering.

### Command-Line Interface

```bash
arpeggia contacts -i structure.pdb -o results/ -g "A/B" -t parquet
arpeggia sasa -i structure.pdb -o results/ --level residue --chains A,B
arpeggia dsasa -i structure.pdb -g "A/B"
arpeggia seq structure.pdb
arpeggia align-seqs ACDEFGHIK ACDEYGHIK
arpeggia rmsd reference.cif query.cif --atoms backbone
arpeggia cluster-structs -i structures/ -o results/ \
  --num-clusters 5 --pairwise-rmsd --num-threads 8
```

Run `arpeggia --help` to list commands, or `arpeggia <command> --help` for options.

### Chain Groups Specification

The `groups` parameter allows you to specify which chains interact with each other:

- `"/"` - All chains interact with all chains (including self) for contacts
- `"A,B/C,D"` - Chains A,B interact with chains C,D
- `"A/"` - Chain A interacts with all other chains
- `"A,B/"` - Chains A,B interact with all remaining chains

dSASA and SC require two disjoint, non-empty groups; `"/"` is therefore
invalid for those calculations.

## Development

[BUILD.md](https://github.com/y1zhou/arpeggia/blob/master/BUILD.md) contains locked
build, test and native-extension rebuild instructions.

## License

GNU General Public License v3.0 - see [LICENSE](https://github.com/y1zhou/arpeggia/blob/master/LICENSE) for details.

Bundled germline data are © 1995–2026 IMGT®, the international ImMunoGeneTics
information system®, Montpellier, France, licensed separately under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Arpeggia filters GENE-DB release 202636-7 and reformats llama protein displays;
see the [source credits, notices and preparation record](https://github.com/y1zhou/arpeggia/blob/master/data/germlines/README.md).

## Credit

- [Arpeggio](https://github.com/PDBeurope/arpeggio/): Original Python library for protein-protein interaction analysis.
- [pdbtbx](https://github.com/douweschulte/pdbtbx/): PDB and mmCIF parsing.
- [RustSASA](https://github.com/maxall41/RustSASA): Library for calculating solvent accessible surface area.
- [sc-rs](https://github.com/cytokineking/sc-rs/): Library for calculating the Shape Complementarity by Lawrence & Colman (1993).
- [Rosetta](https://github.com/RosettaCommons/rosetta): SAP reference definition and SASA polarity conventions.
