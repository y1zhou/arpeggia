"""Test suite for arpeggia Python bindings using pytest."""

# ruff: noqa: S101
from pathlib import Path
from typing import Any, cast

import pytest


@pytest.fixture
def test_pdb_file():
    """Fixture providing path to test PDB file."""
    test_file = Path(__file__).parent.parent.parent / "test-data" / "1ubq.pdb"
    if not test_file.exists():
        pytest.skip(f"Test file not found: {test_file}")
    return str(test_file)


def test_import():
    """Test that the module can be imported."""
    import arpeggia
    from arpeggia import _contract

    assert hasattr(arpeggia, "__version__")
    assert hasattr(arpeggia, "contacts")
    assert hasattr(arpeggia, "sasa")
    assert hasattr(arpeggia, "seq")
    assert hasattr(arpeggia, "rmsd")
    assert hasattr(arpeggia, "cluster_structs")
    assert arpeggia.__all__ == list(_contract.EXPORTED_FUNCTIONS)


def test_stub_consumes_python_contract():
    """Test that the type stub imports the shared exposure contract."""
    stub = Path(__file__).parent.parent / "arpeggia" / "arpeggia.pyi"
    stub_text = stub.read_text()

    assert "ProtonationMode" in stub_text
    assert "level: SasaLevel = ..." in stub_text
    assert "level: SapLevel = ..." in stub_text
    assert "def seq(input_file: str, model_num: int = ...) -> SequenceList" in stub_text
    assert "atoms: AtomSubset = ..." in stub_text


def test_rmsd_pairwise_and_clustering(test_pdb_file, tmp_path):
    """Expose exact-correspondence RMSD and clustering through Python."""
    import arpeggia
    import polars as pl
    from arpeggia import _contract

    assert (
        arpeggia.rmsd(
            test_pdb_file,
            test_pdb_file,
            superpose_residues="A:1-20",
            rmsd_residues="A:1-20",
        )
        < 1e-12
    )
    assert (
        arpeggia.rmsd(
            test_pdb_file,
            test_pdb_file,
            superpose_residues="A:1-20",
        )
        < 1e-12
    )
    with pytest.raises(ValueError, match="RMSD Selection"):
        arpeggia.rmsd(
            "missing-reference.pdb",
            "missing-mobile.pdb",
            rmsd_residues="A:",
        )

    structures = tmp_path / "structures"
    structures.mkdir()
    for structure_id in ("a", "b", "c"):
        (structures / f"{structure_id}.pdb").write_bytes(
            Path(test_pdb_file).read_bytes()
        )
    pairs = arpeggia.pairwise_rmsd(
        str(structures),
        superpose_residues="A:1-20",
        rmsd_residues="A:1-20",
        num_threads=2,
    )
    assert tuple(pairs.columns) == _contract.PAIRWISE_RMSD_COLUMNS
    assert pairs.height == 3
    assert pairs["rmsd"].max() == pytest.approx(0.0, abs=1e-12)

    clusters = arpeggia.cluster_structs(pairwise_rmsd=pairs, num_clusters=1)
    assert tuple(clusters.columns) == _contract.CLUSTER_COLUMNS
    assert clusters.schema == {
        "id": pl.String,
        "cluster_id": pl.UInt32,
        "medoid_id": pl.String,
        "rmsd_to_medoid": pl.Float64,
    }
    duplicate_clusters = arpeggia.cluster_structs(pairwise_rmsd=pairs, num_clusters=2)
    assert set(duplicate_clusters["cluster_id"]) == {0, 1}
    with pytest.raises(ValueError, match="max_clusters must be between 2 and 2"):
        arpeggia.cluster_structs(pairwise_rmsd=pairs, max_clusters=3)
    with pytest.raises(ValueError, match="exactly one"):
        arpeggia.cluster_structs(
            input=str(structures), pairwise_rmsd=pairs, num_clusters=1
        )
    with pytest.raises(ValueError, match="exactly one"):
        arpeggia.cluster_structs(num_clusters=1)
    with pytest.raises(ValueError, match="requires column id_1"):
        arpeggia.cluster_structs(
            pairwise_rmsd=pl.DataFrame({"wrong": [1]}), num_clusters=1
        )
    with pytest.raises(ValueError, match="max_iterations must be positive"):
        arpeggia.cluster_structs(
            pairwise_rmsd=pl.DataFrame({"wrong": [1]}),
            num_clusters=1,
            max_iterations=0,
        )
    with pytest.raises(ValueError, match="atoms must be"):
        arpeggia.rmsd(
            "missing-reference.pdb",
            "missing-mobile.pdb",
            atoms=cast(Any, "invalid"),
        )


def test_contacts(test_pdb_file):
    """Test the contacts function returns expected DataFrame structure."""
    import arpeggia
    from arpeggia import _contract

    df = arpeggia.contacts(test_pdb_file, groups="/", vdw_comp=0.1, dist_cutoff=6.5)

    # Check DataFrame is not empty
    assert df.height > 0, "Contacts DataFrame should not be empty"

    # Check expected columns exist
    expected_columns = _contract.CONTACT_COLUMNS

    for col in expected_columns:
        assert col in df.columns, (
            f"Column '{col}' should be present in contacts DataFrame"
        )

    # Check shape - should have 20 columns (all expected columns)
    assert df.width == 20, f"Expected 20 columns, got {df.width}"

    # Verify some basic properties
    assert df["distance"].dtype.is_float(), "Distance column should be float type"
    assert all(df["distance"] >= 0), "All distances should be non-negative"


def test_contacts_chain_groups(test_pdb_file):
    """Test contacts with specific chain groups."""
    import arpeggia

    # Test with specific chain if available
    df = arpeggia.contacts(test_pdb_file, groups="/", vdw_comp=0.1, dist_cutoff=6.5)

    # Should have some interactions
    assert len(df) > 0


def test_contacts_ignore_zero_occupancy(test_pdb_file):
    """Test contacts with ignore_zero_occupancy parameter."""
    import arpeggia

    # Test with ignore_zero_occupancy=False (default)
    df1 = arpeggia.contacts(
        test_pdb_file,
        groups="/",
        vdw_comp=0.1,
        dist_cutoff=6.5,
        ignore_zero_occupancy=False,
    )

    # Test with ignore_zero_occupancy=True
    df2 = arpeggia.contacts(
        test_pdb_file,
        groups="/",
        vdw_comp=0.1,
        dist_cutoff=6.5,
        ignore_zero_occupancy=True,
    )

    # Both should return valid DataFrames
    assert len(df1) > 0
    assert len(df2) > 0

    # For 1ubq.pdb, all atoms have occupancy 1.0, so results should be identical
    assert len(df1) == len(df2)


def test_sasa(test_pdb_file):
    """Test the sasa function returns expected DataFrame structure."""
    import arpeggia
    from arpeggia import _contract

    df = arpeggia.sasa(test_pdb_file, probe_radius=1.4, n_points=100, model_num=0)

    # Check DataFrame is not empty
    assert df.height == 602, "SASA DataFrame should not be empty"

    # Check expected columns exist
    expected_columns = _contract.SASA_COLUMNS["atom"]

    for col in expected_columns:
        assert col in df.columns, f"Column '{col}' should be present in SASA DataFrame"

    assert df.shape[1] == 9

    # Verify SASA values are reasonable
    assert df["sasa"].dtype.is_float(), "SASA column should be float type"
    assert all(df["sasa"] >= 0), "All SASA values should be non-negative"
    assert any(df["sasa"] > 0), "At least some atoms should have non-zero SASA"


def test_sasa_parameters(test_pdb_file):
    """Test SASA with different parameters."""
    import arpeggia

    # Test with different probe radius
    df1 = arpeggia.sasa(test_pdb_file, probe_radius=1.4, n_points=100)
    df2 = arpeggia.sasa(test_pdb_file, probe_radius=2.0, n_points=100)

    # Both should return data
    assert len(df1) > 0
    assert len(df2) > 0

    # Different probe radius should give different SASA values
    # (though the number of atoms should be the same)
    assert len(df1) == len(df2)


def test_sasa_and_sap_default_model_uses_first_explicit_model(test_pdb_file, tmp_path):
    """Test model_num=0 with explicit model serials that do not start at 1."""
    import arpeggia

    source_lines = Path(test_pdb_file).read_text().splitlines()
    coordinate_lines = [
        line for line in source_lines if line.startswith(("ATOM  ", "HETATM", "TER   "))
    ]
    multimodel_pdb = tmp_path / "first-model-serial-7.pdb"
    multimodel_pdb.write_text(
        "MODEL        7\n"
        + "\n".join(coordinate_lines)
        + "\nENDMDL\n"
        + "MODEL        9\n"
        + "\n".join(coordinate_lines)
        + "\nENDMDL\nEND\n"
    )

    default_sasa = arpeggia.sasa(str(multimodel_pdb), model_num=0)
    explicit_sasa = arpeggia.sasa(str(multimodel_pdb), model_num=7)
    default_sap = arpeggia.sap_score(str(multimodel_pdb), model_num=0)
    explicit_sap = arpeggia.sap_score(str(multimodel_pdb), model_num=7)

    assert default_sasa.height == explicit_sasa.height == 602
    assert default_sap.height == explicit_sap.height > 0


def test_seq(test_pdb_file):
    """Test the seq function returns expected structure."""
    import arpeggia

    seqs = arpeggia.seq(test_pdb_file)

    # Check return type
    assert isinstance(seqs, list), "Sequences should return a list"
    assert len(seqs) > 0, "Should have at least one chain"

    # For 1ubq.pdb, we know it has 1 chain with a specific sequence
    # Chain should be present
    assert len(seqs) == 1, f"Expected 1 chain, got {len(seqs)}"

    # Check sequence properties
    for chain_id, seq in seqs:
        assert isinstance(chain_id, str), "Chain ID should be string"
        assert isinstance(seq, str), "Sequence should be string"
        assert len(seq) > 0, "Sequence should not be empty"

        # For 1ubq, the sequence should be 76 residues
        # This is the known ubiquitin sequence
        assert len(seq) == 76, f"Expected 76 residues for ubiquitin, got {len(seq)}"

        # Check it starts with the expected sequence
        expected_start = "MQIFVKTLTG"
        assert seq.startswith(expected_start), (
            f"Sequence should start with {expected_start}, got {seq[:10]}"
        )


def test_seq_selects_one_model(tmp_path):
    """Select observed sequence from one explicit model."""
    import arpeggia

    structure = tmp_path / "models.pdb"
    structure.write_text(
        "MODEL        7\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "ENDMDL\n"
        "MODEL        9\n"
        "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "ENDMDL\nEND\n"
    )
    assert arpeggia.seq(str(structure)) == [("A", "A")]
    assert arpeggia.seq(str(structure), model_num=9) == [("A", "G")]
    with pytest.raises(ValueError, match="model 3 does not exist"):
        arpeggia.seq(str(structure), model_num=3)


def test_sequences_validity(test_pdb_file):
    """Test that returned sequences contain valid amino acid codes."""
    import arpeggia

    seqs = arpeggia.seq(test_pdb_file)

    # Valid single-letter amino acid codes
    valid_codes = set("ACDEFGHIKLMNPQRSTVWYX")

    for chain_id, seq in seqs:
        # All characters should be valid amino acid codes
        assert all(aa in valid_codes for aa in seq), (
            f"Sequence for chain {chain_id} contains invalid amino acid codes"
        )


def test_python_errors_and_conformer_warning(tmp_path):
    """Map failures to Python errors and report automatic conformer choice."""
    import arpeggia

    with pytest.raises(OSError):
        arpeggia.seq(str(tmp_path / "missing.pdb"))
    with pytest.raises(ValueError, match="finite and positive"):
        arpeggia.sasa(
            str(Path(__file__).parents[2] / "test-data/1ubq.pdb"),
            probe_radius=float("nan"),
        )

    alternate = tmp_path / "alternate.pdb"
    alternate.write_text(
        "ATOM      1  CB AALA A   1       0.000   0.000   0.000  0.50 20.00           C  \n"
        "ATOM      2  CB BALA A   1       1.000   0.000   0.000  0.50 20.00           C  \n"
        "END\n"
    )
    with pytest.warns(UserWarning, match="CONFORMER_SELECTED"):
        assert arpeggia.seq(str(alternate)) == [("A", "A")]

    unsupported_sc = tmp_path / "unsupported-sc-radius.pdb"
    unsupported_sc.write_text(
        "ATOM      1  QQ  ALA A   1       0.000   0.000   0.000  1.00 20.00          RN  \n"
        "ATOM      2  CB  ALA B   1       3.000   0.000   0.000  1.00 20.00           C  \n"
        "END\n"
    )
    with pytest.raises(RuntimeError, match="van der Waals radius"):
        arpeggia.sc(str(unsupported_sc), groups="A/B")
    with pytest.raises(RuntimeError, match="van der Waals radius"):
        arpeggia.sasa(str(unsupported_sc))


def test_seqres_is_declared_and_seq_is_observed(tmp_path):
    """Keep coordinate-observed and declared polymer sequences distinct."""
    import arpeggia

    structure = tmp_path / "declared.pdb"
    structure.write_text(
        "SEQRES   1 A    3  ALA MSE GLY\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "END\n"
    )
    assert arpeggia.seq(str(structure)) == [("A", "A")]
    assert arpeggia.seqres(str(structure)) == [("A", "AMG")]


def test_seqres_warns_when_mapping_unknown_monomers(tmp_path):
    """Map unknown declared monomers to X with a visible warning."""
    import arpeggia

    structure = tmp_path / "unknown-seqres.pdb"
    structure.write_text("SEQRES   1 B    1  UNK\nEND\n")
    with pytest.warns(UserWarning, match="UNSUPPORTED_MONOMER"):
        assert arpeggia.seqres(str(structure)) == [("B", "X")]


def test_histidine_modes_and_potential_ionic_category(tmp_path):
    """Distinguish inferred histidine charge from explicit evidence."""
    import arpeggia

    structure = tmp_path / "histidine.pdb"
    structure.write_text(
        "ATOM      1  CG  HIS A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "ATOM      2  OD1 ASP B   1       3.000   0.000   0.000  1.00 20.00           O  \n"
        "END\n"
    )
    with pytest.warns(UserWarning, match="UNRESOLVED_HISTIDINE"):
        compatible = arpeggia.contacts(str(structure), groups="A/B")
    assert "PotentialIonicBond" in compatible["interaction"].to_list()

    with pytest.warns(UserWarning, match="UNRESOLVED_HISTIDINE"):
        explicit = arpeggia.contacts(
            str(structure), groups="A/B", protonation="explicit-only"
        )
    assert "PotentialIonicBond" not in explicit["interaction"].to_list()
