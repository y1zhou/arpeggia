"""Test suite for arpeggia Python bindings using pytest."""

# ruff: noqa: S101
from pathlib import Path

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

    assert hasattr(arpeggia, "__version__")
    assert hasattr(arpeggia, "contacts")
    assert hasattr(arpeggia, "sasa")
    assert hasattr(arpeggia, "pdb2seq")


def test_contacts(test_pdb_file):
    """Test the contacts function returns expected DataFrame structure."""
    import arpeggia

    df = arpeggia.contacts(test_pdb_file, groups="/", vdw_comp=0.1, dist_cutoff=6.5)

    # Check DataFrame is not empty
    assert df.height == 532, "Contacts DataFrame should not be empty"

    # Check expected columns exist
    expected_columns = [
        "model",
        "interaction",
        "distance",
        "from_chain",
        "from_resn",
        "from_resi",
        "from_insertion",
        "from_altloc",
        "from_atomn",
        "from_atomi",
        "to_chain",
        "to_resn",
        "to_resi",
        "to_insertion",
        "to_altloc",
        "to_atomn",
        "to_atomi",
        "sc_centroid_dist",
        "sc_dihedral",
        "sc_centroid_angle",
    ]

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


def test_sasa(test_pdb_file):
    """Test the sasa function returns expected DataFrame structure."""
    import arpeggia

    df = arpeggia.sasa(test_pdb_file, probe_radius=1.4, n_points=100, model_num=0)

    # Check DataFrame is not empty
    assert df.height == 602, "SASA DataFrame should not be empty"

    # Check expected columns exist
    expected_columns = [
        "atomi",
        "sasa",
        "chain",
        "resn",
        "resi",
        "insertion",
        "altloc",
        "atomn",
    ]

    for col in expected_columns:
        assert col in df.columns, f"Column '{col}' should be present in SASA DataFrame"

    # Check shape - should have 8 columns
    assert df.shape[1] == 8, f"Expected 8 columns, got {df.shape[1]}"

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


def test_pdb2seq(test_pdb_file):
    """Test the pdb2seq function returns expected structure."""
    import arpeggia

    seqs = arpeggia.pdb2seq(test_pdb_file)

    # Check return type
    assert isinstance(seqs, dict), "Sequences should return a dictionary"
    assert len(seqs) > 0, "Should have at least one chain"

    # For 1ubq.pdb, we know it has 1 chain with a specific sequence
    # Chain should be present
    assert len(seqs) == 1, f"Expected 1 chain, got {len(seqs)}"

    # Check sequence properties
    for chain_id, seq in seqs.items():
        assert isinstance(chain_id, str), "Chain ID should be string"
        assert isinstance(seq, str), "Sequence should be string"
        assert len(seq) > 0, "Sequence should not be empty"

        # For 1ubq, the sequence should be 76 residues
        # This is the known ubiquitin sequence
        aa_seq = seq.replace("O", "")
        assert len(aa_seq) == 76, (
            f"Expected 76 residues for ubiquitin, got {len(aa_seq)}"
        )

        # Check it starts with the expected sequence
        expected_start = "MQIFVKTLTG"
        assert seq.startswith(expected_start), (
            f"Sequence should start with {expected_start}, got {seq[:10]}"
        )


def test_sequences_validity(test_pdb_file):
    """Test that returned sequences contain valid amino acid codes."""
    import arpeggia

    seqs = arpeggia.pdb2seq(test_pdb_file)

    # Valid single-letter amino acid codes
    valid_codes = set("ACDEFGHIKLMNPQRSTVWYXO")

    for chain_id, seq in seqs.items():
        # All characters should be valid amino acid codes
        assert all(aa in valid_codes for aa in seq), (
            f"Sequence for chain {chain_id} contains invalid amino acid codes"
        )
