"""Simple test script to verify arpeggia Python bindings work correctly."""

import sys
from pathlib import Path

def test_import():
    """Test that the module can be imported."""
    print("Testing import...")
    try:
        import arpeggia
        print(f"✓ Successfully imported arpeggia version {arpeggia.__version__}")
        return True
    except ImportError as e:
        print(f"✗ Failed to import arpeggia: {e}")
        return False

def test_contacts():
    """Test the contacts function."""
    print("\nTesting contacts function...")
    try:
        import arpeggia
        
        # Path to test data
        test_file = Path(__file__).parent.parent / "test-data" / "1ubq.pdb"
        
        if not test_file.exists():
            print(f"✗ Test file not found: {test_file}")
            return False
        
        # Call the contacts function
        df = arpeggia.contacts(str(test_file), groups="/", vdw_comp=0.1, dist_cutoff=6.5)
        
        print(f"✓ Found {len(df)} contacts")
        print(f"  DataFrame shape: {df.shape}")
        print(f"  Columns: {df.columns[:5]}...")  # Show first 5 columns
        
        return True
    except Exception as e:
        print(f"✗ Error testing contacts: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_sasa():
    """Test the sasa function."""
    print("\nTesting sasa function...")
    try:
        import arpeggia
        
        # Path to test data
        test_file = Path(__file__).parent.parent / "test-data" / "1ubq.pdb"
        
        if not test_file.exists():
            print(f"✗ Test file not found: {test_file}")
            return False
        
        # Call the sasa function
        df = arpeggia.sasa(str(test_file), probe_radius=1.4, n_points=100)
        
        print(f"✓ Calculated SASA for {len(df)} atoms")
        print(f"  DataFrame shape: {df.shape}")
        print(f"  Columns: {df.columns}")
        
        return True
    except Exception as e:
        print(f"✗ Error testing sasa: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_sequences():
    """Test the sequences function."""
    print("\nTesting sequences function...")
    try:
        import arpeggia
        
        # Path to test data
        test_file = Path(__file__).parent.parent / "test-data" / "1ubq.pdb"
        
        if not test_file.exists():
            print(f"✗ Test file not found: {test_file}")
            return False
        
        # Call the sequences function
        seqs = arpeggia.sequences(str(test_file))
        
        print(f"✓ Extracted sequences for {len(seqs)} chains")
        for chain_id, seq in seqs.items():
            print(f"  Chain {chain_id}: {len(seq)} residues")
        
        return True
    except Exception as e:
        print(f"✗ Error testing sequences: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests."""
    print("=" * 60)
    print("Arpeggia Python Bindings Test Suite")
    print("=" * 60)
    
    results = {
        "import": test_import(),
        "contacts": test_contacts(),
        "sasa": test_sasa(),
        "sequences": test_sequences(),
    }
    
    print("\n" + "=" * 60)
    print("Test Results:")
    print("=" * 60)
    
    for test_name, passed in results.items():
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"{test_name:20s}: {status}")
    
    all_passed = all(results.values())
    print("=" * 60)
    
    if all_passed:
        print("All tests passed!")
        return 0
    else:
        print("Some tests failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())
