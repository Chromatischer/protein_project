#!/usr/bin/env python3
"""
Test script to verify KEGG integration in main function.
"""

import sys
import os

# Add the src directory to the path so we can import our modules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from models.protein import Protein
from services.protein_service import fetch_kegg_info_batch
from Bio.SeqRecord import SeqRecord


def test_kegg_integration():
    """Test KEGG fetching integration with sample protein data."""
    print("Testing KEGG integration...")
    
    # Create a sample protein object for testing
    # Using a simple SeqRecord for testing
    sample_seq_record = SeqRecord("")
    sample_seq_record.id = "test_protein"
    sample_seq_record.description = "Test protein description"
    
    test_protein = Protein(
        seq_id="test_seq",
        sci_identifier="NP_000001.1",  # Sample protein ID
        desc="Test protein",
        aa_start=1,
        aa_stop=100,
        peptide="TESTPEPTIDE",
        std_dev=0.5,
        info=sample_seq_record,
        xref_id=3569,  # Real NCBI Gene ID for IL6 gene
        kegg_info=None,
        K004A1="1.0",
        K004A2="1.1",
        K225A1="1.2",
        K225A2="1.3",
        K003B1="1.4",
        K003B2="1.5",
        K237A1="1.6",
        K237A2="1.7",
        K022B1="1.8",
        K022B2="1.9",
        K302A1="2.0",
        K302A2="2.1",
    )
    
    protein_list = [test_protein]
    
    try:
        # Test the KEGG batch fetching
        print("Calling fetch_kegg_info_batch...")
        kegg_info_map = fetch_kegg_info_batch(protein_list)
        
        print(f"KEGG info map returned: {type(kegg_info_map)}")
        print(f"Keys in map: {list(kegg_info_map.keys())}")
        
        # Test populating KEGG info
        for protein in protein_list:
            kegg_data = kegg_info_map.get(protein.sci_identifier)
            if kegg_data:
                protein.kegg_info = kegg_data
                print(f"Successfully populated KEGG info for {protein.sci_identifier}")
                print(f"KEGG info type: {type(kegg_data)}")
            else:
                print(f"No KEGG info found for {protein.sci_identifier}")
        
        print("KEGG integration test completed successfully!")
        return True
        
    except Exception as e:
        print(f"Error during KEGG integration test: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_kegg_integration()
    sys.exit(0 if success else 1)

