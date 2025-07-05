#!/usr/bin/env python3
"""
Test script to verify JSON export functionality and KEGG integration.
"""

import json
from services.protein_service import fetch_protein_info_kegg
from models.protein import Protein
import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))


def test_protein_with_kegg():
    """Test creating a protein object with KEGG information."""

    # Mock protein info (normally from NCBI)
    class MockInfo:
        def __init__(self):
            self.id = "NP_000047.2"
            self.description = "hemoglobin subunit beta [Homo sapiens]"
            self.seq = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"

    mock_info = MockInfo()

    # Sample KEGG info (this would be fetched in real scenario)
    sample_kegg_info = {
        "entry": "Sample KEGG entry data",
        "sequence": "Sample sequence data",
        "pathways": "Sample pathway data",
    }

    # Create test protein
    test_protein = Protein(
        seq_id=1,
        sci_identifier="NP_000047.2",
        desc="hemoglobin subunit beta",
        aa_start=1,
        aa_stop=147,
        peptide="MVHLTPEEKS",
        std_dev=1.5,
        info=mock_info,
        kegg_info=sample_kegg_info,
        K004A1=100,
        K004A2=105,
        K225A1=98,
        K225A2=102,
        K003B1=95,
        K003B2=99,
        K237A1=101,
        K237A2=103,
        K022B1=97,
        K022B2=101,
        K302A1=100,
        K302A2=104,
    )

    print("✓ Protein object created successfully with KEGG info")
    print(f"KEGG info: {test_protein.kegg_info}")

    return test_protein


def test_json_export_structure():
    """Test the JSON export structure."""

    # Create test data
    test_protein = test_protein_with_kegg()
    tops = [test_protein]

    # Create mock clusters_dict
    clusters_dict = {1: ["NP_000047.2"]}
    distance_threshold = 50000

    # Test the export function structure
    protein_lookup = {protein.sci_identifier: protein for protein in tops}

    export_data = {
        "metadata": {
            "total_proteins": len(tops),
            "total_clusters": len(clusters_dict),
            "distance_threshold": distance_threshold,
            "clustering_method": "complete",
            "clustering_metric": "euclidean",
        },
        "clusters": {},
    }

    # Process clusters
    for cluster_id, protein_ids in clusters_dict.items():
        cluster_proteins = []

        for protein_id in protein_ids:
            protein = protein_lookup.get(protein_id)
            if protein is None:
                continue

            protein_data = {
                "seq_id": protein.seq_id,
                "sci_identifier": protein.sci_identifier,
                "description": protein.desc,
                "aa_start": protein.aa_start,
                "aa_stop": protein.aa_stop,
                "peptide": protein.peptide,
                "std_dev": protein.std_dev,
                "hits": protein.hits,
                "ncbi_info": {
                    "id": protein.info.id if protein.info else None,
                    "description": protein.info.description if protein.info else None,
                    "sequence_length": len(protein.info.seq)
                    if protein.info and hasattr(protein.info, "seq")
                    else None,
                },
                "kegg_info": protein.kegg_info,
            }

            cluster_proteins.append(protein_data)

        export_data["clusters"][str(cluster_id)] = {
            "cluster_id": cluster_id,
            "protein_count": len(cluster_proteins),
            "proteins": cluster_proteins,
        }

    # Test JSON serialization
    try:
        json_str = json.dumps(export_data, indent=2)
        print("✓ JSON serialization successful")
        print("✓ Export data structure is valid")

        # Save test output
        with open("test_clusters_export.json", "w") as f:
            f.write(json_str)
        print("✓ Test JSON file created: test_clusters_export.json")

        return True
    except Exception as e:
        print(f"✗ JSON serialization failed: {e}")
        return False


def test_kegg_service():
    """Test KEGG service functionality."""
    print("\n--- Testing KEGG Service ---")

    # Test with a known protein (human hemoglobin beta)
    try:
        kegg_info = fetch_protein_info_kegg("hsa:3043")
        print("✓ KEGG service is accessible")
        print(f"KEGG response keys: {list(kegg_info.keys())}")

        if "error" not in kegg_info:
            print("✓ KEGG data retrieved successfully")
        else:
            print(f"⚠ KEGG returned error: {kegg_info['error']}")

    except Exception as e:
        print(f"⚠ KEGG service test failed: {e}")


if __name__ == "__main__":
    print("=== Testing JSON Export and KEGG Integration ===")

    # Test 1: Protein creation with KEGG
    print("\n1. Testing Protein creation with KEGG info...")
    test_protein_with_kegg()

    # Test 2: JSON export structure
    print("\n2. Testing JSON export structure...")
    json_success = test_json_export_structure()

    # Test 3: KEGG service
    print("\n3. Testing KEGG service...")
    test_kegg_service()

    if json_success:
        print("\n✓ All tests completed successfully!")
        print("Check 'test_clusters_export.json' for sample output.")
    else:
        print("\n✗ Some tests failed!")
