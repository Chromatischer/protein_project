#!/usr/bin/env python3
"""
Test script for KEGG protein query functions.

This script tests the newly implemented KEGG functionality in protein_service.py
"""

from services.protein_service import (
    fetch_protein_info_kegg,
    search_protein_kegg,
    fetch_protein_info_kegg_batch,
)
import sys
import os

# Add src directory to path to import our modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))


def test_fetch_protein_info_kegg():
    """Test fetching protein information from KEGG for a single protein."""
    print("=== Testing fetch_protein_info_kegg ===")

    # Test with human p53 gene
    protein_id = "hsa:7157"
    print(f"Testing with protein ID: {protein_id}")

    result = fetch_protein_info_kegg(protein_id)

    if "error" in result:
        print(f"Error occurred: {result['error']}")
    else:
        print("Success! Retrieved information:")
        if "entry" in result:
            print(f"Entry data (first 200 chars): {result['entry'][:200]}...")
        if "sequence" in result:
            print(
                f"Sequence data (first 100 chars): {result['sequence'][:100]}...")
        if "pathways" in result:
            print(f"Pathway data: {result['pathways']}")

        # Print any errors for individual components
        for key in result:
            if key.endswith("_error"):
                print(f"Warning - {key}: {result[key]}")

    return result


def test_search_protein_kegg():
    """Test searching for proteins in KEGG database."""
    print("\n=== Testing search_protein_kegg ===")

    # Test searching for p53 in human
    query = "p53"
    organism = "hsa"
    print(f"Searching for '{query}' in organism '{organism}'")

    results = search_protein_kegg(query, organism)

    if results:
        print(f"Found {len(results)} results:")
        for i, result in enumerate(results[:3]):  # Show first 3 results
            print(
                f"  {i + 1}. ID: {result['id']}, Description: {result['description'][:60]}..."
            )
    else:
        print("No results found")

    return results


def test_fetch_protein_info_kegg_batch():
    """Test batch fetching of protein information from KEGG."""
    print("\n=== Testing fetch_protein_info_kegg_batch ===")

    # Test with a small batch of human genes
    protein_ids = ["hsa:7157", "hsa:672"]  # p53 and BRCA1
    print(f"Testing batch fetch with protein IDs: {protein_ids}")

    results = fetch_protein_info_kegg_batch(protein_ids)

    for protein_id, result in results.items():
        print(f"\nResults for {protein_id}:")
        if "error" in result:
            print(f"  Error: {result['error']}")
        else:
            print(f"  Successfully retrieved information")
            if "entry" in result:
                print(
                    f"  Entry data available: {len(result['entry'])} characters")

    return results


def main():
    """Run all KEGG tests."""
    print("Starting KEGG protein service tests...\n")

    try:
        # Test individual function
        test_fetch_protein_info_kegg()

        # Test search function
        test_search_protein_kegg()

        # Test batch function
        test_fetch_protein_info_kegg_batch()

        print("\n=== All tests completed ===")

    except Exception as e:
        print(f"Test failed with error: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
