#!/usr/bin/env python3
"""
Script to filter proteins from clusters_export.json that have XP identifiers (but not NP)
and export them to a new JSON file with the same structure.
"""

import json
import sys
from datetime import datetime


def filter_xp_proteins(
    input_file="clusters_export.json", output_file="xp_proteins_export.json"
):
    """
    Filter proteins with XP identifiers (but not NP) from the clusters export JSON.

    Args:
        input_file (str): Path to the input JSON file
        output_file (str): Path to the output JSON file
    """

    try:
        # Load the original data
        with open(input_file, "r") as f:
            data = json.load(f)

        print(f"Loaded data from {input_file}")
        print(
            f"Original data contains {data['metadata']['total_proteins']} proteins in {data['metadata']['total_clusters']} clusters"
        )

        # Create new data structure with filtered proteins
        filtered_data = {"metadata": data["metadata"].copy(), "clusters": {}}

        total_filtered_proteins = 0
        filtered_clusters = 0

        # Process each cluster
        for cluster_id, cluster_data in data["clusters"].items():
            filtered_proteins = []

            # Filter proteins in this cluster
            for protein in cluster_data["proteins"]:
                sci_id = protein["sci_identifier"]

                # Check if protein has XP identifier but not NP
                if "XP_" in sci_id and "NP_" not in sci_id:
                    filtered_proteins.append(protein)
                    total_filtered_proteins += 1

            # Only include cluster if it has filtered proteins
            if filtered_proteins:
                filtered_clusters += 1
                filtered_data["clusters"][cluster_id] = {
                    "cluster_id": cluster_data["cluster_id"],
                    "protein_count": len(filtered_proteins),
                    "proteins": filtered_proteins,
                }

        # Update metadata with filtered counts
        filtered_data["metadata"]["total_proteins"] = total_filtered_proteins
        filtered_data["metadata"]["total_clusters"] = filtered_clusters
        filtered_data["metadata"]["filter_applied"] = (
            "XP identifiers only (excluding NP)"
        )
        filtered_data["metadata"]["filter_timestamp"] = datetime.now().isoformat()
        filtered_data["metadata"]["original_file"] = input_file

        # Save filtered data
        with open(output_file, "w") as f:
            json.dump(filtered_data, f, indent=2)

        print(f"\nFiltering complete!")
        print(
            f"Filtered {total_filtered_proteins} proteins with XP identifiers (no NP)"
        )
        print(f"Retained {filtered_clusters} clusters containing XP proteins")
        print(f"Results saved to {output_file}")

        return filtered_data

    except FileNotFoundError:
        print(f"Error: Could not find input file '{input_file}'")
        return None
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON format in '{input_file}'")
        return None
    except Exception as e:
        print(f"Error: {str(e)}")
        return None


def main():
    """Main function to run the filtering process."""

    # Check command line arguments
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = "clusters_export.json"

    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    else:
        output_file = "xp_proteins_export.json"

    print("XP Protein Filter")
    print("=" * 50)
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    print()

    # Run the filtering
    result = filter_xp_proteins(input_file, output_file)

    if result is None:
        sys.exit(1)

    print("\nDone!")


if __name__ == "__main__":
    main()
