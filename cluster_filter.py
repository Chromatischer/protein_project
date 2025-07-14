#!/usr/bin/env python3
"""
Script to extract clusters from cluster_export.json while filtering out long genome sequences.

This script removes or truncates very long sequence fields like 'aaseq' and 'ntseq'
to make the data more manageable while preserving all other important metadata.
"""

import json
import argparse
import sys
from typing import Dict, Any, List


def truncate_long_sequences(
    data: Dict[str, Any], max_length: int = 100
) -> Dict[str, Any]:
    """
    Recursively traverse the data structure and truncate long sequence fields.

    Args:
        data: The data structure to process
        max_length: Maximum length for sequence fields before truncation

    Returns:
        Processed data with truncated sequences
    """
    if isinstance(data, dict):
        result = {}
        for key, value in data.items():
            if (
                key in ["aaseq", "ntseq"]
                and isinstance(value, str)
                and len(value) > max_length
            ):
                # Truncate long sequences and add indication of truncation
                result[key] = (
                    value[:max_length]
                    + f"... [TRUNCATED - original length: {len(value)}]"
                )
            elif isinstance(value, (dict, list)):
                result[key] = truncate_long_sequences(value, max_length)
            else:
                result[key] = value
        return result
    elif isinstance(data, list):
        return [truncate_long_sequences(item, max_length) for item in data]
    else:
        return data


def remove_sequence_fields(
    data: Dict[str, Any], fields_to_remove: List[str] = None
) -> Dict[str, Any]:
    """
    Recursively traverse the data structure and remove specified sequence fields.

    Args:
        data: The data structure to process
        fields_to_remove: List of field names to remove completely

    Returns:
        Processed data with specified fields removed
    """
    if fields_to_remove is None:
        fields_to_remove = ["aaseq", "ntseq", "references"]

    if isinstance(data, dict):
        result = {}
        for key, value in data.items():
            if key not in fields_to_remove:
                if isinstance(value, (dict, list)):
                    result[key] = remove_sequence_fields(value, fields_to_remove)
                else:
                    result[key] = value
        return result
    elif isinstance(data, list):
        return [remove_sequence_fields(item, fields_to_remove) for item in data]
    else:
        return data


def get_file_size_mb(filepath: str) -> float:
    """Get file size in megabytes."""
    try:
        import os

        return os.path.getsize(filepath) / (1024 * 1024)
    except Exception as e:
        raise e
        return 0


def main():
    parser = argparse.ArgumentParser(
        description="Filter cluster export data by removing or truncating long sequences"
    )
    parser.add_argument("input_file", help="Input JSON file (cluster_export.json)")
    parser.add_argument("output_file", help="Output JSON file for filtered data")
    parser.add_argument(
        "--mode",
        choices=["truncate", "remove"],
        default="truncate",
        help="Mode: truncate sequences or remove them completely (default: truncate)",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=100,
        help="Maximum length for sequences when truncating (default: 100)",
    )
    parser.add_argument(
        "--fields",
        nargs="+",
        default=["aaseq", "ntseq", "references"],
        help="Fields to process (default: aaseq ntseq)",
    )
    parser.add_argument(
        "--pretty", action="store_true", help="Pretty print the output JSON"
    )

    args = parser.parse_args()

    try:
        # Read input file
        print(f"Reading {args.input_file}...")
        original_size = get_file_size_mb(args.input_file)
        print(f"Original file size: {original_size:.2f} MB")

        with open(args.input_file, "r", encoding="utf-8") as f:
            data = json.load(f)

        print(f"Loaded data with {len(data)} clusters")

        # Process data based on mode
        if args.mode == "truncate":
            print(f"Truncating sequence fields to max {args.max_length} characters...")
            filtered_data = truncate_long_sequences(data, args.max_length)
        else:  # remove
            print(f"Removing sequence fields: {args.fields}")
            filtered_data = remove_sequence_fields(data, args.fields)

        # Write output file
        print(f"Writing filtered data to {args.output_file}...")
        with open(args.output_file, "w", encoding="utf-8") as f:
            if args.pretty:
                json.dump(filtered_data, f, indent=2, ensure_ascii=False)
            else:
                json.dump(filtered_data, f, ensure_ascii=False)

        # Show size comparison
        new_size = get_file_size_mb(args.output_file)
        print(f"Filtered file size: {new_size:.2f} MB")
        if original_size > 0:
            reduction = ((original_size - new_size) / original_size) * 100
            print(f"Size reduction: {reduction:.1f}%")

        print("✅ Filtering completed successfully!")

    except FileNotFoundError:
        print(f"❌ Error: Input file '{args.input_file}' not found.")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"❌ Error: Invalid JSON in input file - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
