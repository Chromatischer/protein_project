#!/usr/bin/env python3
"""
CSV Data Grouping Script

This script separates CSV raw data into groups of 4 columns and groups proteins by:
- Proteins in Group 1&2 but not 3&4
- Proteins in Group 3&4 but not 1&2

Groups:
- "Group_1": ["K004A1", "K004A2"],
- "Group_2": ["K225A1", "K225A2"],
- "Group_3": ["K003B1", "K003B2"],
- "Group_4": ["K237A1", "K237A2"],
- "Group_5": ["K022B1", "K022B2"],
- "Group_6": ["K302A1", "K302A2"],

"""

import pandas as pd
import numpy as np
import json
from pathlib import Path


def load_csv_data(csv_path):
    """Load the CSV data and return DataFrame"""
    return pd.read_csv(csv_path)


def define_column_groups():
    """Define the column groups as specified in the task"""
    groups = {
        "Group_1": ["K004A1", "K004A2"],
        "Group_2": ["K225A1", "K225A2"],
        "Group_3": ["K003B1", "K003B2"],
        "Group_4": ["K237A1", "K237A2"],
        "Group_5": ["K022B1", "K022B2"],
        "Group_6": ["K302A1", "K302A2"],
    }
    return groups


def has_hits_in_group(row, group_columns, threshold=0):
    """Check if a protein has hits (non-zero values) in any column of a group"""
    group_values = row[group_columns]
    return any(val > threshold for val in group_values if pd.notna(val))


def separate_proteins_by_groups(df, groups):
    """
    Separate proteins into categories:
    1. Proteins in Group 1&2 but not 3&4
    2. Proteins in Group 3&4 but not 1&2
    3. Proteins in both Group 1&2 and Group 3&4
    4. Proteins with no hits in any group
    """

    results = {
        "group_1_2_only": [],
        "group_3_4_only": [],
        "group_5_6_only": [],
        "multi_1_2_3_4": [],
        "multi_1_2_5_6": [],
        "multi_3_4_5_6": [],
        "multi_1_2_3_4_5_6": [],
        "no_hits": [],
    }

    for idx, row in df.iterrows():
        # Check hits in each group
        has_group_1 = has_hits_in_group(row, groups["Group_1"])
        has_group_2 = has_hits_in_group(row, groups["Group_2"])
        has_group_3 = has_hits_in_group(row, groups["Group_3"])
        has_group_4 = has_hits_in_group(row, groups["Group_4"])
        has_group_5 = has_hits_in_group(row, groups["Group_5"])
        has_group_6 = has_hits_in_group(row, groups["Group_6"])

        # Create protein info
        protein_info = {
            "SeqID": row["SeqID"],
            "Protein": row["Protein"],
            "Description": row["Description"],
            "Peptide": row["Peptide"],
            "Group_1_hits": has_group_1,
            "Group_2_hits": has_group_2,
            "Group_3_hits": has_group_3,
            "Group_4_hits": has_group_4,
            "Group_5_hits": has_group_5,
            "Group_6_hits": has_group_6,
            "Group_1_values": {col: row[col] for col in groups["Group_1"]},
            "Group_2_values": {col: row[col] for col in groups["Group_2"]},
            "Group_3_values": {col: row[col] for col in groups["Group_3"]},
            "Group_4_values": {col: row[col] for col in groups["Group_4"]},
            "Group_5_values": {col: row[col] for col in groups["Group_5"]},
            "Group_6_values": {col: row[col] for col in groups["Group_6"]},
        }

        # Categorize proteins
        has_group_1_or_2 = has_group_1 or has_group_2
        has_group_3_or_4 = has_group_3 or has_group_4
        has_group_5_or_6 = has_group_5 or has_group_6

        if has_group_1_or_2 and not has_group_3_or_4 and not has_group_5_or_6:
            results["group_1_2_only"].append(protein_info)
        elif has_group_3_or_4 and not has_group_1_or_2 and not has_group_5_or_6:
            results["group_3_4_only"].append(protein_info)
        elif has_group_5_or_6 and not has_group_3_or_4 and not has_group_1_or_2:
            results["group_5_6_only"].append(protein_info)
        elif has_group_1_or_2 and has_group_3_or_4:
            results["multi_1_2_3_4"].append(protein_info)
        elif has_group_1_or_2 and has_group_5_or_6:
            results["multi_1_2_5_6"].append(protein_info)
        elif has_group_3_or_4 and has_group_5_or_6:
            results["multi_3_4_5_6"].append(protein_info)
        else:
            results["no_hits"].append(protein_info)

    return results


def create_separated_csv_files(df, groups, output_dir="output"):
    """Create separate CSV files for each group"""
    Path(output_dir).mkdir(exist_ok=True)

    metadata_cols = [
        "SeqID",
        "Protein",
        "Description",
        "AA_Start",
        "AA_Stop",
        "Peptide",
        "std_dev",
    ]

    # Group 1 CSV
    group_1_cols = metadata_cols + groups["Group_1"]
    df_group_1 = df[group_1_cols].copy()
    df_group_1.to_csv(f"{output_dir}/group_1_K004.csv", index=False)

    # Group 2 CSV
    group_2_cols = metadata_cols + groups["Group_2"]
    df_group_2 = df[group_2_cols].copy()
    df_group_2.to_csv(f"{output_dir}/group_2_K225.csv", index=False)

    # Group 3 CSV
    group_3_cols = metadata_cols + groups["Group_3"]
    df_group_3 = df[group_3_cols].copy()
    df_group_3.to_csv(f"{output_dir}/group_3_K003.csv", index=False)

    # Group 4 CSV
    group_4_cols = metadata_cols + groups["Group_4"]
    df_group_4 = df[group_4_cols].copy()
    df_group_4.to_csv(f"{output_dir}/group_4_K237.csv", index=False)

    # Group 5 CSV
    group_5_cols = metadata_cols + groups["Group_5"]
    df_group_5 = df[group_5_cols].copy()
    df_group_5.to_csv(f"{output_dir}/group_5_K022.csv", index=False)

    # Group 6 CSV
    group_6_cols = metadata_cols + groups["Group_6"]
    df_group_6 = df[group_6_cols].copy()
    df_group_6.to_csv(f"{output_dir}/group_6_K302.csv", index=False)

    return {
        "group_1": f"{output_dir}/group_1_K004.csv",
        "group_2": f"{output_dir}/group_2_K225.csv",
        "group_3": f"{output_dir}/group_3_K003.csv",
        "group_4": f"{output_dir}/group_4_K237.csv",
        "group_5": f"{output_dir}/group_5_K022.csv",
        "group_6": f"{output_dir}/group_6_K302.csv",
    }


def print_summary(results):
    """Print a summary of the protein grouping results"""
    print("=== PROTEIN GROUPING SUMMARY ===")
    print(
        f"Proteins in Group 1&2 but NOT in Group 3&4: {len(results['group_1_2_only'])}"
    )
    print(
        f"Proteins in Group 3&4 but NOT in Group 1&2: {len(results['group_3_4_only'])}"
    )
    print(f"Proteins in BOTH Group 1&2 AND Group 3&4: {len(results['multi_1_2_3_4'])}")
    print(f"Proteins in BOTH Group 3&4 AND Group 5&6: {len(results['multi_3_4_5_6'])}")
    print(f"Proteins in BOTH Group 1&2 AND Group 5&6: {len(results['multi_1_2_5_6'])}")
    print(f"Proteins with NO hits in any group: {len(results['no_hits'])}")
    print()

    # Show some examples
    if results["group_1_2_only"]:
        print("Examples of proteins in Group 1&2 only:")
        for i, protein in enumerate(results["group_1_2_only"][:3]):
            print(f"  {i + 1}. {protein['Protein']} - {protein['Description'][:50]}...")
        print()

    if results["group_3_4_only"]:
        print("Examples of proteins in Group 3&4 only:")
        for i, protein in enumerate(results["group_3_4_only"][:3]):
            print(f"  {i + 1}. {protein['Protein']} - {protein['Description'][:50]}...")
        print()


def main():
    """Main function to execute the CSV grouping task"""
    # Load data
    csv_path = "resources/WWOP230228_PHSN_C12-2_Top100_Variant_Peptides_blinded.csv"
    df = load_csv_data(csv_path)

    print(f"Loaded CSV with {len(df)} proteins")
    print("Columns:", list(df.columns))
    print()

    # Define groups
    groups = define_column_groups()
    print("Column Groups:")
    for group_name, cols in groups.items():
        print(f"  {group_name}: {cols}")
    print()

    # Separate proteins by groups
    results = separate_proteins_by_groups(df, groups)

    # Print summary
    print_summary(results)

    # Create separate CSV files
    output_files = create_separated_csv_files(df, groups)
    print("Created separate CSV files:")
    for group_name, file_path in output_files.items():
        print(f"  {group_name}: {file_path}")
    print()

    # Save detailed results to JSON
    output_dir = Path("output")
    output_dir.mkdir(exist_ok=True)

    with open("output/protein_grouping_results.json", "w") as f:
        json.dump(results, f, indent=2)

    print("Detailed results saved to: output/protein_grouping_results.json")

    return results


if __name__ == "__main__":
    results = main()
