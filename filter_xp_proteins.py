#!/usr/bin/env python3
"""
Enhanced script to filter proteins that have XP identifiers (but not NP) using Protein struct
and reading data from the cache. This version includes full NCBI and KEGG information.
"""

import json
import sys
from datetime import datetime
from pathlib import Path

# Add src to path to import modules
sys.path.insert(0, "src")

try:
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from models.protein import Protein
    from services.protein_service import fetch_protein_info_batch, fetch_kegg_info_batch
    from utils.cache import get_cache
except ImportError as e:
    print(f"Import error: {e}")
    print("Make sure you're running this from the project root directory")
    sys.exit(1)


def load_proteins_from_cache() -> list:
    """
    Load all proteins from cache and reconstruct Protein objects.

    Returns:
        list: List of Protein objects loaded from cache
    """
    print("Loading proteins from cache...")

    # Look for processed protein data
    protein_files = [
        "proteins_export.json",
        "filtered_II.json",
        "clusters_filtered.json",
        "clusters_export.json",
    ]

    for filename in protein_files:
        file_path = Path(filename)
        if file_path.exists():
            print(f"Loading from {filename}...")
            try:
                with open(file_path, "r") as f:
                    data = json.load(f)

                proteins = []

                # Handle different file formats
                if "proteins" in data:
                    # Direct protein list format
                    for protein_data in data["proteins"]:
                        protein = create_protein_from_dict(protein_data)
                        if protein:
                            proteins.append(protein)
                elif "clusters" in data:
                    # Cluster format - extract proteins from clusters
                    for cluster_id, cluster_data in data["clusters"].items():
                        for protein_data in cluster_data.get("proteins", []):
                            protein = create_protein_from_dict(protein_data)
                            if protein:
                                proteins.append(protein)

                return proteins

            except Exception as e:
                print(f"Error loading {filename}: {e}")
                continue

    return []


def create_protein_from_dict(data: dict) -> Protein:
    """
    Create a Protein object from dictionary data.

    Args:
        data (dict): Protein data dictionary

    Returns:
        Protein: Protein object or None if creation failed
    """
    try:
        # Handle both flat and nested structures
        if "protein" in data:
            # Nested structure (from clusters)
            protein_data = data["protein"]
        else:
            # Flat structure
            protein_data = data

        # Extract required fields
        seq_id = protein_data.get("seq_id", 0)
        sci_identifier = protein_data.get("sci_identifier", "")
        desc = protein_data.get("desc", "")
        aa_start = protein_data.get("aa_start", 0)
        aa_stop = protein_data.get("aa_stop", 0)
        peptide = protein_data.get("peptide", "")
        std_dev = float(protein_data.get("std_dev", 0.0))
        xref_id = int(protein_data.get("xref_id", 0))

        # Extract hit values
        hits = protein_data.get("hits", {})
        if not hits:
            # Try to extract from flat structure
            hit_fields = [
                "K004A1",
                "K004A2",
                "K225A1",
                "K225A2",
                "K003B1",
                "K003B2",
                "K237A1",
                "K237A2",
                "K022B1",
                "K022B2",
                "K302A1",
                "K302A2",
            ]
            for field in hit_fields:
                if field in protein_data:
                    hits[field] = int(str(protein_data[field]).replace(".", ""))

        # Extract additional info
        info = protein_data.get("info", {})
        ncbi_info = protein_data.get("ncbi_info", None)
        kegg_info = protein_data.get("kegg_info", None)

        # Create protein with all necessary data
        protein_kwargs = {}
        for key, value in hits.items():
            if value is not None:
                protein_kwargs[key] = value

        protein = Protein(
            seq_id=seq_id,
            sci_identifier=sci_identifier,
            desc=desc,
            aa_start=aa_start,
            aa_stop=aa_stop,
            peptide=peptide,
            std_dev=std_dev,
            info=info,
            xref_id=xref_id,
            **protein_kwargs,
        )

        # Add NCBI and KEGG info if available
        if ncbi_info:
            protein.ncbi_info = ncbi_info
        if kegg_info:
            protein.kegg_info = kegg_info

        return protein

    except Exception as e:
        print(f"Error creating protein: {e}")
        return None


def load_proteins_from_csv() -> list:
    """
    Load proteins from CSV file using cache.

    Returns:
        list: List of Protein objects
    """
    print("Loading proteins from CSV file...")

    # Find CSV file
    csv_files = list(Path(".").glob("*.csv"))
    if not csv_files:
        print("No CSV files found")
        return []

    csv_file = csv_files[0]
    print(f"Using {csv_file}")

    try:
        df = pd.read_csv(csv_file)

        # Check if required columns exist
        required_cols = [
            "Protein",
            "Description",
            "AA_Start",
            "AA_Stop",
            "Peptide",
            "std_dev",
        ]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Missing columns: {missing_cols}")
            return []

        # Get unique protein identifiers
        protein_ids = df["Protein"].unique().tolist()
        print(f"Found {len(protein_ids)} unique proteins")

        # Use cache/service to get protein info
        proteins = []

        for _, row in df.iterrows():
            try:
                # Create basic protein structure
                protein_kwargs = {}
                hit_fields = [
                    "K004A1",
                    "K004A2",
                    "K225A1",
                    "K225A2",
                    "K003B1",
                    "K003B2",
                    "K237A1",
                    "K237A2",
                    "K022B1",
                    "K022B2",
                    "K302A1",
                    "K302A2",
                ]

                for field in hit_fields:
                    if field in row and pd.notna(row[field]):
                        try:
                            protein_kwargs[field] = int(
                                str(row[field]).replace(".", "")
                            )
                        except:
                            protein_kwargs[field] = 0

                protein = Protein(
                    seq_id=int(row.get("SeqID", 0)),
                    sci_identifier=str(row["Protein"]),
                    desc=str(row["Description"]),
                    aa_start=int(row["AA_Start"]),
                    aa_stop=int(row["AA_Stop"]),
                    peptide=str(row["Peptide"]),
                    std_dev=float(row["std_dev"]),
                    info={},
                    xref_id=0,  # Will be populated from cache if available
                    **protein_kwargs,
                )

                proteins.append(protein)

            except Exception as e:
                print(f"Error processing row: {e}")
                continue

        return proteins

    except Exception as e:
        print(f"Error loading CSV: {e}")
        return []


def filter_xp_proteins(proteins: list) -> list:
    """
    Filter proteins with XP identifiers but not NP identifiers.

    Args:
        proteins (list): List of Protein objects to filter

    Returns:
        list: Filtered list of proteins with XP identifiers
    """
    filtered = []

    for protein in proteins:
        sci_id = str(protein.sci_identifier).upper()

        # Check for XP identifier without NP
        if "XP_" in sci_id and "NP_" not in sci_id:
            filtered.append(protein)

    return filtered


def fetch_additional_info(proteins: list) -> list:
    """
    Fetch NCBI and KEGG information for proteins missing this data.

    Args:
        proteins (list): List of Protein objects

    Returns:
        list: List of proteins with updated information
    """
    print("Fetching additional NCBI and KEGG information...")

    # Get protein identifiers that need info
    identifiers = [
        p.sci_identifier
        for p in proteins
        if not hasattr(p, "ncbi_info") or not hasattr(p, "kegg_info")
    ]

    if not identifiers:
        print("All proteins already have NCBI/KEGG information")
        return proteins

    try:
        from services.protein_service import (
            fetch_protein_info_batch,
            fetch_kegg_info_batch,
        )

        # Fetch NCBI info
        print(f"Fetching NCBI info for {len(identifiers)} proteins...")
        ncbi_results = fetch_protein_info_batch(identifiers)

        # Fetch KEGG info
        print(f"Fetching KEGG info for {len(identifiers)} proteins...")
        kegg_results = fetch_kegg_info_batch(identifiers)

        # Update proteins with fetched data
        for protein in proteins:
            if protein.sci_identifier in ncbi_results:
                protein.ncbi_info = ncbi_results[protein.sci_identifier]
            if protein.sci_identifier in kegg_results:
                protein.kegg_info = kegg_results[protein.sci_identifier]

    except ImportError as e:
        print(f"Could not fetch additional info: {e}")
        print("Proceeding with available information...")

    return proteins


def save_filtered_proteins(proteins: list, output_file: str = "xp_proteins_cache.json"):
    """
    Save filtered proteins to JSON file with full NCBI and KEGG information.

    Args:
        proteins (list): List of filtered Protein objects
        output_file (str): Output filename
    """
    protein_list = []

    for protein in proteins:
        # Create detailed protein entry with full database information
        protein_entry = {
            "seq_id": protein.seq_id,
            "sci_identifier": protein.sci_identifier,
            "desc": protein.desc,
            "aa_start": protein.aa_start,
            "aa_stop": protein.aa_stop,
            "peptide": protein.peptide,
            "std_dev": protein.std_dev,
            "xref_id": protein.xref_id,
            "hits": protein.hits,
            "basic_info": protein.info,
            "full_ncbi_info": None,
            "full_kegg_info": None,
            "database_cross_references": {},
        }

        # Add full NCBI information if available
        if hasattr(protein, "ncbi_info") and protein.ncbi_info:
            if protein.ncbi_info is not None:
                if hasattr(protein.ncbi_info, "__dict__"):
                    protein_entry["full_ncbi_info"] = protein.ncbi_info.__dict__
                else:
                    protein_entry["full_ncbi_info"] = protein.ncbi_info

        # Add full KEGG information if available
        if hasattr(protein, "kegg_info") and protein.kegg_info:
            if protein.kegg_info is not None:
                if hasattr(protein.kegg_info, "__dict__"):
                    protein_entry["full_kegg_info"] = protein.kegg_info.__dict__
                else:
                    protein_entry["full_kegg_info"] = protein.kegg_info

        # Add cross-references from info
        if protein.info:
            # Extract common database cross-references
            db_refs = {}
            for key, value in protein.info.items():
                if any(
                    db_key in str(key).lower()
                    for db_key in ["geneid", "hgnc", "mim", "taxon", "pdb", "go"]
                ):
                    db_refs[key] = value
                elif (
                    "cross_reference" in str(key).lower()
                    or "db_xref" in str(key).lower()
                ):
                    db_refs[key] = value

            if db_refs:
                protein_entry["database_cross_references"] = db_refs

        protein_list.append(protein_entry)

    output = {
        "metadata": {
            "total_proteins": len(proteins),
            "filter_applied": "XP identifiers only (excluding NP)",
            "filter_timestamp": datetime.now().isoformat(),
            "created_by": "filter_xp_proteins_cache.py (Enhanced with NCBI/KEGG)",
            "data_sources": ["Protein cache", "NCBI GenBank", "KEGG Database"],
            "fields_included": [
                "basic_protein_info",
                "experimental_hits",
                "full_ncbi_info",
                "full_kegg_info",
                "database_cross_references",
            ],
        },
        "proteins": protein_list,
    }

    with open(output_file, "w") as f:
        json.dump(output, f, indent=2, default=str)

    return output_file


def main():
    """Main function"""
    print("Enhanced XP Protein Filter (With NCBI/KEGG Information)")
    print("=" * 65)

    # Try different sources in order
    proteins = []

    # 1. Try loading from existing protein files
    proteins = load_proteins_from_cache()

    # 2. If empty, try CSV-based loading
    if not proteins:
        try:
            proteins = load_proteins_from_csv()
        except ImportError:
            print("Pandas not available, trying direct JSON approach...")

    # 3. Final fallback - try loading from JSON files directly
    if not proteins:
        print("Attempting direct JSON file loading...")
        try:
            with open("output/clusters_export.json", "r") as f:
                data = json.load(f)

            # Extract proteins from the JSON structure
            if isinstance(data, dict):
                if "proteins" in data:
                    for p_data in data["proteins"]:
                        protein = create_protein_from_dict(p_data)
                        if protein:
                            proteins.append(protein)
                elif "clusters" in data:
                    for cluster_data in data["clusters"].values():
                        for p_data in cluster_data.get("proteins", []):
                            protein = create_protein_from_dict(p_data)
                            if protein:
                                proteins.append(protein)

            if proteins:
                print(f"Successfully loaded {len(proteins)} proteins from JSON file")

        except Exception as e:
            raise (e)

    if not proteins:
        print("Error: Could not load any proteins from available sources")
        print("Please ensure you have:")
        print("1. A CSV file with protein data, OR")
        print(
            "2. A JSON file with protein information (filtered_II.json, clusters_filtered.json, etc.)"
        )
        sys.exit(1)

    print(f"Loaded {len(proteins)} total proteins")

    # Fetch additional NCBI and KEGG information
    proteins = fetch_additional_info(proteins)

    # Filter for XP proteins
    filtered_proteins = filter_xp_proteins(proteins)

    print(f"Found {len(filtered_proteins)} proteins with XP identifiers (no NP)")

    # Save results
    output_file = "xp_proteins_enhanced.json"
    saved_file = save_filtered_proteins(filtered_proteins, output_file)

    print(f"\nEnhanced filtering complete!")
    print(f"Results saved to: {saved_file}")
    print(f"- Total proteins filtered: {len(filtered_proteins)}")
    print(f"- Source data: {len(proteins)} total proteins")
    print(f"- Enhanced with full NCBI and KEGG information")


if __name__ == "__main__":
    main()
