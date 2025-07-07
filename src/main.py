"""
Main module for protein data analysis.

This module orchestrates the analysis of protein data including:
- Loading and parsing CSV data
- Fetching protein information from NCBI
- Processing and grouping proteins
- Visualizing results with heatmaps
- Performing hierarchical clustering
"""

import pandas as pd
import json
import numpy as np
from datetime import datetime
from scipy.cluster import hierarchy
from collections import defaultdict

# Import our custom modules
from models.protein import Protein
from services.protein_service import (
    fetch_protein_info_batch,
    fetch_protein_info_kegg_batch,
)
from utils.data_processing import extract_gene_info, grouping
from visualization.plotting import renderPlot
from utils.json_utils import export_clusters_to_json


if __name__ == "__main__":
    df = pd.read_csv(
        "resources/WWOP230228_PHSN_C12-2_Top100_Variant_Peptides_blinded.csv",
        dtype={
            "id": str,
            "protein": str,
            "description": str,
            "aa_start": int,
            "aa_stop": int,
            "peptide": str,
            "K004A1": str,
            "K004A2": str,
            "K225A1": str,
            "K225A2": str,
            "K003B1": str,
            "K003B2": str,
            "K237A1": str,
            "K237A2": str,
            "K022B1": str,
            "K022B2": str,
            "K302A1": str,
            "K302A2": str,
            "std_dev": float,
        },
    )
    print(df.to_string())
    print(
        "{:18}".format("Protein")
        + "{:9}".format("X-RefID")
        + "{:70}".format("Extracted")
        + "{:100}".format("From")
        + "Values"
    )

    # --- REFACTORED LOGIC ---
    # 1. Get a list of unique protein IDs from the DataFrame to avoid redundant fetches.
    unique_protein_ids = df["Protein"].unique().tolist()
    print(f"Found {len(unique_protein_ids)} unique protein IDs to fetch.")

    # 2. Fetch all protein information in one batched operation.
    # This returns a dictionary: {'protein_id': SeqRecord, ...}
    protein_info_map = fetch_protein_info_batch(unique_protein_ids)
    print("\n--- Protein fetching complete. Processing DataFrame rows. ---\n")

    protein_list: list[Protein] = []
    # Collect all xref_ids for KEGG batch fetching
    xref_ids_for_kegg = []

    for row in df.itertuples():
        # Get the pre-fetched info for the protein in the current row.
        # Use .get() to safely handle cases where a protein was not found.
        info = protein_info_map.get(row.Protein)

        # Skip rows where protein info could not be fetched.
        if info is None:
            print(
                f"Skipping row for protein '{row.Protein}' as its info could not be fetched."
            )
            continue

        cds: dict = None
        xref_id: int = None

        for feature in info.features:
            if feature.type == "CDS":
                cds = feature
                break

        print(cds.qualifiers.get("db_xref"))
        for xref in cds.qualifiers.get("db_xref"):
            if "GeneID" in xref:
                xref_id = int(xref.split("GeneID:")[1])
                break

        assert xref_id is not None, "X-RefID is None, should be int"
        xref_ids_for_kegg.append(xref_id)

        prot = Protein(
            seq_id=row.SeqID,
            sci_identifier=row.Protein,
            desc=row.Description,
            aa_start=row.AA_Start,
            aa_stop=row.AA_Stop,
            peptide=row.Peptide,
            std_dev=row.std_dev,
            info=info,
            xref_id=xref_id,
            kegg_info=None,  # Will be populated after KEGG batch fetch
            K004A1=row.K004A1,
            K004A2=row.K004A2,
            K225A1=row.K225A1,
            K225A2=row.K225A2,
            K003B1=row.K003B1,
            K003B2=row.K003B2,
            K237A1=row.K237A1,
            K237A2=row.K237A2,
            K022B1=row.K022B1,
            K022B2=row.K022B2,
            K302A1=row.K302A1,
            K302A2=row.K302A2,
        )

        extracted_gene = extract_gene_info(prot.info.description)
        print(
            "{:18}".format(prot.sci_identifier)
            + "{:9}".format(str(prot.xref_id))
            + "{:70}".format(extracted_gene)
            # Truncate for display
            + "{:100}".format(prot.info.description[:100])
            + str(list(prot.hits.values()))
        )
        protein_list.append(prot)

    # 3. Fetch KEGG information in a batch using the collected xref_ids.
    print(f"Found {len(xref_ids_for_kegg)} unique xref_ids for KEGG fetching.")
    # Pass the list of Protein objects directly to the batch function
    kegg_info_map = fetch_protein_info_kegg_batch(protein_list)
    print("\n--- KEGG fetching complete. Populating KEGG info for proteins. ---\n")

    # 4. Populate kegg_info for each Protein object.
    for protein in protein_list:
        # The key in kegg_info_map is the 'hsa:xref_id' string, not just the xref_id
        kegg_data = kegg_info_map.get(f"hsa:{protein.xref_id}")
        if kegg_data:
            protein.kegg_info = kegg_data
        else:
            print(
                f"No KEGG info found for protein {protein.sci_identifier} (X-RefID: {protein.xref_id})."
            )

    print("\n--- KEGG info population complete. ---\n")

    print("Render plot")

    row_names = []
    column_names = []
    matrix = []

    for protein_sciID_from_cluster in protein_list:
        protein_sciID_from_cluster: Protein = protein_sciID_from_cluster
        row = []
        for col_name, value in protein_sciID_from_cluster.hits.items():
            row.append(value)
            if col_name not in column_names:
                column_names.append(col_name)
        matrix.append(row)
        row_names.append(protein_sciID_from_cluster.sci_identifier)

    matrix_frame = pd.DataFrame(data=matrix, index=row_names, columns=column_names)

    print(matrix_frame.to_string())

    # renderPlot(matrix_frame, column_names, row_names)

    print("--- Original DataFrame Head ---")
    print(matrix_frame.head())

    # --- 2. Perform Hierarchical Clustering on ROWS ONLY ---

    # The linkage function performs the clustering
    row_linkage = hierarchy.linkage(
        matrix_frame.values, method="complete", metric="euclidean"
    )

    # The dendrogram function calculates the optimal leaf ordering for rows
    row_dendrogram = hierarchy.dendrogram(row_linkage, no_plot=True)

    # The 'leaves' key gives the indices of the rows in their new, clustered order
    clustered_row_indices = row_dendrogram["leaves"]

    # Get the actual row names in this new order
    clustered_row_names = matrix_frame.index[clustered_row_indices]

    # --- 3. Create the New DataFrame with Clustered Rows ---
    # Use .reindex() with only the 'index' argument to reorder rows
    # and leave the columns in their original order.
    clustered_rows_df = matrix_frame.reindex(index=clustered_row_names)

    print("\n--- DataFrame with Clustered Rows ---")
    print("Rows are reordered based on similarity; columns are unchanged.")
    print(clustered_rows_df.head())

    # renderPlot(matrix_frame, column_names, clustered_row_names)

    # --- 2. Get Cluster Labels ---
    # Define a distance threshold to cut the dendrogram. You may need to adjust this.
    # A smaller 't' will result in more, smaller clusters.
    distance_threshold = 50000
    cluster_labels = hierarchy.fcluster(
        row_linkage, t=distance_threshold, criterion="distance"
    )
    # `cluster_labels` is a NumPy array like: array([2, 1, 1, 3, 2, ...])

    # --- 3. Group Row Names into a Dictionary ---
    # defaultdict(list) creates a dictionary where each new key automatically
    # gets an empty list as its value.
    clusters_dict = defaultdict(list)

    # zip pairs each row name with its corresponding cluster label
    for row_name, label in zip(matrix_frame.index, cluster_labels):
        clusters_dict[label].append(row_name)

    # --- 4. View the Result ---
    print("--- Clusters stored in a dictionary ---")
    for cluster_id, clustered_proteins in clusters_dict.items():
        print(f"\nCluster {cluster_id}:")
        # Using np.array for cleaner printing of the list
        print(str(clustered_proteins))
        if len(clustered_proteins) > 1:
            for protein_sciID_from_cluster in clustered_proteins:
                protein = next(
                    (
                        p
                        for p in protein_list
                        if p.sci_identifier == protein_sciID_from_cluster
                    ),
                    None,
                )

                print(protein.sci_identifier)
                print(protein.xref_id)

    # Export the clusters
    export_filename = export_clusters_to_json(
        clusters_dict, protein_list, distance_threshold
    )
