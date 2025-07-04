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
from datetime import datetime
from scipy.cluster import hierarchy
from collections import defaultdict

# Import our custom modules
from models.protein import Protein
from services.protein_service import (
    fetch_protein_info_batch,
    fetch_protein_info_kegg_batch,
    search_protein_kegg
)
from utils.data_processing import extract_gene_info, grouping
from visualization.plotting import renderPlot


if __name__ == "__main__":
    df = pd.read_csv(
        "../resources/WWOP230228_PHSN_C12-2_Top100_Variant_Peptides_blinded.csv",
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
        + "{:80}".format("Extracted")
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
    
    # 3. Fetch KEGG information for all unique proteins
    print("--- Fetching KEGG information for proteins ---")
    kegg_info_map = fetch_protein_info_kegg_batch(unique_protein_ids)
    print("--- KEGG fetching complete. ---\n")
    # --- END OF REFACTORED LOGIC ---

    tops: list[Protein] = []
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

        # Get KEGG info for this protein
        kegg_info = kegg_info_map.get(row.Protein, {})
        
        prot = Protein(
            seq_id=row.SeqID,
            sci_identifier=row.Protein,
            desc=row.Description,
            aa_start=row.AA_Start,
            aa_stop=row.AA_Stop,
            peptide=row.Peptide,
            std_dev=row.std_dev,
            info=info,
            kegg_info=kegg_info,
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
            + "{:80}".format(extracted_gene)
            # Truncate for display
            + "{:100}".format(prot.info.description[:100])
            + str(list(prot.hits.values()))
        )
        tops.append(prot)

    for protein in tops:
        print(protein.sci_identifier, str(protein.hits.values()))

    grouped = grouping(tops)

    print("=== ADVANCED GROUPING ===")
    print("Count, Base Name, Proteins")
    for base_name, proteins in grouped.items():
        print(f"{'{:13}'.format(len(proteins))}, {'{:40}'.format(base_name)}")
        for protein in proteins:
            print(f"{protein['description']}")

    print("Render plot")

    row_names = []
    column_names = []
    matrix = []

    for protein in tops:
        protein: Protein = protein
        row = []
        for col_name, value in protein.hits.items():
            row.append(value)
            if col_name not in column_names:
                column_names.append(col_name)
        matrix.append(row)
        row_names.append(protein.sci_identifier)

    matrix_frame = pd.DataFrame(
        data=matrix, index=row_names, columns=column_names)

    print(matrix_frame.to_string())

    renderPlot(matrix_frame, column_names, row_names)

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

    renderPlot(matrix_frame, column_names, clustered_row_names)

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
    for cluster_id, peptides in clusters_dict.items():
        print(f"\nCluster {cluster_id}:")
        # Using np.array for cleaner printing of the list
        print(str(peptides))
        if len(peptides) > 1:
            for peptide in peptides:
                res_top = None
                for top in tops:
                    if peptide == top.sci_identifier:
                        res_top = top
                        break

                print(res_top.info)

    # --- 5. Export clusters to JSON ---
    def export_clusters_to_json(clusters_dict, tops, distance_threshold, filename="clusters_export.json"):
        """
        Export cluster information to a JSON file with complete protein data and KEGG information.
        
        Args:
            clusters_dict: Dictionary of cluster_id -> list of protein identifiers
            tops: List of Protein objects with complete data
            distance_threshold: The distance threshold used for clustering
            filename: Output JSON filename
        """
        # Create a lookup dictionary for quick protein access
        protein_lookup = {protein.sci_identifier: protein for protein in tops}
        
        # Prepare the export data structure
        export_data = {
            "metadata": {
                "export_timestamp": datetime.now().isoformat(),
                "total_proteins": len(tops),
                "total_clusters": len(clusters_dict),
                "distance_threshold": distance_threshold,
                "clustering_method": "complete",
                "clustering_metric": "euclidean"
            },
            "clusters": {}
        }
        
        # Process each cluster
        for cluster_id, protein_ids in clusters_dict.items():
            cluster_proteins = []
            
            for protein_id in protein_ids:
                protein = protein_lookup.get(protein_id)
                if protein is None:
                    continue
                    
                # Prepare protein data for JSON serialization
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
                        "sequence_length": len(protein.info.seq) if protein.info and hasattr(protein.info, 'seq') else None
                    },
                    "kegg_info": protein.kegg_info
                }
                
                cluster_proteins.append(protein_data)
            
            export_data["clusters"][str(cluster_id)] = {
                "cluster_id": cluster_id,
                "protein_count": len(cluster_proteins),
                "proteins": cluster_proteins
            }
        
        # Write to JSON file
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(export_data, f, indent=2, ensure_ascii=False)
        
        print(f"\n--- Clusters exported to {filename} ---")
        print(f"Total clusters: {len(clusters_dict)}")
        print(f"Total proteins: {len(tops)}")
        
        return filename

    # Export the clusters
    export_filename = export_clusters_to_json(clusters_dict, tops, distance_threshold)
