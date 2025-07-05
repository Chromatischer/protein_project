import json
import numpy as np
from datetime import datetime
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from models.protein import Protein


class NumpyEncoder(json.JSONEncoder):
    """Custom JSON encoder that handles NumPy data types."""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NumpyEncoder, self).default(obj)


def export_clusters_to_json(
    clusters_dict: dict[Any, Any],
    tops: list["Protein"],
    distance_threshold: float,
    filename: str = "clusters_export.json",
) -> str:
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
            "clustering_metric": "euclidean",
        },
        "clusters": {},
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
                "xref_id": protein.xref_id,
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

    # Write to JSON file
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(export_data, f, indent=2, ensure_ascii=False, cls=NumpyEncoder)

    print(f"\n--- Clusters exported to {filename} ---")
    print(f"Total clusters: {len(clusters_dict)}")
    print(f"Total proteins: {len(tops)}")

    return filename
