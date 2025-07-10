"""
Evaluation module for analyzing protein clusters and comparing protein families.
This module focuses on comparing proteins based on the KEGG BRITE classification,
specifically the "09182 Protein families: genetic information processing" field.
"""

import json
from collections import Counter
from typing import Dict, List, Set


def load_filtered_clusters(filepath: str) -> List[Dict]:
    """
    Load the filtered clusters JSON file.

    Args:
        filepath: Path to the filtered clusters JSON file

    Returns:
        List of cluster dictionaries
    """
    with open(filepath, "r") as f:
        data = json.load(f)
        # Extract clusters from the structure - clusters are stored as a dictionary
        if isinstance(data, dict) and "clusters" in data:
            return list(data["clusters"].values())
        else:
            # If it's already a list, return as is
            return data


def extract_protein_families_09182(brite_entries: List[str]) -> Set[str]:
    """
    Extract protein family information from BRITE entries that are under
    "09182 Protein families: genetic information processing".

    Args:
        brite_entries: List of BRITE classification entries

    Returns:
        Set of protein family classifications under 09182
    """
    protein_families = set()
    in_09182_section = False

    for entry in brite_entries:
        # Check if we're entering the 09182 section
        if "Protein families" in entry:
            in_09182_section = True
            protein_families.add(entry.split(":")[1])
            continue

        # Check if we're leaving the 09182 section (next main category)
        if in_09182_section and entry.startswith("0"):
            in_09182_section = False
            continue

        # If we're in the 09182 section, collect the protein family info
        if in_09182_section:
            # Skip empty entries and general hierarchy markers
            if entry.strip() and not entry.startswith("KEGG"):
                protein_families.add(entry.strip())

    return protein_families


def analyze_cluster_protein_families(cluster: Dict) -> Dict[str, Set[str]]:
    """
    Analyze protein families for all proteins in a cluster.

    Args:
        cluster: Dictionary containing cluster information

    Returns:
        Dictionary mapping protein IDs to their 09182 protein families
    """
    protein_families_map = {}

    for protein in cluster.get("proteins", []):
        protein_id = protein.get("sci_identifier", "Unknown")
        kegg_info = protein.get("kegg_info", {})
        brite_entries = kegg_info.get("brite", [])

        families = extract_protein_families_09182(brite_entries)
        protein_families_map[protein_id] = families

    return protein_families_map


def compare_protein_families_in_cluster(cluster: Dict) -> Dict:
    """
    Compare protein families within a single cluster and provide statistics.

    Args:
        cluster: Dictionary containing cluster information

    Returns:
        Dictionary with comparison results and statistics
    """
    protein_families_map = analyze_cluster_protein_families(cluster)

    # Count occurrences of each protein family
    family_counter = Counter()
    for families in protein_families_map.values():
        family_counter.update(families)

    # Calculate statistics
    total_proteins = len(protein_families_map)
    proteins_with_families = sum(
        1 for families in protein_families_map.values() if families
    )
    unique_families = set()
    for families in protein_families_map.values():
        unique_families.update(families)

    # Find common families (present in multiple proteins)
    common_families = {
        family: count for family, count in family_counter.items() if count > 1
    }

    # Calculate homogeneity score (percentage of proteins sharing at least one family)
    if total_proteins > 1:
        shared_family_count = 0
        for i, (protein1, families1) in enumerate(protein_families_map.items()):
            for j, (protein2, families2) in enumerate(protein_families_map.items()):
                if i < j and families1.intersection(families2):
                    shared_family_count += 1

        max_pairs = (total_proteins * (total_proteins - 1)) // 2
        homogeneity_score = shared_family_count / max_pairs if max_pairs > 0 else 0
    else:
        homogeneity_score = 1.0 if proteins_with_families > 0 else 0.0

    return {
        "cluster_id": cluster.get("cluster_id", "Unknown"),
        "total_proteins": total_proteins,
        "proteins_with_families": proteins_with_families,
        "unique_families_count": len(unique_families),
        "unique_families": list(unique_families),
        "common_families": common_families,
        "family_distribution": dict(family_counter),
        "homogeneity_score": homogeneity_score,
        "protein_families_map": {k: list(v) for k, v in protein_families_map.items()},
    }


def evaluate_all_clusters(filepath: str) -> List[Dict]:
    """
    Evaluate protein family comparisons for all clusters.

    Args:
        filepath: Path to the filtered clusters JSON file

    Returns:
        List of evaluation results for each cluster
    """
    clusters = load_filtered_clusters(filepath)
    results = []
    skipped = []

    for cluster in clusters:
        if int(cluster.get("protein_count")) < 2:
            skipped.append(cluster.get("cluster_id"))
            continue
        result = compare_protein_families_in_cluster(cluster)
        results.append(result)

    return results, skipped


def generate_summary_statistics(evaluation_results: List[Dict]) -> Dict:
    """
    Generate summary statistics across all clusters.

    Args:
        evaluation_results: List of evaluation results from evaluate_all_clusters

    Returns:
        Dictionary with summary statistics
    """
    total_clusters = len(evaluation_results)
    total_proteins = sum(result["total_proteins"]
                         for result in evaluation_results)
    total_proteins_with_families = sum(
        result["proteins_with_families"] for result in evaluation_results
    )

    # Homogeneity statistics
    homogeneity_scores = [result["homogeneity_score"]
                          for result in evaluation_results]
    avg_homogeneity = (
        sum(homogeneity_scores) /
        len(homogeneity_scores) if homogeneity_scores else 0
    )

    # Family distribution across all clusters
    all_families = Counter()
    for result in evaluation_results:
        all_families.update(result["family_distribution"].keys())

    # Clusters with high homogeneity (>= 0.5)
    high_homogeneity_clusters = sum(
        1 for score in homogeneity_scores if score >= 0.5)

    return {
        "total_clusters": total_clusters,
        "total_proteins": total_proteins,
        "total_proteins_with_families": total_proteins_with_families,
        "coverage_percentage": (total_proteins_with_families / total_proteins * 100)
        if total_proteins > 0
        else 0,
        "average_homogeneity_score": avg_homogeneity,
        "high_homogeneity_clusters": high_homogeneity_clusters,
        "high_homogeneity_percentage": (
            high_homogeneity_clusters / total_clusters * 100
        )
        if total_clusters > 0
        else 0,
        "unique_families_across_all_clusters": len(all_families),
        "most_common_families": dict(all_families.most_common(10)),
    }


def print_evaluation_report(
    evaluation_results: List[Dict], summary_stats: Dict, skipped: list[int]
):
    """
    Print a formatted evaluation report.

    Args:
        evaluation_results: List of evaluation results from evaluate_all_clusters
        summary_stats: Summary statistics from generate_summary_statistics
    """
    print("=" * 80)
    print("PROTEIN FAMILY EVALUATION REPORT")
    print("=" * 80)

    print("\nSUMMARY STATISTICS:")
    print(f"Total clusters analyzed: {summary_stats['total_clusters']}")
    print(f"Total proteins: {summary_stats['total_proteins']}")
    print(
        f"Proteins with family annotations: {summary_stats['total_proteins_with_families']}"
    )
    print(f"Coverage: {summary_stats['coverage_percentage']:.1f}%")
    print(
        f"Average homogeneity score: {summary_stats['average_homogeneity_score']:.3f}"
    )
    print(
        f"Clusters with high homogeneity (≥0.5): {summary_stats['high_homogeneity_clusters']} ({summary_stats['high_homogeneity_percentage']:.1f}%)"
    )
    print(
        f"Unique protein families found: {summary_stats['unique_families_across_all_clusters']}"
    )

    print("\nMOST COMMON PROTEIN FAMILIES:")
    for family, count in summary_stats["most_common_families"].items():
        print(f"  {family}: {count} clusters")

    print("\nDETAILED CLUSTER ANALYSIS:")
    print("-" * 80)

    for result in evaluation_results:
        print(f"\nCluster ID: {result['cluster_id']}")
        print(
            f"  Proteins: {result['total_proteins']} (with families: {result['proteins_with_families']})"
        )
        print(f"  Unique families: {result['unique_families_count']}")
        print(f"  Homogeneity score: {result['homogeneity_score']:.3f}")

        if result["common_families"]:
            print("  Common families:")
            for family, count in result["common_families"].items():
                print(f"    {family}: {count} proteins")
        else:
            print("  No common families found")

    print("\nSKIPPED:")
    print(f"{len(skipped)} Entries: [{str(skipped)}]")


if __name__ == "__main__":
    # Example usage
    filepath = "../clusters_filtered.json"

    print("Loading and evaluating clusters...")
    evaluation_results, skipped = evaluate_all_clusters(filepath)
    summary_stats = generate_summary_statistics(evaluation_results)

    print_evaluation_report(evaluation_results, summary_stats, skipped)
