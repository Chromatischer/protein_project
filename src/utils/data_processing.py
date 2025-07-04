"""
Data processing utilities module.

This module provides functions for processing protein data, including
gene name extraction, protein grouping, and other data manipulation tasks.
"""

import re
from collections import defaultdict


def extract_gene_info(description):
    """
    Extract more detailed gene information for better grouping.
    
    This function cleans protein descriptions by removing organism information,
    prediction prefixes, and isoform/variant details to extract the base gene name.
    This enables better grouping of related proteins.
    
    Args:
        description (str): Raw protein description from database
        
    Returns:
        str: Cleaned base gene name suitable for grouping
        
    Example:
        >>> desc = "PREDICTED: hemoglobin subunit alpha isoform X1 [Homo sapiens]"
        >>> extract_gene_info(desc)
        'hemoglobin subunit alpha'
    """
    # Remove organism info (content within square brackets)
    desc_clean = re.sub(r"\[.*?\]", "", description).strip()

    # Remove PREDICTED: prefix
    desc_clean = re.sub(r"^PREDICTED:\s*", "", desc_clean, flags=re.IGNORECASE)

    # Extract potential gene name patterns
    gene_patterns = [
        r"(\w+)\s+isoform",  # Gene name before isoform
        r"(\w+)\s+variant",  # Gene name before variant
        r"^([A-Z0-9]+)\s",  # Capital gene name at start
    ]

    # Remove isoform information to get base gene name
    isoform_patterns = [
        r"\s+isoform\s+[A-Za-z0-9]+.*$",
        r"\s+variant\s+[A-Za-z0-9]+.*$",
        r"\s+form\s+[A-Za-z0-9]+.*$",
        r"\s+type\s+[A-Za-z0-9]+.*$",
    ]

    base_name = desc_clean
    for pattern in isoform_patterns:
        base_name = re.sub(pattern, "", base_name, flags=re.IGNORECASE)

    return base_name.strip()


def grouping(protein_data):
    """
    Advanced grouping with better gene name extraction.
    
    This function groups proteins by their extracted gene names, allowing
    for better organization of related proteins (e.g., different isoforms
    of the same gene).
    
    Args:
        protein_data (list): List of Protein objects to group
        
    Returns:
        dict: Dictionary where keys are base gene names and values are lists
              of dictionaries containing protein information
              
    Example:
        >>> proteins = [protein1, protein2, protein3]
        >>> grouped = grouping(proteins)
        >>> for gene_name, protein_list in grouped.items():
        ...     print(f"{gene_name}: {len(protein_list)} proteins")
    """
    grouped = defaultdict(list)

    for protein in protein_data:
        info = protein.info.description
        base_name = extract_gene_info(info)
        grouped[base_name].append(
            {
                "description": info,
                "original_desc": protein.info,
            }
        )

    return dict(grouped)

