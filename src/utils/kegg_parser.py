"""
KEGG data parser utilities.

This module provides functions to parse KEGG database entries into structured data.
"""

from models.kegg_model import KeggEntry


def parse_kegg_entry(text_entry: str) -> KeggEntry:
    """
    Parse a KEGG text entry into a KeggEntry object.
    
    Args:
        text_entry (str): Raw KEGG text entry from KEGG API
        
    Returns:
        KeggEntry: Parsed KEGG entry object
    """
    return KeggEntry.from_kegg_text(text_entry)

