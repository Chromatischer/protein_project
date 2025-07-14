"""
Models package for protein data structures.
"""

from .protein import Protein
from .kegg_model import KeggEntry
from .ncbi_model import NCBIEntry

__all__ = ["Protein", "KeggEntry", "NCBIEntry"]
