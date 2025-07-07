"""
Models package for protein data structures.
"""

from .protein import Protein
from .kegg_model import KeggEntry

__all__ = ["Protein", "KeggEntry"]
