"""
Services package for data fetching functionality.
"""

from .protein_service import fetch_protein_info, fetch_protein_info_batch

__all__ = ['fetch_protein_info', 'fetch_protein_info_batch']

