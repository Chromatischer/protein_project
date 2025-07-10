"""
Caching module for protein data.

This module provides caching functionality for NCBI and KEGG protein data
to avoid redundant API calls and improve performance. It supports both
SeqRecord objects from NCBI and KeggEntry objects from KEGG with automatic
serialization and deserialization.
"""

import json
import pickle
import hashlib
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, Any, Optional, Union, List
from dataclasses import asdict

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from io import StringIO
from models.kegg_model import KeggEntry


class ProteinCache:
    """
    A caching system for protein data from NCBI and KEGG databases.
    
    This class provides methods to cache and retrieve protein information,
    with automatic expiration and efficient serialization for both NCBI
    SeqRecord objects and KEGG KeggEntry objects.
    """
    
    def __init__(self, cache_dir: str = "cache", cache_expiry_days: int = 30):
        """
        Initialize the protein cache.
        
        Args:
            cache_dir (str): Directory to store cache files
            cache_expiry_days (int): Number of days after which cache expires
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        
        self.ncbi_cache_dir = self.cache_dir / "ncbi"
        self.kegg_cache_dir = self.cache_dir / "kegg"
        self.ncbi_cache_dir.mkdir(exist_ok=True)
        self.kegg_cache_dir.mkdir(exist_ok=True)
        
        self.cache_expiry = timedelta(days=cache_expiry_days)
        
        # Metadata files to track cache info
        self.ncbi_metadata_file = self.ncbi_cache_dir / "metadata.json"
        self.kegg_metadata_file = self.kegg_cache_dir / "metadata.json"
        
        # Load existing metadata
        self.ncbi_metadata = self._load_metadata(self.ncbi_metadata_file)
        self.kegg_metadata = self._load_metadata(self.kegg_metadata_file)
    
    def _load_metadata(self, metadata_file: Path) -> Dict[str, Any]:
        """Load metadata from file or create empty metadata."""
        if metadata_file.exists():
            try:
                with open(metadata_file, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError):
                return {}
        return {}
    
    def _save_metadata(self, metadata: Dict[str, Any], metadata_file: Path):
        """Save metadata to file."""
        try:
            with open(metadata_file, 'w') as f:
                json.dump(metadata, f, indent=2, default=str)
        except IOError as e:
            print(f"Warning: Could not save metadata to {metadata_file}: {e}")
    
    def _get_cache_key(self, identifier: str) -> str:
        """Generate a cache key from an identifier."""
        return hashlib.md5(identifier.encode()).hexdigest()
    
    def _is_cache_valid(self, cache_time: str) -> bool:
        """Check if cache entry is still valid based on expiry time."""
        try:
            cached_datetime = datetime.fromisoformat(cache_time)
            return datetime.now() - cached_datetime < self.cache_expiry
        except (ValueError, TypeError):
            return False
    
    def _serialize_seqrecord(self, seq_record: SeqRecord) -> str:
        """Serialize a SeqRecord object to GenBank format string."""
        output = StringIO()
        SeqIO.write(seq_record, output, "genbank")
        return output.getvalue()
    
    def _deserialize_seqrecord(self, data: str) -> SeqRecord:
        """Deserialize a GenBank format string to SeqRecord object."""
        input_stream = StringIO(data)
        return SeqIO.read(input_stream, "genbank")
    
    def _serialize_kegg_entry(self, kegg_entry: KeggEntry) -> Dict[str, Any]:
        """Serialize a KeggEntry object to dictionary."""
        return asdict(kegg_entry)
    
    def _deserialize_kegg_entry(self, data: Dict[str, Any]) -> KeggEntry:
        """Deserialize a dictionary to KeggEntry object."""
        return KeggEntry(**data)
    
    # NCBI Cache Methods
    
    def cache_ncbi_protein(self, protein_id: str, seq_record: SeqRecord):
        """
        Cache an NCBI protein SeqRecord.
        
        Args:
            protein_id (str): Protein accession ID
            seq_record (SeqRecord): The protein record to cache
        """
        if seq_record is None:
            return
        
        cache_key = self._get_cache_key(protein_id)
        cache_file = self.ncbi_cache_dir / f"{cache_key}.gb"
        
        try:
            # Serialize and save the SeqRecord
            serialized_data = self._serialize_seqrecord(seq_record)
            with open(cache_file, 'w') as f:
                f.write(serialized_data)
            
            # Update metadata
            self.ncbi_metadata[protein_id] = {
                'cache_key': cache_key,
                'cached_time': datetime.now().isoformat(),
                'file_path': str(cache_file)
            }
            self._save_metadata(self.ncbi_metadata, self.ncbi_metadata_file)
            
        except Exception as e:
            print(f"Warning: Could not cache NCBI protein {protein_id}: {e}")
    
    def get_ncbi_protein(self, protein_id: str) -> Optional[SeqRecord]:
        """
        Retrieve a cached NCBI protein SeqRecord.
        
        Args:
            protein_id (str): Protein accession ID
            
        Returns:
            SeqRecord or None: The cached protein record if found and valid
        """
        if protein_id not in self.ncbi_metadata:
            return None
        
        metadata = self.ncbi_metadata[protein_id]
        
        # Check if cache is still valid
        if not self._is_cache_valid(metadata['cached_time']):
            self._remove_ncbi_cache_entry(protein_id)
            return None
        
        try:
            cache_file = Path(metadata['file_path'])
            if cache_file.exists():
                with open(cache_file, 'r') as f:
                    serialized_data = f.read()
                return self._deserialize_seqrecord(serialized_data)
        except Exception as e:
            print(f"Warning: Could not load cached NCBI protein {protein_id}: {e}")
            self._remove_ncbi_cache_entry(protein_id)
        
        return None
    
    def cache_ncbi_proteins_batch(self, protein_records: Dict[str, SeqRecord]):
        """
        Cache multiple NCBI protein records.
        
        Args:
            protein_records (Dict[str, SeqRecord]): Dictionary mapping protein IDs to records
        """
        for protein_id, seq_record in protein_records.items():
            self.cache_ncbi_protein(protein_id, seq_record)
    
    def get_ncbi_proteins_batch(self, protein_ids: List[str]) -> Dict[str, Optional[SeqRecord]]:
        """
        Retrieve multiple cached NCBI protein records.
        
        Args:
            protein_ids (List[str]): List of protein accession IDs
            
        Returns:
            Dict[str, Optional[SeqRecord]]: Dictionary mapping protein IDs to cached records
        """
        results = {}
        for protein_id in protein_ids:
            results[protein_id] = self.get_ncbi_protein(protein_id)
        return results
    
    def _remove_ncbi_cache_entry(self, protein_id: str):
        """Remove a specific NCBI cache entry."""
        if protein_id in self.ncbi_metadata:
            metadata = self.ncbi_metadata[protein_id]
            cache_file = Path(metadata['file_path'])
            if cache_file.exists():
                cache_file.unlink()
            del self.ncbi_metadata[protein_id]
            self._save_metadata(self.ncbi_metadata, self.ncbi_metadata_file)
    
    # KEGG Cache Methods
    
    def cache_kegg_entry(self, identifier: str, kegg_entry: KeggEntry):
        """
        Cache a KEGG entry.
        
        Args:
            identifier (str): Protein identifier or KEGG ID
            kegg_entry (KeggEntry): The KEGG entry to cache
        """
        if kegg_entry is None:
            return
        
        cache_key = self._get_cache_key(identifier)
        cache_file = self.kegg_cache_dir / f"{cache_key}.json"
        
        try:
            # Serialize and save the KeggEntry
            serialized_data = self._serialize_kegg_entry(kegg_entry)
            with open(cache_file, 'w') as f:
                json.dump(serialized_data, f, indent=2)
            
            # Update metadata
            self.kegg_metadata[identifier] = {
                'cache_key': cache_key,
                'cached_time': datetime.now().isoformat(),
                'file_path': str(cache_file)
            }
            self._save_metadata(self.kegg_metadata, self.kegg_metadata_file)
            
        except Exception as e:
            print(f"Warning: Could not cache KEGG entry {identifier}: {e}")
    
    def get_kegg_entry(self, identifier: str) -> Optional[KeggEntry]:
        """
        Retrieve a cached KEGG entry.
        
        Args:
            identifier (str): Protein identifier or KEGG ID
            
        Returns:
            KeggEntry or None: The cached KEGG entry if found and valid
        """
        if identifier not in self.kegg_metadata:
            return None
        
        metadata = self.kegg_metadata[identifier]
        
        # Check if cache is still valid
        if not self._is_cache_valid(metadata['cached_time']):
            self._remove_kegg_cache_entry(identifier)
            return None
        
        try:
            cache_file = Path(metadata['file_path'])
            if cache_file.exists():
                with open(cache_file, 'r') as f:
                    serialized_data = json.load(f)
                return self._deserialize_kegg_entry(serialized_data)
        except Exception as e:
            print(f"Warning: Could not load cached KEGG entry {identifier}: {e}")
            self._remove_kegg_cache_entry(identifier)
        
        return None
    
    def cache_kegg_entries_batch(self, kegg_entries: Dict[str, KeggEntry]):
        """
        Cache multiple KEGG entries.
        
        Args:
            kegg_entries (Dict[str, KeggEntry]): Dictionary mapping identifiers to KEGG entries
        """
        for identifier, kegg_entry in kegg_entries.items():
            self.cache_kegg_entry(identifier, kegg_entry)
    
    def get_kegg_entries_batch(self, identifiers: List[str]) -> Dict[str, Optional[KeggEntry]]:
        """
        Retrieve multiple cached KEGG entries.
        
        Args:
            identifiers (List[str]): List of protein identifiers or KEGG IDs
            
        Returns:
            Dict[str, Optional[KeggEntry]]: Dictionary mapping identifiers to cached entries
        """
        results = {}
        for identifier in identifiers:
            results[identifier] = self.get_kegg_entry(identifier)
        return results
    
    def _remove_kegg_cache_entry(self, identifier: str):
        """Remove a specific KEGG cache entry."""
        if identifier in self.kegg_metadata:
            metadata = self.kegg_metadata[identifier]
            cache_file = Path(metadata['file_path'])
            if cache_file.exists():
                cache_file.unlink()
            del self.kegg_metadata[identifier]
            self._save_metadata(self.kegg_metadata, self.kegg_metadata_file)
    
    # Cache Management Methods
    
    def clear_expired_cache(self):
        """Remove all expired cache entries."""
        # Clear expired NCBI entries
        expired_ncbi = []
        for protein_id, metadata in self.ncbi_metadata.items():
            if not self._is_cache_valid(metadata['cached_time']):
                expired_ncbi.append(protein_id)
        
        for protein_id in expired_ncbi:
            self._remove_ncbi_cache_entry(protein_id)
        
        # Clear expired KEGG entries
        expired_kegg = []
        for identifier, metadata in self.kegg_metadata.items():
            if not self._is_cache_valid(metadata['cached_time']):
                expired_kegg.append(identifier)
        
        for identifier in expired_kegg:
            self._remove_kegg_cache_entry(identifier)
        
        print(f"Removed {len(expired_ncbi)} expired NCBI entries and {len(expired_kegg)} expired KEGG entries")
    
    def clear_all_cache(self):
        """Remove all cache entries."""
        # Clear NCBI cache
        for cache_file in self.ncbi_cache_dir.glob("*.gb"):
            cache_file.unlink()
        self.ncbi_metadata.clear()
        self._save_metadata(self.ncbi_metadata, self.ncbi_metadata_file)
        
        # Clear KEGG cache
        for cache_file in self.kegg_cache_dir.glob("*.json"):
            cache_file.unlink()
        self.kegg_metadata.clear()
        self._save_metadata(self.kegg_metadata, self.kegg_metadata_file)
        
        print("All cache entries cleared")
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """Get statistics about the cache."""
        return {
            'ncbi_entries': len(self.ncbi_metadata),
            'kegg_entries': len(self.kegg_metadata),
            'cache_dir': str(self.cache_dir),
            'expiry_days': self.cache_expiry.days
        }


# Global cache instance
_protein_cache = None

def get_cache() -> ProteinCache:
    """Get the global protein cache instance."""
    global _protein_cache
    if _protein_cache is None:
        _protein_cache = ProteinCache()
    return _protein_cache

