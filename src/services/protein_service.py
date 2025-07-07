"""
Protein data fetching service module.

This module provides functions to fetch protein information from NCBI databases
using Biopython's Entrez API and KEGG database using Bio.KEGG.REST API.
It includes both single protein fetching and batch processing capabilities
for efficient data retrieval with JSON-based caching.
"""

from Bio import Entrez, SeqIO
from Bio.KEGG import REST
import ssl
import time
import math
import json
import urllib.request
from pathlib import Path
from datetime import datetime
from models import Protein, KeggEntry

# Configure Entrez email for NCBI API access
Entrez.email = "dominik@hildania.de"

# Configure SSL context to handle certificate issues
ssl_context = ssl.create_default_context()
ssl_context.check_hostname = False
ssl_context.verify_mode = ssl.CERT_NONE
opener = urllib.request.build_opener(urllib.request.HTTPSHandler(context=ssl_context))
urllib.request.install_opener(opener)

# Configuration constants for batch processing
BATCH_SIZE = 50
REQUEST_DELAY_SECONDS = 0.5

# Cache configuration
CACHE_DIR = Path("cache")
CACHE_EXPIRY_DAYS = 7


def _ensure_cache_dir():
    """Ensure the cache directory exists."""
    CACHE_DIR.mkdir(exist_ok=True)


def _get_cache_filename(source: str) -> Path:
    """Get the cache filename for a given source (ncbi or kegg)."""
    return CACHE_DIR / f"{source}_protein_cache.json"


def _is_cache_valid(cache_file: Path) -> bool:
    """Check if cache file exists and is not expired."""
    if not cache_file.exists():
        return False

    try:
        cache_age = datetime.now() - datetime.fromtimestamp(cache_file.stat().st_mtime)
        return cache_age.days < CACHE_EXPIRY_DAYS
    except Exception:
        return False


def _load_cached_results(protein_ids: list, source: str) -> tuple[dict, list]:
    """
    Load cached results and return both cached data and remaining IDs to fetch.

    Args:
        protein_ids: List of protein IDs to check
        source: Cache source ('ncbi' or 'kegg')

    Returns:
        Tuple of (cached_results_dict, remaining_ids_list)
    """
    _ensure_cache_dir()
    cache_file = _get_cache_filename(source)

    cached_results = {}
    remaining_ids = protein_ids.copy()

    if _is_cache_valid(cache_file):
        try:
            with open(cache_file, "r") as f:
                cache_data = json.load(f)

            for protein_id in protein_ids:
                if protein_id in cache_data:
                    cached_results[protein_id] = cache_data[protein_id]
                    remaining_ids.remove(protein_id)

        except Exception as e:
            print(f"Error loading cache: {e}")

    return cached_results, remaining_ids


def _cache_batch_results(results: dict, source: str):
    """
    Cache batch results to JSON file.

    Args:
        results: Dictionary of protein_id -> result data
        source: Cache source ('ncbi' or 'kegg')
    """
    _ensure_cache_dir()
    cache_file = _get_cache_filename(source)

    # Load existing cache
    existing_cache = {}
    if cache_file.exists():
        try:
            with open(cache_file, "r") as f:
                existing_cache = json.load(f)
        except Exception as e:
            print(f"Error loading existing cache: {e}")

    # Update cache with new results
    # For NCBI, we need to serialize SeqRecord objects
    if source == "ncbi":
        for protein_id, record in results.items():
            if record is not None:
                # Convert SeqRecord to dict for JSON serialization
                existing_cache[protein_id] = {
                    "id": record.id,
                    "description": record.description,
                    "sequence": str(record.seq),
                    "features": len(record.features),
                    "annotations": record.annotations,
                    "_cached_at": datetime.now().isoformat(),
                }
            else:
                existing_cache[protein_id] = None
    else:
        # For KEGG, results are already JSON-serializable
        for protein_id, data in results.items():
            if data is not None:
                data["_cached_at"] = datetime.now().isoformat()
            existing_cache[protein_id] = data

    # Save updated cache
    try:
        with open(cache_file, "w") as f:
            json.dump(existing_cache, f, indent=2)
    except Exception as e:
        print(f"Error saving cache: {e}")


def fetch_protein_info(protein_id):
    """
    Fetch protein information using protein accession ID.

    This function searches for a single protein by its accession ID and
    retrieves detailed information from the NCBI protein database.

    Args:
        protein_id (str): Protein accession ID (e.g., 'NP_000047.2')

    Returns:
        Bio.SeqRecord.SeqRecord or None: Protein record if found, None if not found
        or if an error occurred during fetching

    Example:
        >>> record = fetch_protein_info('NP_000047.2')
        >>> if record:
        ...     print(record.description)
    """
    try:
        # Search for the protein ID in the protein database
        search_handle = Entrez.esearch(db="protein", term=protein_id)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if not search_results["IdList"]:
            print(f"No results found for {protein_id}")
            return None

        # Get the first ID from search results
        protein_uid = search_results["IdList"][0]

        # Fetch detailed information
        fetch_handle = Entrez.efetch(
            db="protein", id=protein_uid, rettype="gb", retmode="text"
        )
        protein_record = SeqIO.read(fetch_handle, "genbank")
        fetch_handle.close()

        return protein_record

    except Exception as e:
        print(f"Error fetching {protein_id}: {e}")
        return None


def fetch_protein_info_batch(protein_ids: list) -> dict:
    """
    Fetch protein information for a list of protein accession IDs in batches.

    This function efficiently retrieves protein information for multiple proteins
    by processing them in batches to reduce API calls and respect NCBI's rate limits.
    It includes automatic retry handling and progress reporting.

    Args:
        protein_ids (list): A list of protein accession IDs (e.g., ['NP_000047.2'])

    Returns:
        dict: A dictionary mapping each protein ID to its Bio.SeqRecord.SeqRecord
              object, or None if fetching failed for that ID

    Raises:
        TypeError: If protein_ids is not a list

    Example:
        >>> ids = ['NP_000047.2', 'NP_000048.1']
        >>> results = fetch_protein_info_batch(ids)
        >>> for protein_id, record in results.items():
        ...     if record:
        ...         print(f"{protein_id}: {record.description}")
    """
    if not isinstance(protein_ids, list):
        raise TypeError("protein_ids must be a list")
    if not protein_ids:
        return {}

    all_protein_records = {}

    # Process proteins in batches to respect API limits
    for i in range(0, len(protein_ids), BATCH_SIZE):
        batch_ids = protein_ids[i : i + BATCH_SIZE]
        print(
            f"Processing batch {i // BATCH_SIZE + 1}/{(len(protein_ids) + BATCH_SIZE - 1) // BATCH_SIZE}..."
        )

        try:
            # Fetch multiple records in a single API call
            fetch_handle = Entrez.efetch(
                db="protein",
                id=",".join(batch_ids),
                # idtype="acc" is crucial for telling Entrez these are accessions
                idtype="acc",
                rettype="gb",
                retmode="text",
            )
            # Use SeqIO.parse for multiple records; SeqIO.read is for a single record
            batch_records = SeqIO.parse(fetch_handle, "genbank")
            current_batch_results = {rec.id: rec for rec in batch_records}
            fetch_handle.close()

            # Map each requested protein ID to its corresponding record
            for protein_id in batch_ids:
                # The returned record ID (rec.id) includes the version, e.g., "NP_000047.2"
                # Match it against the original requested ID
                found_record = None
                for rec_id_in_batch, record_obj in current_batch_results.items():
                    if rec_id_in_batch.startswith(protein_id):
                        found_record = record_obj
                        break
                all_protein_records[protein_id] = found_record

        except Exception as e:
            print(f"Error fetching batch starting with {batch_ids[0]}: {e}")
            # Mark all proteins in this batch as failed
            for protein_id in batch_ids:
                all_protein_records[protein_id] = None

        # Add delay between batches to respect API rate limits
        if i + BATCH_SIZE < len(protein_ids):
            time.sleep(REQUEST_DELAY_SECONDS)

    return all_protein_records


def _is_ncbi_accession(protein_id: str) -> bool:
    """
    Check if a protein ID is an NCBI accession number.

    NCBI accession numbers typically follow patterns like:
    - NP_123456.1 (RefSeq proteins)
    - XP_123456.1 (RefSeq predicted proteins)
    - YP_123456.1 (RefSeq proteins from complete genomes)
    - WP_123456.1 (RefSeq non-redundant proteins)

    Args:
        protein_id (str): Protein identifier to check

    Returns:
        bool: True if it appears to be an NCBI accession number
    """
    import re

    # Pattern for NCBI accession numbers
    ncbi_pattern = r"^[NXYW]P_\d+\.\d+$"
    return bool(re.match(ncbi_pattern, protein_id))


def fetch_protein_info_kegg(protein_id: str) -> dict:
    if "hsa:" not in protein_id:
        raise ValueError("Human Genome has to start with hsa")

    result = None

    try:
        result = REST.kegg_get(protein_id).read()
    except Exception as e:
        print(e)

    assert result is not None

    return result


def fetch_kegg_info_batch(proteins: list[Protein] | list[str]) -> dict:
    """
    This function fetches a batch of proteins at once.
    Takes as input a list[Proteins] or a list[str] in the format: hsa:xref_id
    Returns a dict {"hsa:"|"NP_XYZ": Info}
    """
    result: dict = {}

    queries = (
        [f"hsa:{p.xref_id}" for p in proteins]
        if type(proteins) is not list[str]
        else proteins
    )

    response: str = ""
    try:
        for num in range(math.ceil(len(queries) / 10)):
            low = 0 + num * 10
            high = 10 + num * 10
            print(f"Fetching Kegg [{low} to {high}]")
            querie = queries[low:high]
            response += REST.kegg_get(querie).read()

    except Exception as e:
        print(e)

    assert response != ""

    i = 0
    for part in response.split("///"):
        if part.strip() != "":
            if type(proteins[i]) is Protein:
                index = proteins[i].sci_identifier
            else:
                index = proteins[i]

            result[index] = KeggEntry.from_kegg_text(part)

            i += 1

    return result
