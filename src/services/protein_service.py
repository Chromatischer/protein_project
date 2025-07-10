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
from utils.cache import get_cache

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
    It includes automatic retry handling, progress reporting, and caching support.

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

    cache = get_cache()
    all_protein_records = {}

    # First, try to get cached results
    cached_results = cache.get_ncbi_proteins_batch(protein_ids)
    proteins_to_fetch = []

    for protein_id in protein_ids:
        cached_record = cached_results.get(protein_id)
        if cached_record is not None:
            all_protein_records[protein_id] = cached_record
            print(f"Using cached NCBI data for {protein_id}")
        else:
            proteins_to_fetch.append(protein_id)

    if not proteins_to_fetch:
        print("All NCBI protein data found in cache!")
        return all_protein_records

    print(f"Need to fetch {len(proteins_to_fetch)} proteins from NCBI API")
    newly_fetched = {}

    # Process proteins in batches to respect API limits
    for i in range(0, len(proteins_to_fetch), BATCH_SIZE):
        batch_ids = proteins_to_fetch[i : i + BATCH_SIZE]
        print(
            f"Processing batch {i // BATCH_SIZE + 1}/{(len(proteins_to_fetch) + BATCH_SIZE - 1) // BATCH_SIZE}..."
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
                newly_fetched[protein_id] = found_record

        except Exception as e:
            print(f"Error fetching batch starting with {batch_ids[0]}: {e}")
            # Mark all proteins in this batch as failed
            for protein_id in batch_ids:
                all_protein_records[protein_id] = None
                newly_fetched[protein_id] = None

        # Add delay between batches to respect API rate limits
        if i + BATCH_SIZE < len(proteins_to_fetch):
            time.sleep(REQUEST_DELAY_SECONDS)

    # Cache the newly fetched results
    if newly_fetched:
        cache.cache_ncbi_proteins_batch(newly_fetched)
        print(f"Cached {len(newly_fetched)} new NCBI protein records")

    return all_protein_records


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
    cache = get_cache()
    result: dict = {}

    # Create identifier mapping for caching
    identifier_map = {}
    queries_to_fetch = []

    if type(proteins) is list[str]:
        # Handle list of strings
        for query in proteins:
            cached_entry = cache.get_kegg_entry(query)
            if cached_entry is not None:
                result[query] = cached_entry
                print(f"Using cached KEGG data for {query}")
            else:
                queries_to_fetch.append(query)
                identifier_map[query] = query
    else:
        # Handle list of Protein objects
        for protein in proteins:
            kegg_query = f"hsa:{protein.xref_id}"
            cached_entry = cache.get_kegg_entry(protein.sci_identifier)
            if cached_entry is not None:
                result[protein.sci_identifier] = cached_entry
                print(f"Using cached KEGG data for {protein.sci_identifier}")
            else:
                queries_to_fetch.append(kegg_query)
                identifier_map[kegg_query] = protein.sci_identifier

    if not queries_to_fetch:
        print("All KEGG data found in cache!")
        return result

    print(f"Need to fetch {len(queries_to_fetch)} entries from KEGG API")
    newly_fetched = {}

    response: str = ""
    try:
        for num in range(math.ceil(len(queries_to_fetch) / 10)):
            low = 0 + num * 10
            high = 10 + num * 10
            print(f"Fetching Kegg [{low} to {high}]")
            batch_queries = queries_to_fetch[low:high]
            response += REST.kegg_get(batch_queries).read()

    except Exception as e:
        print(e)

    if response == "":
        print("Warning: No response from KEGG API")
        return result

    # Parse response and map back to original identifiers
    i = 0
    for part in response.split("///"):
        if part.strip() != "":
            if i < len(queries_to_fetch):
                kegg_query = queries_to_fetch[i]
                identifier = identifier_map[kegg_query]
                kegg_entry = KeggEntry.from_kegg_text(part)

                result[identifier] = kegg_entry
                newly_fetched[identifier] = kegg_entry

                i += 1

    # Cache the newly fetched results
    if newly_fetched:
        cache.cache_kegg_entries_batch(newly_fetched)
        print(f"Cached {len(newly_fetched)} new KEGG entries")

    return result
