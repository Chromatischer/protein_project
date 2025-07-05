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
import urllib.request
import time
import json
from pathlib import Path
from datetime import datetime


# Configure Entrez email for NCBI API access
Entrez.email = "dominik@hildania.de"

# Configure SSL context to handle certificate issues
ssl_context = ssl.create_default_context()
ssl_context.check_hostname = False
ssl_context.verify_mode = ssl.CERT_NONE
opener = urllib.request.build_opener(
    urllib.request.HTTPSHandler(context=ssl_context))
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
        batch_ids = protein_ids[i: i + BATCH_SIZE]
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


def convert_refseq_to_uniprot_batch(refseq_ids: list) -> dict:
    """
    Convert a batch of RefSeq accession numbers to UniProt IDs using UniProt ID Mapping API.

    Args:
        refseq_ids: A list of RefSeq accession numbers.

    Returns:
        A dictionary mapping RefSeq IDs to UniProt IDs.
    """
    import requests
    import time
    import json

    if not refseq_ids:
        return {}

    # UniProt ID Mapping API endpoints
    jobId: int = None
    BASE_URL = "https://rest.uniprot.org/"
    RUN_URL = f"{BASE_URL}/idmapping/run"
    RESULTS_URL = f"{BASE_URL}/idmapping/uniref/results/{jobId}"

    # Prepare the request payload
    payload = {
        "from": "RefSeq_NT",  # Adjusted to match UniProt's accepted nomenclature
        "to": "UniProtKB",
        "ids": ",".join(list(set(refseq_ids))),
    }

    try:
        # Submit ID mapping job
        headers = {"Content-Type": "application/json",
                   "Accept": "application/json"}
        response = requests.post(RUN_URL, json=payload, headers=headers)

        # Detailed error logging
        print("Request Payload:", json.dumps(payload, indent=2))
        print(f"API Response Status: {response.status_code}")
        print(f"API Response Content: {response.text}")

        # Raise an exception for HTTP errors
        response.raise_for_status()

        # Extract job ID
        response_json = response.json()
        job_id = response_json.get("jobId")
        if not job_id:
            print("Failed to obtain job ID from UniProt API")
            print("Full response:", response_json)
            return {}

        # Poll for job completion
        max_attempts = 10
        poll_interval = 2  # seconds
        for attempt in range(max_attempts):
            print(
                f"Polling job {job_id} (Attempt {attempt + 1}/{max_attempts})")

            results_url = f"{BASE_URL}/results/{job_id}"
            results_response = requests.get(results_url)
            results_response.raise_for_status()

            results_data = results_response.json()

            # Check if job is complete
            if results_data.get("results"):
                # Create mapping dictionary
                mapping = {}
                for result in results_data["results"]:
                    # Extract source and target IDs
                    source_id = result.get("from")
                    target_id = result.get("to", {}).get("primaryAccession")
                    if source_id and target_id:
                        mapping[source_id] = target_id

                print(f"Successfully mapped {len(mapping)} IDs")
                return mapping

            # If not complete, wait and retry
            time.sleep(poll_interval)

        print(f"ID mapping job {job_id} did not complete within expected time")
        return {}

    except requests.RequestException as e:
        print(f"Error in UniProt ID mapping: {e}")
        # Print any additional error details
        print(
            "Error details:",
            e.response.text if hasattr(
                e, "response") else "No additional details",
        )
        return {}


def fetch_protein_info_kegg(protein_id: str) -> dict:
    """
    Fetch protein information from KEGG database using a protein ID.

    This function queries the KEGG database for protein information using
    the Bio.KEGG.REST API. It can handle various KEGG protein identifiers
    including UniProt IDs and KEGG gene/protein IDs.

    Note: NCBI accession numbers (like NP_123456.1, XP_123456.1) are not
    compatible with KEGG API and will be skipped to avoid errors.

    Args:
        protein_id (str): Protein identifier (e.g., 'hsa:7157' for human p53,
                         or UniProt ID like 'P04637')

    Returns:
        dict: A dictionary containing protein information from KEGG including:
              - 'entry': KEGG entry information
              - 'sequence': protein sequence if available
              - 'pathway': associated pathways
              - 'error': error message if query failed
              - 'skipped': True if ID format is incompatible with KEGG

    Example:
        >>> info = fetch_protein_info_kegg('hsa:7157')
        >>> if 'error' not in info:
        ...     print(info['entry'])
    """
    try:
        result = {}

        # Try to get protein entry information
        try:
            entry_handle = REST.kegg_get(protein_id)
            entry_data = entry_handle.read()
            entry_handle.close()
            result["entry"] = (
                entry_data.decode("utf-8")
                if isinstance(entry_data, bytes)
                else entry_data
            )
        except Exception as e:
            result["entry_error"] = f"Could not fetch entry for {protein_id}: {e}"

        # Try to get protein sequence
        try:
            seq_handle = REST.kegg_get(protein_id, "aaseq")
            seq_data = seq_handle.read()
            seq_handle.close()
            result["sequence"] = (
                seq_data.decode("utf-8") if isinstance(seq_data,
                                                       bytes) else seq_data
            )
        except Exception as e:
            result["sequence_error"] = f"Could not fetch sequence for {protein_id}: {e}"

        # Try to get pathway information
        try:
            # Extract gene ID if it's in the format 'organism:gene'
            if ":" in protein_id:
                gene_id = protein_id
            else:
                gene_id = protein_id

            pathway_handle = REST.kegg_link("pathway", gene_id)
            pathway_data = pathway_handle.read()
            pathway_handle.close()
            result["pathways"] = (
                pathway_data.decode("utf-8")
                if isinstance(pathway_data, bytes)
                else pathway_data
            )
        except Exception as e:
            result["pathway_error"] = f"Could not fetch pathways for {protein_id}: {e}"

        return result

    except Exception as e:
        return {"error": f"General error fetching {protein_id} from KEGG: {e}"}


def search_protein_kegg(query: str, organism: str = None) -> list:
    """
    Search for proteins in KEGG database using a query string.

    This function searches the KEGG database for proteins matching a query.
    It can search by gene name, protein name, or other identifiers.

    Args:
        query (str): Search query (e.g., 'p53', 'insulin', 'BRCA1')
        organism (str, optional): Organism code (e.g., 'hsa' for human, 'mmu' for mouse)
                                 If None, searches across all organisms

    Returns:
        list: A list of dictionaries containing search results, each with:
              - 'id': KEGG identifier
              - 'description': protein description
              - 'organism': organism code

    Example:
        >>> results = search_protein_kegg('p53', 'hsa')
        >>> for result in results:
        ...     print(f"{result['id']}: {result['description']}")
    """
    try:
        search_results = []

        # Construct search query
        if organism:
            search_query = f"{organism}:{query}"
        else:
            search_query = query

        # Search in genes database
        try:
            search_handle = REST.kegg_find("genes", search_query)
            search_data = search_handle.read()
            search_handle.close()

            if isinstance(search_data, bytes):
                search_data = search_data.decode("utf-8")

            # Parse search results
            for line in search_data.strip().split("\n"):
                if line.strip():
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        gene_id = parts[0].strip()
                        description = parts[1].strip()

                        # Extract organism from gene_id if present
                        if ":" in gene_id:
                            org_code = gene_id.split(":")[0]
                        else:
                            org_code = organism or "unknown"

                        search_results.append(
                            {
                                "id": gene_id,
                                "description": description,
                                "organism": org_code,
                            }
                        )

        except Exception as e:
            print(f"Error searching KEGG genes: {e}")

        return search_results

    except Exception as e:
        print(f"General error searching KEGG: {e}")
        return []


def fetch_protein_info_kegg_batch(protein_ids: list) -> dict:
    """
    Fetch protein information from KEGG database for multiple protein IDs.

    This function efficiently retrieves protein information for multiple proteins
    from KEGG database with JSON-based caching and rate limiting to respect API constraints.
    Since KEGG API doesn't support true batch operations, this uses optimized sequential
    fetching with comprehensive caching to avoid repeated requests.

    Args:
        protein_ids (list): A list of KEGG protein/gene identifiers

    Returns:
        dict: A dictionary mapping each protein ID to its KEGG information
              or error message if fetching failed

    Raises:
        TypeError: If protein_ids is not a list

    Example:
        >>> ids = ['hsa:7157', 'hsa:672']
        >>> results = fetch_protein_info_kegg_batch(ids)
        >>> for protein_id, info in results.items():
        ...     if 'error' not in info:
        ...         print(f"{protein_id}: Found protein info")
    """
    if not isinstance(protein_ids, list):
        raise TypeError("protein_ids must be a list")
    if not protein_ids:
        return {}

    # Separate NCBI IDs for conversion
    ncbi_ids = [pid for pid in protein_ids if _is_ncbi_accession(pid)]
    other_ids = [pid for pid in protein_ids if not _is_ncbi_accession(pid)]

    # Convert NCBI IDs to UniProt IDs
    conversion_map = {}
    reverse_conversion_map = {}
    if ncbi_ids:
        print(f"--- Converting {len(ncbi_ids)} NCBI IDs to UniProt IDs ---")
        conversion_map = convert_refseq_to_uniprot_batch(ncbi_ids)
        converted_uniprot_ids = list(conversion_map.values())
        print(
            f"--- Successfully converted {len(converted_uniprot_ids)} IDs ---")

        # Create a reverse map to link UniProt IDs back to original NCBI IDs
        reverse_conversion_map = {v: k for k, v in conversion_map.items()}

        # Combine original non-NCBI IDs with newly converted UniProt IDs
        protein_ids_to_fetch = other_ids + converted_uniprot_ids
    else:
        protein_ids_to_fetch = other_ids

    # Check cache first with the IDs that are ready for KEGG (converted UniProt IDs + other IDs)
    cached_results, remaining_ids = _load_cached_results(
        protein_ids_to_fetch, "kegg")
    all_protein_info = {}

    # Process cached results and map them back to original IDs
    for protein_id, cached_info in cached_results.items():
        # Check if this is a converted UniProt ID that should be mapped back to NCBI ID
        original_id = reverse_conversion_map.get(protein_id, protein_id)
        all_protein_info[original_id] = cached_info

    # Handle NCBI IDs that couldn't be converted to UniProt
    for ncbi_id in ncbi_ids:
        if ncbi_id not in conversion_map:
            all_protein_info[ncbi_id] = {
                "skipped": True,
                "reason": f"Could not convert NCBI ID {ncbi_id} to UniProt ID",
                "note": "Conversion failed - unable to find UniProt cross-reference",
            }

    if not remaining_ids:
        print(
            f"All {len(protein_ids)} proteins found in cache or could not be converted"
        )
        return all_protein_info

    print(
        f"Found {len(cached_results)} cached, fetching {len(remaining_ids)} from KEGG"
    )

    # Fetch remaining proteins with caching
    batch_results = {}
    for i, protein_id in enumerate(remaining_ids):
        print(
            f"Fetching KEGG info for {protein_id} ({i + 1}/{len(remaining_ids)})")

        try:
            info = fetch_protein_info_kegg(protein_id)
            batch_results[protein_id] = info

            # Add delay between requests to respect API rate limits
            if i < len(remaining_ids) - 1:
                time.sleep(REQUEST_DELAY_SECONDS)

        except Exception as e:
            error_info = {"error": f"Failed to fetch {protein_id}: {e}"}
            batch_results[protein_id] = error_info

    # Cache the newly fetched results
    if batch_results:
        _cache_batch_results(batch_results, "kegg")

    # Map newly fetched results back to original IDs
    for protein_id, info in batch_results.items():
        original_id = reverse_conversion_map.get(protein_id, protein_id)
        all_protein_info[original_id] = info

    return all_protein_info


def fetch_protein_info_kegg_batch_optimized(protein_ids: list) -> dict:
    """
    Attempt to fetch KEGG protein information using any available batch operations.

    This function explores KEGG API capabilities for batch operations like kegg_link
    and kegg_list to reduce the number of API calls where possible, falling back to
    cached sequential fetching when true batch operations aren't available.

    Args:
        protein_ids (list): A list of KEGG protein/gene identifiers

    Returns:
        dict: A dictionary mapping each protein ID to its KEGG information

    Example:
        >>> ids = ['hsa:7157', 'hsa:672']
        >>> results = fetch_protein_info_kegg_batch_optimized(ids)
    """
    if not isinstance(protein_ids, list):
        raise TypeError("protein_ids must be a list")
    if not protein_ids:
        return {}

    # Check cache first
    cached_results, remaining_ids = _load_cached_results(protein_ids, "kegg")
    all_protein_info = cached_results.copy()

    if not remaining_ids:
        print(f"All {len(protein_ids)} proteins found in cache")
        return all_protein_info

    print(
        f"Found {len(cached_results)} cached, attempting batch fetch for {len(remaining_ids)} from KEGG"
    )

    # Try to use KEGG batch operations where possible
    batch_results = {}

    try:
        # Attempt to get basic info for all IDs at once using kegg_list
        # This works for some KEGG database entries
        ids_string = " ".join(remaining_ids)
        list_handle = REST.kegg_list("genes", ids_string)
        list_data = list_handle.read()
        list_handle.close()

        if isinstance(list_data, bytes):
            list_data = list_data.decode("utf-8")

        # Parse batch results if any
        batch_basic_info = {}
        for line in list_data.strip().split("\n"):
            if line.strip() and "\t" in line:
                parts = line.split("\t", 1)
                if len(parts) >= 2:
                    gene_id = parts[0].strip()
                    description = parts[1].strip()
                    batch_basic_info[gene_id] = {
                        "entry": f"{gene_id}\t{description}"}

        print(f"Got batch basic info for {len(batch_basic_info)} proteins")

        # For proteins that got basic info, fetch additional details individually
        # For others, fetch everything individually
        for protein_id in remaining_ids:
            try:
                if protein_id in batch_basic_info:
                    # Start with batch info and add details
                    info = batch_basic_info[protein_id].copy()
                else:
                    info = {}

                # Add sequence and pathway info (these typically need individual calls)
                try:
                    seq_handle = REST.kegg_get(protein_id, "aaseq")
                    seq_data = seq_handle.read()
                    seq_handle.close()
                    info["sequence"] = (
                        seq_data.decode("utf-8")
                        if isinstance(seq_data, bytes)
                        else seq_data
                    )
                except Exception as e:
                    info["sequence_error"] = (
                        f"Could not fetch sequence for {protein_id}: {e}"
                    )

                try:
                    pathway_handle = REST.kegg_link("pathway", protein_id)
                    pathway_data = pathway_handle.read()
                    pathway_handle.close()
                    info["pathways"] = (
                        pathway_data.decode("utf-8")
                        if isinstance(pathway_data, bytes)
                        else pathway_data
                    )
                except Exception as e:
                    info["pathway_error"] = (
                        f"Could not fetch pathways for {protein_id}: {e}"
                    )

                all_protein_info[protein_id] = info
                batch_results[protein_id] = info

                # Small delay between individual requests
                time.sleep(REQUEST_DELAY_SECONDS)

            except Exception as e:
                error_info = {"error": f"Failed to fetch {protein_id}: {e}"}
                all_protein_info[protein_id] = error_info
                batch_results[protein_id] = error_info

    except Exception as e:
        print(
            f"Batch operation failed, falling back to individual fetching: {e}")
        # Fallback to individual fetching for all remaining IDs
        return fetch_protein_info_kegg_batch(protein_ids)

    # Cache the newly fetched results
    if batch_results:
        _cache_batch_results(batch_results, "kegg")

    return all_protein_info
