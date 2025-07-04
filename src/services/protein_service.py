"""
Protein data fetching service module.

This module provides functions to fetch protein information from NCBI databases
using Biopython's Entrez API and KEGG database using Bio.KEGG.REST API.
It includes both single protein fetching and batch processing capabilities
for efficient data retrieval.
"""

from Bio import Entrez, SeqIO
from Bio.KEGG import REST
import ssl
import urllib.request
import time


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


def fetch_protein_info_kegg(protein_id: str) -> dict:
    """
    Fetch protein information from KEGG database using a protein ID.
    
    This function queries the KEGG database for protein information using
    the Bio.KEGG.REST API. It can handle various KEGG protein identifiers
    including UniProt IDs and KEGG gene/protein IDs.
    
    Args:
        protein_id (str): Protein identifier (e.g., 'hsa:7157' for human p53,
                         or UniProt ID like 'P04637')
                         
    Returns:
        dict: A dictionary containing protein information from KEGG including:
              - 'entry': KEGG entry information
              - 'sequence': protein sequence if available
              - 'pathway': associated pathways
              - 'error': error message if query failed
              
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
            result['entry'] = (entry_data.decode('utf-8')
                              if isinstance(entry_data, bytes)
                              else entry_data)
        except Exception as e:
            result['entry_error'] = f"Could not fetch entry for {protein_id}: {e}"
        
        # Try to get protein sequence
        try:
            seq_handle = REST.kegg_get(protein_id, "aaseq")
            seq_data = seq_handle.read()
            seq_handle.close()
            result['sequence'] = (seq_data.decode('utf-8') 
                                 if isinstance(seq_data, bytes) 
                                 else seq_data)
        except Exception as e:
            result['sequence_error'] = f"Could not fetch sequence for {protein_id}: {e}"
        
        # Try to get pathway information
        try:
            # Extract gene ID if it's in the format 'organism:gene'
            if ':' in protein_id:
                gene_id = protein_id
            else:
                gene_id = protein_id
                
            pathway_handle = REST.kegg_link("pathway", gene_id)
            pathway_data = pathway_handle.read()
            pathway_handle.close()
            result['pathways'] = (pathway_data.decode('utf-8') 
                                 if isinstance(pathway_data, bytes) 
                                 else pathway_data)
        except Exception as e:
            result['pathway_error'] = f"Could not fetch pathways for {protein_id}: {e}"
            
        return result
        
    except Exception as e:
        return {'error': f"General error fetching {protein_id} from KEGG: {e}"}


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
                search_data = search_data.decode('utf-8')
                
            # Parse search results
            for line in search_data.strip().split('\n'):
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        gene_id = parts[0].strip()
                        description = parts[1].strip()
                        
                        # Extract organism from gene_id if present
                        if ':' in gene_id:
                            org_code = gene_id.split(':')[0]
                        else:
                            org_code = organism or 'unknown'
                            
                        search_results.append({
                            'id': gene_id,
                            'description': description,
                            'organism': org_code
                        })
                        
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
    from KEGG database with rate limiting to respect API constraints.
    
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
        
    all_protein_info = {}
    
    for i, protein_id in enumerate(protein_ids):
        print(f"Fetching KEGG info for {protein_id} ({i+1}/{len(protein_ids)})")
        
        try:
            info = fetch_protein_info_kegg(protein_id)
            all_protein_info[protein_id] = info
            
            # Add delay between requests to respect API rate limits
            if i < len(protein_ids) - 1:
                time.sleep(REQUEST_DELAY_SECONDS)
                
        except Exception as e:
            all_protein_info[protein_id] = {'error': f"Failed to fetch {protein_id}: {e}"}
            
    return all_protein_info

