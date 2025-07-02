import pandas as pd
from Bio import Entrez, SeqIO
import ssl
import urllib.request
import re
from collections import defaultdict
import plotly.graph_objects as go
import time
from scipy.cluster import hierarchy
import numpy as np


Entrez.email = "dominik@hildania.de"
# Fix SSL issues
ssl_context = ssl.create_default_context()
ssl_context.check_hostname = False
ssl_context.verify_mode = ssl.CERT_NONE
opener = urllib.request.build_opener(
    urllib.request.HTTPSHandler(context=ssl_context))
urllib.request.install_opener(opener)


class Protein:
    seq_id: int
    sci_identifier: str
    desc: str
    aa_start: int
    aa_stop: int
    peptide: str
    hits: dict
    std_dev: float
    info: dict

    def __init__(
        self,
        seq_id: int,
        sci_identifier: str,
        desc: str,
        aa_start: int,
        aa_stop: int,
        peptide: str,
        std_dev: float,
        info: dict,
        **kwargs,
    ):
        self.seq_id = seq_id
        self.sci_identifier = sci_identifier
        self.desc = desc
        self.aa_start = aa_start
        self.aa_stop = aa_stop
        self.peptide = peptide
        self.std_dev = std_dev
        self.info = info
        self.hits = {
            "K004A1": None,
            "K004A2": None,
            "K225A1": None,
            "K225A2": None,
            "K003B1": None,
            "K003B2": None,
            "K237A1": None,
            "K237A2": None,
            "K022B1": None,
            "K022B2": None,
            "K302A1": None,
            "K302A2": None,
        }

        for name, arg in kwargs.items():
            if str(name) not in self.hits.keys():
                raise ValueError(
                    "Arg: "
                    + str(name)
                    + " is not part of the hits prefab: "
                    + str(self.hits.keys())
                )
            self.hits[str(name)] = int(str(arg).replace(".", ""))

        for key in self.hits.keys():
            if self.hits.get(key) == None:
                raise ValueError("Wrong")


def fetch_protein_info(protein_id):
    """
    Fetch protein information using protein accession ID
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


BATCH_SIZE = 50
REQUEST_DELAY_SECONDS = 0.5


def fetch_protein_info_batch(protein_ids: list) -> dict:
    """
    Fetch protein information for a list of protein accession IDs in batches.

    Args:
        protein_ids: A list of protein accession IDs (e.g., ['NP_000047.2']).

    Returns:
        A dictionary mapping each protein ID to its Bio.SeqRecord.SeqRecord
        object, or None if fetching failed for that ID.
    """
    if not isinstance(protein_ids, list):
        raise TypeError("protein_ids must be a list")
    if not protein_ids:
        return {}

    all_protein_records = {}

    for i in range(0, len(protein_ids), BATCH_SIZE):
        batch_ids = protein_ids[i: i + BATCH_SIZE]
        print(
            f"Processing batch {i // BATCH_SIZE + 1}/{(len(protein_ids) + BATCH_SIZE - 1) // BATCH_SIZE}..."
        )

        try:
            fetch_handle = Entrez.efetch(
                db="protein",
                id=",".join(batch_ids),
                # idtype="acc" is crucial for telling Entrez these are accessions.
                idtype="acc",
                rettype="gb",
                retmode="text",
            )
            # Use SeqIO.parse for multiple records; SeqIO.read is for a single record.
            batch_records = SeqIO.parse(fetch_handle, "genbank")
            current_batch_results = {rec.id: rec for rec in batch_records}
            fetch_handle.close()

            for protein_id in batch_ids:
                # The returned record ID (rec.id) includes the version, e.g., "NP_000047.2".
                # Match it against the original requested ID.
                found_record = None
                for rec_id_in_batch, record_obj in current_batch_results.items():
                    if rec_id_in_batch.startswith(protein_id):
                        found_record = record_obj
                        break
                all_protein_records[protein_id] = found_record

        except Exception as e:
            print(f"Error fetching batch starting with {batch_ids[0]}: {e}")
            for protein_id in batch_ids:
                all_protein_records[protein_id] = None

        if i + BATCH_SIZE < len(protein_ids):
            time.sleep(REQUEST_DELAY_SECONDS)

    return all_protein_records


def extract_gene_info(description):
    """
    Extract more detailed gene information for better grouping
    """
    # Remove organism info
    desc_clean = re.sub(r"\[.*?\]", "", description).strip()

    # Remove PREDICTED: prefix
    desc_clean = re.sub(r"^PREDICTED:\s*", "", desc_clean, flags=re.IGNORECASE)

    # Extract potential gene name patterns
    gene_patterns = [
        r"(\w+)\s+isoform",  # Gene name before isoform
        r"(\w+)\s+variant",  # Gene name before variant
        r"^([A-Z0-9]+)\s",  # Capital gene name at start
    ]

    # Remove isoform information
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


def grouping(protein_data: [Protein]):
    """
    Advanced grouping with better gene name extraction
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


def renderPlot(data, x, y):
    num_rows = len(data.index)
    num_cols = len(data.columns)
    cell_size = 30  # pixels per cell
    margin = 200  # space for labels, colorbar, etc.

    width = num_cols * cell_size + margin
    height = num_rows * cell_size + margin

    fig = go.Figure(
        go.Heatmap(
            z=data,
            x=x,
            y=y,
        )
    )

    fig.update_layout(
        width=width,
        height=height,
        xaxis=dict(
            scaleanchor="y",
            scaleratio=1,
        ),
        yaxis=dict(
            scaleanchor="x",
        ),
    )

    fig.show()


if __name__ == "__main__":
    df = pd.read_csv(
        "../resources/WWOP230228_PHSN_C12-2_Top100_Variant_Peptides_blinded.csv",
        dtype={
            "id": str,
            "protein": str,
            "description": str,
            "aa_start": int,
            "aa_stop": int,
            "peptide": str,
            "K004A1": str,
            "K004A2": str,
            "K225A1": str,
            "K225A2": str,
            "K003B1": str,
            "K003B2": str,
            "K237A1": str,
            "K237A2": str,
            "K022B1": str,
            "K022B2": str,
            "K302A1": str,
            "K302A2": str,
            "std_dev": float,
        },
    )
    print(df.to_string())
    print(
        "{:18}".format("Protein")
        + "{:80}".format("Extracted")
        + "{:100}".format("From")
        + "Values"
    )

    # --- REFACTORED LOGIC ---
    # 1. Get a list of unique protein IDs from the DataFrame to avoid redundant fetches.
    unique_protein_ids = df["Protein"].unique().tolist()
    print(f"Found {len(unique_protein_ids)} unique protein IDs to fetch.")

    # 2. Fetch all protein information in one batched operation.
    # This returns a dictionary: {'protein_id': SeqRecord, ...}
    protein_info_map = fetch_protein_info_batch(unique_protein_ids)
    print("\n--- Protein fetching complete. Processing DataFrame rows. ---\n")
    # --- END OF REFACTORED LOGIC ---

    tops: list[Protein] = []
    for row in df.itertuples():
        # Get the pre-fetched info for the protein in the current row.
        # Use .get() to safely handle cases where a protein was not found.
        info = protein_info_map.get(row.Protein)

        # Skip rows where protein info could not be fetched.
        if info is None:
            print(
                f"Skipping row for protein '{row.Protein}' as its info could not be fetched."
            )
            continue

        prot = Protein(
            seq_id=row.SeqID,
            sci_identifier=row.Protein,
            desc=row.Description,
            aa_start=row.AA_Start,
            aa_stop=row.AA_Stop,
            peptide=row.Peptide,
            std_dev=row.std_dev,
            info=info,
            K004A1=row.K004A1,
            K004A2=row.K004A2,
            K225A1=row.K225A1,
            K225A2=row.K225A2,
            K003B1=row.K003B1,
            K003B2=row.K003B2,
            K237A1=row.K237A1,
            K237A2=row.K237A2,
            K022B1=row.K022B1,
            K022B2=row.K022B2,
            K302A1=row.K302A1,
            K302A2=row.K302A2,
        )

        extracted_gene = extract_gene_info(prot.info.description)
        print(
            "{:18}".format(prot.sci_identifier)
            + "{:80}".format(extracted_gene)
            # Truncate for display
            + "{:100}".format(prot.info.description[:100])
            + str(list(prot.hits.values()))
        )
        tops.append(prot)

    for protein in tops:
        print(protein.sci_identifier, str(protein.hits.values()))

    grouped = grouping(tops)

    print("=== ADVANCED GROUPING ===")
    print("Count, Base Name, Proteins")
    for base_name, proteins in grouped.items():
        print(f"{'{:13}'.format(len(proteins))}, {'{:40}'.format(base_name)}")
        for protein in proteins:
            print(f"{protein['description']}")

    print("Render plot")

    row_names = []
    column_names = []
    matrix = []

    for protein in tops:
        protein: Protein = protein
        row = []
        for col_name, value in protein.hits.items():
            row.append(value)
            if col_name not in column_names:
                column_names.append(col_name)
        matrix.append(row)
        row_names.append(protein.sci_identifier)

    matrix_frame = pd.DataFrame(
        data=matrix, index=row_names, columns=column_names)

    print(matrix_frame.to_string())

    renderPlot(matrix_frame, column_names, row_names)

    print("--- Original DataFrame Head ---")
    print(matrix_frame.head())

    # --- 2. Perform Hierarchical Clustering on ROWS ONLY ---

    # The linkage function performs the clustering
    row_linkage = hierarchy.linkage(
        matrix_frame.values, method="complete", metric="euclidean"
    )

    # The dendrogram function calculates the optimal leaf ordering for rows
    row_dendrogram = hierarchy.dendrogram(row_linkage, no_plot=True)

    # The 'leaves' key gives the indices of the rows in their new, clustered order
    clustered_row_indices = row_dendrogram["leaves"]

    # Get the actual row names in this new order
    clustered_row_names = matrix_frame.index[clustered_row_indices]

    # --- 3. Create the New DataFrame with Clustered Rows ---
    # Use .reindex() with only the 'index' argument to reorder rows
    # and leave the columns in their original order.
    clustered_rows_df = matrix_frame.reindex(index=clustered_row_names)

    print("\n--- DataFrame with Clustered Rows ---")
    print("Rows are reordered based on similarity; columns are unchanged.")
    print(clustered_rows_df.head())

    renderPlot(matrix_frame, column_names, clustered_row_names)

    # --- 2. Get Cluster Labels ---
    # Define a distance threshold to cut the dendrogram. You may need to adjust this.
    # A smaller 't' will result in more, smaller clusters.
    distance_threshold = 50000
    cluster_labels = hierarchy.fcluster(
        row_linkage, t=distance_threshold, criterion="distance"
    )
    # `cluster_labels` is a NumPy array like: array([2, 1, 1, 3, 2, ...])

    # --- 3. Group Row Names into a Dictionary ---
    # defaultdict(list) creates a dictionary where each new key automatically
    # gets an empty list as its value.
    clusters_dict = defaultdict(list)

    # zip pairs each row name with its corresponding cluster label
    for row_name, label in zip(matrix_frame.index, cluster_labels):
        clusters_dict[label].append(row_name)

    # --- 4. View the Result ---
    print("--- Clusters stored in a dictionary ---")
    for cluster_id, peptides in clusters_dict.items():
        print(f"\nCluster {cluster_id}:")
        # Using np.array for cleaner printing of the list
        print(str(peptides))
