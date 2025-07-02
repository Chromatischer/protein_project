import pandas as pd
from Bio import Entrez, SeqIO
import ssl
import urllib.request
import re
from collections import defaultdict
import plotly.graph_objects as go


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
    tops: [Protein] = []
    for row in df.itertuples():
        (
            index,
            seq_id,
            protein,
            description,
            aa_start,
            aa_stop,
            peptide,
            k004a1,
            k004a2,
            k225a1,
            k225a2,
            k003b1,
            k003b2,
            k237a1,
            k237a2,
            k022b1,
            k022b2,
            k302a1,
            k302a2,
            std_dev,
        ) = row

        info = fetch_protein_info(protein)
        prot = Protein(
            seq_id=seq_id,
            sci_identifier=protein,
            desc=description,
            aa_start=aa_start,
            aa_stop=aa_stop,
            peptide=peptide,
            std_dev=std_dev,
            info=info,
            K004A1=k004a1,
            K004A2=k004a2,
            K225A1=k225a1,
            K225A2=k225a2,
            K003B1=k003b1,
            K003B2=k003b2,
            K237A1=k237a1,
            K237A2=k237a2,
            K022B1=k022b1,
            K022B2=k022b2,
            K302A1=k302a1,
            K302A2=k302a2,
        )
        print(
            "{:18}".format(prot.sci_identifier)
            + "{:80}".format(extract_gene_info(prot.info.description))
            + "{:100}".format(prot.info.description)
            + str(prot.hits.values())
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
