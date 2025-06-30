import pandas as pd


class Protein:
    seq_id: int
    sci_identifier: str
    desc: str
    aa_start: int
    aa_stop: int
    peptide: str
    hits: dict(
        K004A1=str,
        K004A2=str,
        K225A1=str,
        K225A2=str,
        K003B1=str,
        K003B2=str,
        K237A1=str,
        K237A2=str,
        K002B1=str,
        K002B2=str,
        K302A1=str,
        K302A2=str,
    )
    std_dev: float

    def __init__(
        self, seq_id, sci_identifier, desc, aa_start, aa_stop, peptide, hits, std_dev
    ):
        self.seq_id = seq_id
        self.sci_identifier = sci_identifier
        self.desc = desc
        self.aa_start = aa_start
        self.aa_stop = aa_stop
        self.peptide = peptide
        self.hits = hits
        self.std_dev = std_dev


if __name__ == "__main__":
    print("Hello World")
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
    tops: list[Protein] = []
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
        print(index, protein)
