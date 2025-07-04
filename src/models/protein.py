"""
Protein data model module.

This module contains the Protein class which represents a protein with its
sequence information, peptide data, and experimental hit values.
"""


class Protein:
    """
    Represents a protein with sequence information and experimental data.
    
    This class stores protein metadata including sequence identifiers,
    peptide information, and experimental hit values from various samples.
    
    Attributes:
        seq_id (int): Sequence identifier
        sci_identifier (str): Scientific protein identifier (e.g., NP_000047.2)
        desc (str): Description of the protein
        aa_start (int): Starting amino acid position
        aa_stop (int): Ending amino acid position
        peptide (str): Peptide sequence
        hits (dict): Dictionary containing experimental hit values for different samples
        std_dev (float): Standard deviation of measurements
        info (dict): Additional protein information from database
        kegg_info (dict): KEGG database information for the protein
    """
    
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
        kegg_info: dict = None,
        **kwargs,
    ):
        """
        Initialize a Protein instance.
        
        Args:
            seq_id (int): Sequence identifier
            sci_identifier (str): Scientific protein identifier
            desc (str): Description of the protein
            aa_start (int): Starting amino acid position
            aa_stop (int): Ending amino acid position
            peptide (str): Peptide sequence
            std_dev (float): Standard deviation of measurements
            info (dict): Additional protein information from database
            **kwargs: Experimental hit values for samples (K004A1, K004A2, etc.)
            
        Raises:
            ValueError: If an unknown sample name is provided or if any required
                       hit value is missing
        """
        self.seq_id = seq_id
        self.sci_identifier = sci_identifier
        self.desc = desc
        self.aa_start = aa_start
        self.aa_stop = aa_stop
        self.peptide = peptide
        self.std_dev = std_dev
        self.info = info
        self.kegg_info = kegg_info or {}
        
        # Initialize predefined hit samples with None values
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

        # Validate and assign hit values from kwargs
        for name, arg in kwargs.items():
            if str(name) not in self.hits.keys():
                raise ValueError(
                    "Arg: "
                    + str(name)
                    + " is not part of the hits prefab: "
                    + str(self.hits.keys())
                )
            # Convert hit values to integers, removing decimal points
            self.hits[str(name)] = int(str(arg).replace(".", ""))

        # Ensure all hit values are provided
        for key in self.hits.keys():
            if self.hits.get(key) == None:
                raise ValueError("Wrong")

