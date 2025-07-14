from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
import re


@dataclass
class NCBIFeature:
    """Represents a feature annotation from NCBI GenBank file."""
    feature_type: str
    location: str
    qualifiers: Dict[str, Any] = field(default_factory=dict)


@dataclass
class NCBIEntry:
    """
    Comprehensive model for NCBI GenBank data.
    
    This class represents the complete structure of an NCBI GenBank file,
    including header information, features, and sequence data.
    """
    # Header information
    locus: Optional[str] = None
    locus_length: Optional[int] = None
    locus_type: Optional[str] = None
    locus_division: Optional[str] = None
    locus_date: Optional[str] = None
    definition: Optional[str] = None
    accession: Optional[str] = None
    version: Optional[str] = None
    dblink: List[str] = field(default_factory=list)
    keywords: Optional[str] = None
    source: Optional[str] = None
    organism: Optional[str] = None
    taxonomy: List[str] = field(default_factory=list)
    dbsource: Optional[str] = None
    comment: Optional[str] = None
    
    # Feature annotations
    features: List[NCBIFeature] = field(default_factory=list)
    
    # Sequence data
    origin_sequence: Optional[str] = None
    
    # Parsed feature-specific data for easy access
    gene_name: Optional[str] = None
    gene_synonyms: List[str] = field(default_factory=list)
    product: Optional[str] = None
    calculated_mol_wt: Optional[str] = None
    chromosome: Optional[str] = None
    db_xrefs: List[str] = field(default_factory=list)
    regions: List[Dict[str, str]] = field(default_factory=list)
    sites: List[Dict[str, str]] = field(default_factory=list)

    def __str__(self):
        return str(self.__dict__)

    @classmethod
    def from_genbank_text(cls, genbank_text: str) -> "NCBIEntry":
        """Parse a GenBank text entry into an NCBIEntry object."""
        ncbi_entry = cls()
        lines = genbank_text.strip().split("\n")
        
        current_section = None
        current_feature = None
        current_qualifier = None
        
        i = 0
        while i < len(lines):
            line = lines[i]
            
            # Handle empty lines
            if not line.strip():
                i += 1
                continue
            
            # Check for main sections
            if line.startswith("LOCUS"):
                current_section = "LOCUS"
                ncbi_entry._parse_locus_line(line)
            elif line.startswith("DEFINITION"):
                current_section = "DEFINITION"
                ncbi_entry.definition = line[12:].strip()
            elif line.startswith("ACCESSION"):
                current_section = "ACCESSION"
                ncbi_entry.accession = line[12:].strip()
            elif line.startswith("VERSION"):
                current_section = "VERSION"
                ncbi_entry.version = line[12:].strip()
            elif line.startswith("DBLINK"):
                current_section = "DBLINK"
                ncbi_entry.dblink.append(line[12:].strip())
            elif line.startswith("KEYWORDS"):
                current_section = "KEYWORDS"
                ncbi_entry.keywords = line[12:].strip()
            elif line.startswith("SOURCE"):
                current_section = "SOURCE"
                ncbi_entry.source = line[12:].strip()
            elif line.startswith("  ORGANISM"):
                current_section = "ORGANISM"
                ncbi_entry.organism = line[12:].strip()
            elif line.startswith("DBSOURCE"):
                current_section = "DBSOURCE"
                ncbi_entry.dbsource = line[12:].strip()
            elif line.startswith("COMMENT"):
                current_section = "COMMENT"
                ncbi_entry.comment = line[12:].strip()
            elif line.startswith("FEATURES"):
                current_section = "FEATURES"
            elif line.startswith("ORIGIN"):
                current_section = "ORIGIN"
                ncbi_entry.origin_sequence = ""
            elif line.startswith("//"):
                break
            elif current_section == "FEATURES":
                # Parse features section
                if line.startswith("     ") and not line.startswith("                     "):
                    # New feature
                    feature_match = re.match(r"     (\w+)\s+(.+)", line)
                    if feature_match:
                        if current_feature:
                            ncbi_entry.features.append(current_feature)
                        current_feature = NCBIFeature(
                            feature_type=feature_match.group(1),
                            location=feature_match.group(2).strip()
                        )
                        current_qualifier = None
                elif line.startswith("                     /"):
                    # Feature qualifier
                    qualifier_match = re.match(r'\s+/([^=]+)=?"?([^"]*)"?', line)
                    if qualifier_match and current_feature:
                        qualifier_name = qualifier_match.group(1)
                        qualifier_value = qualifier_match.group(2).strip('"')
                        # Handle multiple db_xref entries
                        if qualifier_name == "db_xref":
                            if qualifier_name not in current_feature.qualifiers:
                                current_feature.qualifiers[qualifier_name] = []
                            current_feature.qualifiers[qualifier_name].append(qualifier_value)
                        else:
                            current_feature.qualifiers[qualifier_name] = qualifier_value
                        current_qualifier = qualifier_name
                elif line.startswith("                     ") and current_qualifier and current_feature:
                    # Continuation of qualifier value
                    continuation = line.strip().strip('"')
                    if current_qualifier in current_feature.qualifiers:
                        current_feature.qualifiers[current_qualifier] += " " + continuation
            elif current_section == "ORIGIN":
                # Parse sequence data
                sequence_match = re.findall(r'[a-zA-Z]+', line)
                if sequence_match:
                    ncbi_entry.origin_sequence += "".join(sequence_match)
            elif current_section in ["DEFINITION", "COMMENT", "DBLINK", "KEYWORDS"] and line.startswith("            "):
                # Continuation lines for multi-line fields
                continuation = line[12:].strip()
                if current_section == "DEFINITION":
                    ncbi_entry.definition += " " + continuation
                elif current_section == "COMMENT":
                    ncbi_entry.comment += " " + continuation
                elif current_section == "DBLINK":
                    ncbi_entry.dblink.append(continuation)
                elif current_section == "KEYWORDS":
                    ncbi_entry.keywords += " " + continuation
            elif current_section == "ORGANISM" and line.startswith("            "):
                # Taxonomy information
                taxonomy_line = line[12:].strip()
                if taxonomy_line.endswith(";"):
                    taxonomy_line = taxonomy_line[:-1]
                ncbi_entry.taxonomy.extend([t.strip() for t in taxonomy_line.split(";")])
            
            i += 1
        
        # Add the last feature if exists
        if current_feature:
            ncbi_entry.features.append(current_feature)
        
        # Extract commonly used information from features
        ncbi_entry._extract_feature_data()
        
        return ncbi_entry

    def _parse_locus_line(self, line: str):
        """Parse the LOCUS line to extract various components."""
        parts = line.split()
        if len(parts) >= 2:
            self.locus = parts[1]
        if len(parts) >= 3:
            try:
                length_part = parts[2]
                if "aa" in length_part or "bp" in length_part:
                    self.locus_length = int(re.findall(r'\d+', length_part)[0])
                    self.locus_type = "aa" if "aa" in length_part else "bp"
            except (ValueError, IndexError):
                pass
        if len(parts) >= 5:
            self.locus_division = parts[4]
        if len(parts) >= 6:
            self.locus_date = parts[5]

    def _extract_feature_data(self):
        """Extract commonly used data from features for easy access."""
        for feature in self.features:
            if feature.feature_type == "source":
                if "organism" in feature.qualifiers:
                    if not self.organism:  # Only set if not already set from header
                        self.organism = feature.qualifiers["organism"]
                if "chromosome" in feature.qualifiers:
                    self.chromosome = feature.qualifiers["chromosome"]
                if "db_xref" in feature.qualifiers:
                    if isinstance(feature.qualifiers["db_xref"], list):
                        self.db_xrefs.extend(feature.qualifiers["db_xref"])
                    else:
                        self.db_xrefs.append(feature.qualifiers["db_xref"])
            
            elif feature.feature_type == "Protein":
                if "product" in feature.qualifiers:
                    self.product = feature.qualifiers["product"]
                if "calculated_mol_wt" in feature.qualifiers:
                    self.calculated_mol_wt = feature.qualifiers["calculated_mol_wt"]
            
            elif feature.feature_type == "CDS":
                if "gene" in feature.qualifiers:
                    self.gene_name = feature.qualifiers["gene"]
                if "gene_synonym" in feature.qualifiers:
                    synonyms = feature.qualifiers["gene_synonym"]
                    if isinstance(synonyms, str):
                        self.gene_synonyms = [s.strip() for s in synonyms.split(";")]
                    elif isinstance(synonyms, list):
                        self.gene_synonyms = synonyms
                if "db_xref" in feature.qualifiers:
                    db_xref_val = feature.qualifiers["db_xref"]
                    if isinstance(db_xref_val, list):
                        self.db_xrefs.extend(db_xref_val)
                    else:
                        self.db_xrefs.append(db_xref_val)
            
            elif feature.feature_type == "Region":
                region_info = {
                    "location": feature.location,
                    "region_name": feature.qualifiers.get("region_name", ""),
                    "note": feature.qualifiers.get("note", ""),
                    "db_xref": feature.qualifiers.get("db_xref", "")
                }
                self.regions.append(region_info)
            
            elif feature.feature_type == "Site":
                site_info = {
                    "location": feature.location,
                    "site_type": feature.qualifiers.get("site_type", ""),
                    "note": feature.qualifiers.get("note", ""),
                    "db_xref": feature.qualifiers.get("db_xref", "")
                }
                self.sites.append(site_info)

    def get_gene_id(self) -> Optional[str]:
        """Extract GeneID from db_xrefs."""
        for xref in self.db_xrefs:
            if xref.startswith("GeneID:"):
                return xref.split(":", 1)[1]
        return None

    def get_hgnc_id(self) -> Optional[str]:
        """Extract HGNC ID from db_xrefs."""
        for xref in self.db_xrefs:
            if xref.startswith("HGNC:"):
                return xref.split(":", 1)[1]
        return None

    def get_mim_id(self) -> Optional[str]:
        """Extract MIM ID from db_xrefs."""
        for xref in self.db_xrefs:
            if xref.startswith("MIM:"):
                return xref.split(":", 1)[1]
        return None

    def get_taxon_id(self) -> Optional[str]:
        """Extract taxon ID from db_xrefs."""
        for xref in self.db_xrefs:
            if xref.startswith("taxon:"):
                return xref.split(":", 1)[1]
        return None

    def has_domain(self, domain_name: str) -> bool:
        """Check if the protein contains a specific domain."""
        for region in self.regions:
            if domain_name.lower() in region.get("region_name", "").lower():
                return True
            if domain_name.lower() in region.get("note", "").lower():
                return True
        return False

    def get_domains(self) -> List[str]:
        """Get all domain names from regions."""
        domains = []
        for region in self.regions:
            if region.get("region_name"):
                domains.append(region["region_name"])
        return domains

