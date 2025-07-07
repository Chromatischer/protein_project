from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class KeggEntry:
    entry: Optional[str] = None
    symbol: Optional[str] = None
    name: Optional[str] = None
    position: Optional[str] = None
    motif: Optional[str] = None
    dblinks: List[str] = field(default_factory=list)
    aaseq: Optional[str] = None
    ntseq: Optional[str] = None
    orthology: Optional[str] = None
    organism: Optional[str] = None
    pathway: Optional[str] = None
    brite: List[str] = field(default_factory=list)
    network: Optional[str] = None
    disease: Optional[str] = None
    drug_target: Optional[str] = None

    def __str__(self):
        return str(self.__dict__)

    @classmethod
    def from_kegg_text(cls, text_entry: str) -> "KeggEntry":
        """Parse a KEGG text entry into a KeggEntry object."""
        kegg_entry = cls()
        lines = text_entry.strip().split("\n")
        current_field = None

        for line in lines:
            if not line.strip():
                continue

            if line.startswith(" "):
                # Continuation line
                if current_field == "BRITE":
                    kegg_entry.brite.append(line.strip())
                elif current_field == "DBLINKS":
                    kegg_entry.dblinks.append(line.strip())
                elif current_field == "AASEQ":
                    kegg_entry.aaseq = (kegg_entry.aaseq or "") + line.strip()
                elif current_field == "NTSEQ":
                    kegg_entry.ntseq = (kegg_entry.ntseq or "") + line.strip()
                elif current_field:
                    attr_value = getattr(kegg_entry, current_field.lower())
                    setattr(
                        kegg_entry,
                        current_field.lower(),
                        f"{attr_value} {line.strip()}" if attr_value else line.strip(),
                    )
            else:
                # New field
                parts = line.split(" ", 1)
                field_name = parts[0].strip()
                value = parts[1].strip() if len(parts) > 1 else ""
                current_field = field_name

                if field_name == "ENTRY":
                    kegg_entry.entry = value.split(" ", 1)[0]
                elif field_name in [
                    "SYMBOL",
                    "NAME",
                    "ORTHOLOGY",
                    "ORGANISM",
                    "PATHWAY",
                    "POSITION",
                    "MOTIF",
                ]:
                    setattr(kegg_entry, field_name.lower(), value)
                elif field_name == "DBLINKS":
                    kegg_entry.dblinks.append(value)
                elif field_name == "BRITE":
                    kegg_entry.brite.append(value)
                elif field_name == "AASEQ":
                    kegg_entry.aaseq = value
                elif field_name == "NTSEQ":
                    kegg_entry.ntseq = value

        return kegg_entry
