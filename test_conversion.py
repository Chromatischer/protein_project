from src.services.protein_service import convert_refseq_to_uniprot_batch


def test_conversion():
    refSeqIDs = ["NP_055522.1", "XP_006723469.1"]

    result = convert_refseq_to_uniprot_batch(refSeqIDs)

    print(str(result))


if __name__ == "__main__":
    test_conversion()
