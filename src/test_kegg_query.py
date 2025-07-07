from services.protein_service import fetch_protein_info_kegg, fetch_kegg_info_batch


def test_query():
    id = "hsa:153579"

    print(str(fetch_protein_info_kegg(id)))


def test_batch():
    proteins = ["hsa:57101", "hsa:64062"]
    print(fetch_kegg_info_batch(proteins))


if __name__ == "__main__":
    test_batch()
