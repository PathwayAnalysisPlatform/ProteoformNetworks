import pandas as pd

from lib.dictionaries import convert_tab_to_dict



def map_ids(ids, from_database_name='GENENAME', to_database_name='ACC'):
    """Map id list using UniProt identifier mapping service (https://www.uniprot.org/help/api_idmapping)\n
    Returns dictionary with mapping."""
    import urllib.parse
    import urllib.request

    url = 'https://www.uniprot.org/uploadlists/'

    params = {
        'from': from_database_name,
        'to': to_database_name,
        'format': 'tab',
        'query': ' '.join(ids),
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
        #print(str(response.decode('utf-8')))
        return convert_tab_to_dict(str(response.decode('utf-8')))


def create_gene_to_protein_mapping(path_file_genes,
                                   path_result,
                                   file_name_result,
                                   batch_size):
    unique_genes = pd.read_csv(path_file_genes)
    # Convert gene to protein ids by batches
    with open(path_result + file_name_result, "w") as file_gene_to_protein:
        start, end = 0, batch_size
        total = len(unique_genes.gene)
        while True:
            end = min(end, total)
            print(f"  Converting genes {start} to {end}")
            selected_ids = list(unique_genes.gene)[start:end]
            print("Converting: \n", selected_ids)
            mapping = map_ids(selected_ids, from_database_name="GENENAME", to_database_name="ACC")
            for key, values in mapping.items():
                for value in values:
                    file_gene_to_protein.write(f"{key}\t{value}\n")

            if end == total:
                break
            start += batch_size
            end += batch_size
    print("Gene --> Protein mapping READY")