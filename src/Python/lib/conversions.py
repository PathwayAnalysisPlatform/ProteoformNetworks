import os
from os import path

import networkx as nx
import pandas as pd

from config import LEVELS
from interaction_network import add_nodes, connect_reaction_participants, connect_complex_components
from lib.dictionaries import convert_tab_to_dict
from network_topology_queries import get_reaction_participants, get_complex_components


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
    unique_genes = pd.read_csv(path_file_genes, header=None)
    unique_genes.columns = ['gene']
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


def read_or_create_mapping_proteins_to_level(path_result, level):

    if level == "proteins":
        print("Wrong level to map, can not map proteins to proteins.")
        return
    elif level not in LEVELS:
        print(f"Unkown entity level {level}, expected {LEVELS[0]}, {LEVELS[1]} or {LEVELS[2]}")
        return

    file_name = "mapping_proteins_to_" + level + ".tsv"
    mapping_file = path.join(str(path_result), file_name)
    map = {}
    if not path.exists(mapping_file):

        print("Creating file")
        # Get the full list of proteins
        G = nx.Graph()
        to_level = "proteoforms" if level == "proteoforms" else "proteins"
        participants = get_reaction_participants(to_level, True, True)
        components = get_complex_components(to_level, True, True)
        add_nodes(G, participants, to_level)
        add_nodes(G, components, to_level)

        # For each node, write the protein and the predecesor gene
        with open(mapping_file, 'w') as fh:
            for e, p in G.nodes(data='prevId'):
                if level == "genes":
                    fh.write(f"{e}\t{p}\n")
                    map[e] = p
                else:
                    fh.write(f"{p}\t{e}\n")
                    map[p] = e
    else:
        with open(mapping_file, 'r') as fh:
            for line in fh:
                e, p = line.split("\t")
                map[e] = p
    return map


if __name__ == '__main__':
    print(f"Working directory: {os.getcwd()}")
    read_or_create_mapping_proteins_to_level("../../../resources/UniProt/", "genes")
    read_or_create_mapping_proteins_to_level("../../../resources/UniProt/", "proteoforms")
    # graphs_no_sm = [read_or_create_full_graph(level, False, graphs_path, v=False) for level in LEVELS]

    print("Finished!!")