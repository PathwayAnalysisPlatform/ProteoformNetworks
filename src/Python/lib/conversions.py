import os
from pathlib import Path
from sys import path

import networkx as nx
import pandas as pd

import config
from config import LEVELS
from networks import add_nodes, add_edges_reaction_participants, add_edges_complex_components
from lib.dictionaries import convert_tab_to_dict
from queries import get_reaction_participants, get_complex_components


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
    """
    Creates a two column file with the mapping from protein accessions to proteoforms in Reactome.

    :param path_file_genes:
    :param path_result:
    :param file_name_result:
    :param batch_size:
    :return:
    """
    if not Path(path_result + file_name_result).exists():

        print(f"Creating mapping file {path_result + file_name_result}")
        if not Path(path_file_genes).exists():
            raise FileNotFoundError

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


def read_or_create_mapping_proteins_to_level(mapping_file, G, level):
    """
    Create file with a two column table with the mapping from protein accessions to either genes or proteoforms.

    :param path_result: Result file path
    :param level: {"genes", "proteoforms"}
    :return: Dictionary with the mapping from protein accessions to either genes or proteoforms.
    """

    if level == "proteins":
        print("Wrong level to map, can not map proteins to proteins.")
        return
    elif level not in LEVELS:
        print(f"Unkown entity level {level}, expected {LEVELS[0]} or {LEVELS[2]}")
        return

    map = {}
    if not os.path.exists(mapping_file):
        print(f"Creating file {mapping_file}")
        with open(mapping_file, 'w', encoding='utf-8') as fh:
            # For each node, write the protein and the predecessor gene
            for n, p in G.nodes.data('prevId'):
                if level == "genes":
                    fh.write(f"{n}\t{p}\n")
                    map[n] = p
                else:
                    fh.write(f"{p}\t{n}\n")
                    map[p] = n

    else:
        print(f"Reading file {mapping_file}")
        with open(mapping_file, 'r') as fh:
            for line in fh:
                n, p = line.split("\t")
                map[n] = p
    return map


if __name__ == '__main__':
    print(f"Working directory: {os.getcwd()}")

    # create_gene_to_protein_mapping("../../../" + config.GRAPHS_PATH + "genes_vertices.tsv", "../../../resources/UniProt/", "mapping_proteins_to_genes.tsv", 100)

    # read_or_create_mapping_proteins_to_level("../../../resources/UniProt/", "genes")
    # read_or_create_mapping_proteins_to_level("../../../resources/UniProt/", "proteoforms")
    # graphs_no_sm = [read_or_create_full_graph(level, False, graphs_path, v=False) for level in LEVELS]