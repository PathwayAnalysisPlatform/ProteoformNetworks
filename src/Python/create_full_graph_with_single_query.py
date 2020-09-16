from pathlib import Path

import networkx as nx
import pandas as pd

from interaction_network import add_edges_reaction_participants, add_edges_complex_components, add_nodes, save_graph, \
    read_graph
from lib.graph_database import get_query_result
from network_topology_queries import get_reaction_participants, get_complex_components, QUERIES_PARTICIPANTS, \
    fix_neo4j_values, QUERIES_COMPONENTS


def read_or_create_full_graph(level, sm=True, graphs_path="", v=False):
    """
    When any of the files does not exist it creates the full interaction network (graph) and stores it in files.
    If all the files exist it reads the graph from the json file.
    - A json file with all vertices and edges attributes.
    - A .tsv file with the accessioned entity nodes (genes, proteins or proteoforms)
    - A .tsv file with the simple entity nodes (small molecules such as O2, H2O or ATP).
    - A .tsv file for the edges.
    - When the level is "proteins" it also creates a file with the mapping from proteins to genes: proteins_to_genes.tsv
    - When the level is "proteoforms" it creates a file with mapping from proteins to protoeforms: proteins_to_proteoforms.tsv

    :param level: One of the LEVEL enum
    :param sm: bool to add small molecules to the interaction network. The file is created independently if sm = true or false
    :param graphs_path: directory to store the files
    :param v: show extra messages to the console
    :return: networkx complete graph using all genes, proteins or proteoforms
    """

    json_file = Path(graphs_path + level + ".json")
    vertices_file = Path(graphs_path + level + "_vertices.tsv")
    small_molecules_file = Path(graphs_path + level + "_small_molecules.tsv")
    edges_file = Path(graphs_path + level + "_interactions.tsv")

    # Check if files exist
    if Path(json_file).exists() \
            and Path(edges_file).exists() \
            and Path(vertices_file).exists() \
            and Path(small_molecules_file).exists():
        G = read_graph(json_file)
    else:
        print("Creating graph")
        G = nx.Graph()
        G.graph['num_' + level] = 0
        G.graph['num_small_molecules'] = 0

        participants = get_reaction_participants(level, True)
        components = get_complex_components(level, True)

        if sm:
            sm_participants = get_query_result(QUERIES_PARTICIPANTS['sm'])
            sm_participants = fix_neo4j_values(sm_participants, level)
            participants = pd.concat([sm_participants, participants])

            sm_components = get_query_result(QUERIES_COMPONENTS['sm'])
            sm_components = fix_neo4j_values(sm_components, level)
            components = pd.concat([sm_components, components])

        records = pd.concat([participants, components])
        add_nodes(G, records, level)

        add_edges_reaction_participants(G, participants)
        add_edges_complex_components(G, components)

        save_graph(G, level, graphs_path)

    G.graph["level"] = level
    G.graph["sm"] = sm

    if v:
        print(f"Graph edges: {G.size()}")
        print(f"Graph nodes: {len(G.nodes)}")
        print(f"Graph {level} nodes: {G.graph['num_' + level]}")
        print(f"Graph small molecule nodes: {G.graph['num_small_molecules']}")

    return G


if __name__ == '__main__':
    graphs_path = "../../resources/Reactome/"
    G = read_or_create_full_graph("genes", True, graphs_path)
    print(f"The final result is: {len(list(G.nodes))}")
