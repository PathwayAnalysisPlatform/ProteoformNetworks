import os
from os import path

import networkx as nx

from config import LEVELS
from interaction_network import connect_reaction_participants, connect_complex_components, add_nodes
from network_topology_queries import get_reaction_participants, get_complex_components


def read_or_create_full_graph(level, sm=True, graphs_path="", v=True):
    out_file = level
    if not sm:
        out_file += "_no_small_molecules"
    out_file += "_edge_list"
    f = path.join(str(graphs_path), out_file)
    if path.exists(f):
        print("Reading graph")
        G = nx.read_adjlist(f)
        return G

    else:
        print("Creating graph")
        G = nx.Graph()
        G.graph["level"] = level
        G.graph["sm"] = sm

        participants = get_reaction_participants(level, sm, True)
        components = get_complex_components(level, sm, True)

        add_nodes(G, participants, level)
        add_nodes(G, components, level)

        connect_reaction_participants(G, participants, level, False)
        connect_complex_components(G, components, level, False)

        if v:
            print(f"Storing network level: {level} - small molecules: {sm}")

        if len(str(graphs_path)) > 0:
            if not os.path.exists(graphs_path):
                os.mkdir(graphs_path)
        fh = open(os.path.join(str(graphs_path), out_file), 'wb')
        nx.write_adjlist(G, fh)

        print(f"Created graph level: {level} - small molecules: {sm}")
        fh.close()
    print(f"Graph size: {G.size()}")
    return G


if __name__ == '__main__':
    graphs_path = "full_graphs/"
    graphs_path_no_sm = "full_graphs_no_sm/"

    graphs = [read_or_create_full_graph(level, True, graphs_path, v=False) for level in LEVELS]
    graphs_no_sm = [read_or_create_full_graph(level, False, graphs_path, v=False) for level in LEVELS]

    print("Finished!!")
