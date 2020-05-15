import networkx as nx
import os

from interaction_network import connect_reaction_participants, connect_complex_components
from network_topology_queries import get_reaction_participants, get_complex_components


def create_full_graph(level, sm=True, data=False, graphs_path="", v=True):
    G = nx.Graph()
    G.graph["level"] = level
    G.graph["sm"] = sm

    participants = get_reaction_participants(level, sm, True)
    components = get_complex_components(level, sm, True)

    connect_reaction_participants(G, participants, level, False)
    connect_complex_components(G, components, level, False)

    if v:
        print(f"Storing network level: {level} - small molecules: {sm} - data: {data}")

    out_file = level
    if not sm:
        out_file += "_no_small_molecules"
    out_file += "_edge_list"

    if len(str(graphs_path)) > 0:
        if not os.path.exists(graphs_path):
            os.mkdir(graphs_path)
    fh = open(os.path.join(str(graphs_path), out_file), 'wb')
    nx.write_edgelist(G, fh, data=data)

    print(f"Created graph level: {level} - small molecules: {sm} - data: {data}")
    fh.close()
    return G


print("Started creating the graph for proteins")
graphs_path="full_graphs_with_data/"
create_full_graph("genes", True, True, graphs_path)
create_full_graph("genes", False, True, graphs_path)
create_full_graph("proteins", True, True, graphs_path)
create_full_graph("proteins", False, True, graphs_path)
create_full_graph("proteoforms", True, True, graphs_path)
create_full_graph("proteoforms", False, True, graphs_path)
print("Finished!!")