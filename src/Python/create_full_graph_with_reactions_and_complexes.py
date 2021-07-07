import networkx as nx
import os

from lib.networks import add_edges_reaction_participants, add_edges_complex_components
from queries import get_reactions, get_complexes, get_reaction_participants_by_reaction, \
    get_complex_components_by_complex


def create_full_graph_with_reactions_and_complexes(level, showSmallMolecules, graphs_path="", verbose=False):
    G = nx.Graph()
    G.graph["level"] = level
    G.graph["sm"] = showSmallMolecules

    reactions = get_reactions()
    complexes = get_complexes()

    # reactions = reactions[:5]
    # complexes = complexes[:5]

    index = 0
    for reaction in reactions['stId']:
        df_participants = get_reaction_participants_by_reaction(reaction, level, showSmallMolecules, False)
        add_edges_reaction_participants(G, df_participants, level, False)
        index += 1
        if verbose:
            print(
                f"Added reaction {reaction}\t{index}/{len(reactions)}\t\tNodes: {len(G.nodes)}, \t\t Edges: {len(G.edges)}",
                flush=True)

    index = 0
    for complex in complexes['stId']:
        df_components = get_complex_components_by_complex(complex, level, showSmallMolecules, False)
        add_edges_complex_components(G, df_components, level, False)
        index += 1
        if verbose:
            print(
                f"Added complex {complex}\t{index}/{len(complexes)}\t\tNodes: {len(G.nodes)}, \t\t Edges: {len(G.edges)}",
                flush=True)

    if verbose:
        print(f"Storing network")

    out_file = level
    if not showSmallMolecules:
        out_file += "_no_small_molecules"
    out_file += "_edge_list"

    if len(str(graphs_path)) > 0:
        if not os.path.exists(graphs_path):
            os.mkdir(graphs_path)
    fh = open(os.path.join(str(graphs_path), out_file), 'wb')
    nx.write_edgelist(G, fh, data=False)

    print(f"Created graph {level} - {showSmallMolecules}", flush=True)
    fh.close()
    return G

print("Started creating the graph:")
create_full_graph_with_reactions_and_complexes("genes", True, "", True)
print("Finished!!")