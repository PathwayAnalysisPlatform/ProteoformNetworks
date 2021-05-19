import networkx as nx

from config import LEVELS
from lib.graph_database import get_pathway_name


# def create_pathway_graphs(pathway, graphs_path="../../reports/pathways/", v=False):
#     """
#     Creates interactomes for the pathway in all three levels, with and wihout small molecules
#
#     If pathway does not exists, then returns empty interactomes.
#     :param pathway: Pathway stId string
#     :param graphs_path:
#     :return: Get two lists of 3 pathways.
#     """
#
#     name = get_pathway_name(pathway)
#     if len(name) == 0:
#         return [nx.Graph(), nx.Graph(), nx.Graph()], [nx.Graph(), nx.Graph(), nx.Graph()]
#     else:
#         graphs = [create_graph(pathway, level, True, graphs_path, v, True) for level in LEVELS]
#         graphs_no_small_molecules = [create_graph(pathway, level, False, graphs_path, v, True) for level in LEVELS]
#         return graphs, graphs_no_small_molecules


# def create_graphs(pathways, graphs_path="../../reports/pathways/"):
#     """
#     Create the interaction networks of each pathway of the arguments
#     :param pathways: pandas dataframe with column "stId"
#     :param graphs_path: location for output files
#     :return: List of namedtuples()
#     """
#     result = list()
#     if type(pathways) == pd.DataFrame and len(pathways) > 0:
#         for pathway in pathways['stId']:
#             name = get_pathway_name(pathway)
#             name = name.iloc[0]['Name']
#             print(f"Creating networks for pathway {pathway}")
#             graphs, graphs_no_small_molecules = create_pathway_graphs(pathway, graphs_path)
#
#             gs = Pathway_graphs(graphs[0], graphs_no_small_molecules[0], graphs[1], graphs_no_small_molecules[1],
#                                 graphs[2],
#                                 graphs_no_small_molecules[2])
#             result.append(gs)
#     return result


def merge_graphs(graphs):
    # How does the resulting graph look like?
    # - Vertices: Composition of nodes in all interactomes
    #    - Merge: sets of Reactions, Pathways, Complexes
    #    - Copy value of: Id, Type, Entity Color
    # - Edges: Composition of edges in all interactomes
    full_graph = nx.compose_all(graphs)  # Add all nodes setting  Id, Type, and entity_color

    for graph in graphs:
        for node in graph.nodes:
            full_graph.nodes[node]['reactions'].update(graph.nodes[node]['reactions'])
            full_graph.nodes[node]['pathways'].update(graph.nodes[node]['pathways'])
            full_graph.nodes[node]['roles'].update(graph.nodes[node]['roles'])
            full_graph.nodes[node]['complexes'].update(graph.nodes[node]['complexes'])


if __name__ == '__main__':
    print(f"hello from interaction_networks.py")
