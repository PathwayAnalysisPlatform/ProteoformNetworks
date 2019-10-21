import networkx as nx
import pandas as pd

LEVELS = {"gene", "protein", "proteoform"}


def read_scores(file_name):
    data = pd.read_csv(file_name, sep='\t')
    return data.SCORE


def get_graph(trait, level, path_to_root="../../../"):
    if level not in LEVELS:
        raise ValueError("level must be one of %r." % LEVELS)

    G = nx.Graph()

    # Traverse file with module members. Get set of members for the trait
    print("Reading members of module ", trait)
    path_file_modules = path_to_root + "reports/modules/" + level + '_modules.tsv'
    with open(path_file_modules) as file_modules:
        file_modules.readline()  # Read header
        line = file_modules.readline()
        while line:
            columns = line.split("\t")
            if columns[0] == trait:
                G.add_node(columns[1].strip())
            line = file_modules.readline()

    # Traverse file with interactions of the level. Get the set of edges for this trait
    print("Reading edges in module ", trait)
    path_file_interactions = path_to_root + "resources/Reactome/v70/"
    if level == "gene":
        path_file_interactions += "Genes/geneEdges.tsv"
    elif level == "protein":
        path_file_interactions += "Proteins/proteinEdges.tsv"
    else:
        path_file_interactions += "Proteoforms/proteoformEdges.tsv"

    nodes = set(G.nodes)
    with open(path_file_interactions) as file_interactions:
        file_interactions.readline()  # Read header
        line = file_interactions.readline()
        while line:
            columns = line.split("\t")
            if columns[0] in nodes and columns[1] in nodes:
                G.add_edge(columns[0], columns[1])
            line = file_interactions.readline()

    return G
