import re

import networkx as nx
import pandas as pd

LEVELS = {"genes", "proteins", "proteoforms"}


def read_scores(file_name):
    data = pd.read_csv(file_name, sep='\t')
    return data.SCORE


def get_trait_file_name(trait):
    return re.sub("[ ,]", "_", trait.replace("\"", ""))


def get_graph(trait, level, path_to_modules):
    """ Create a networkx graph instance for a trait module at the level specified """
    if level not in LEVELS:
        raise ValueError("level must be one of %r." % LEVELS)

    G = nx.Graph()

    # Traverse file with module members. Get set of members for the trait
    print("Reading members of module ", trait)
    prefix = path_to_modules + get_trait_file_name(trait) + "_" + level
    path_file_module_vertices = prefix + "_vertices.tsv"

    with open(path_file_module_vertices) as file_vertices:
        file_vertices.readline()  # Read header
        line = file_vertices.readline()
        while line:
            G.add_node(line.strip())
            line = file_vertices.readline()

    # Traverse file with interactions of the level. Get the set of edges for this trait
    print("Reading interactions in module ", trait)
    path_file_module_edges = prefix + "_edges.tsv"

    with open(path_file_module_edges) as file_edges:
        file_edges.readline()  # Read header
        line = file_edges.readline()
        while line:
            columns = line.split("\t")
            G.add_edge(columns[0], columns[1].strip())
            line = file_edges.readline()

    return G
