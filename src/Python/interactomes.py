import codecs
from bisect import bisect_left
from pathlib import Path

import networkx as nx
import pandas as pd

import config
from lib.graph_database import get_query_result
from networks import add_edges_reaction_participants, add_edges_complex_components, add_nodes, save_graph, \
    read_graph
from queries import QUERIES_PARTICIPANTS, fix_neo4j_values, QUERIES_COMPONENTS


def read_or_create_interactome(level, sm=True, graphs_path="", v=False):
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
    mapping_proteins_to_genes = graphs_path + config.MAPPING_FILE.replace('level', 'genes')
    mapping_proteins_to_proteoforms = graphs_path + config.MAPPING_FILE.replace('level', 'proteoforms')

    # Check if files exist
    if Path(json_file).exists():
        G = read_graph(json_file)
    else:
        print(f"Not found {json_file}")
        print("Creating graph")
        G = nx.Graph()
        G.graph['num_' + level] = 0
        G.graph['num_small_molecules'] = 0

        participants = get_query_result(QUERIES_PARTICIPANTS[level])
        participants = fix_neo4j_values(participants, level)

        components = get_query_result(QUERIES_COMPONENTS[level])
        components = fix_neo4j_values(components, level)

        if sm:
            sm_participants = get_query_result(QUERIES_PARTICIPANTS['sm'])
            sm_participants = fix_neo4j_values(sm_participants, 'sm')
            participants = pd.concat([sm_participants, participants])

            sm_components = get_query_result(QUERIES_COMPONENTS['sm'])
            sm_components = fix_neo4j_values(sm_components, 'sm')
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


def index(l, x, lo, hi):
    """Locate the leftmost value exactly equal to x in list l"""
    i = bisect_left(l, x, lo, hi)
    if i != len(l) and l[i] == x:
        return i
    print(f"Could not find value: {x} in range [{lo} -- {hi}]")
    raise ValueError


def index_all_vertices(interactomes):
    """
    Merges all gene, protein, proteoform and small molecule vertices in a single indexed sorted list.
    The complete vertices list is sorted section by section.

    :param interactomes: dict with a networkx graph for each level
    :return: List of all vertices, dict with start indexes, dict with end indexes
    """
    vertices = []

    vertices_genes = [n for n, t in interactomes["genes"].nodes(data='type') if t == "genes"]
    vertices_proteins = [n for n, t in interactomes["proteins"].nodes(data='type') if t == "proteins"]
    vertices_proteoforms = [n for n, t in interactomes["proteoforms"].nodes(data='type') if t == "proteoforms"]
    vertices_small_molecules = [n for n, t in interactomes["genes"].nodes(data='type') if t == "SimpleEntity"]

    vertices_genes.sort()
    vertices_proteins.sort()
    vertices_proteoforms.sort()
    vertices_small_molecules.sort()

    vertices = [*vertices_genes, *vertices_proteins, *vertices_proteoforms, *vertices_small_molecules]

    start_index = {}
    end_index = {}

    start_index["genes"] = 0
    end_index["genes"] = len(vertices_genes) - 1
    start_index["proteins"] = end_index["genes"] + 1
    end_index["proteins"] = start_index["proteins"] + len(vertices_proteins) - 1
    start_index["proteoforms"] = end_index["proteins"] + 1
    end_index["proteoforms"] = start_index["proteoforms"] + len(vertices_proteoforms) - 1
    start_index["SimpleEntity"] = end_index["proteoforms"] + 1
    end_index["SimpleEntity"] = start_index["SimpleEntity"] + len(vertices_small_molecules) - 1

    print("Start indexes:")
    print(start_index)

    print("End indexes:")
    print(end_index)

    return vertices, start_index, end_index


def save_indexed_vertices(indexed_vertices, graphs_path):
    vertices_file = Path(graphs_path + "interactome_vertices.tsv")
    with codecs.open(vertices_file, 'w', "utf-8") as fh:
        for v in indexed_vertices:
            fh.write(f"{v}\n")
    print(f"Created indexed vertices file")
    return


def save_edges_with_indexed_vertices(interactomes, output_path, indexed_vertices, start_indexes, end_indexes):
    """
    Save interactome to a file called "interactome_edges.tsv" which has one edge per row. An edge is a pair of ints.

    :param interactomes: three networkx graphs with string names of molecules as vertices
    :param output_path:
    :param indexed_vertices: list of vertices in the order in which they are indexed
    :param start_indexes: starting position in the indexed list for each type of molecule
    :param end_indexes: ending position in the indexed list for each type of molecule
    :return: void
    """
    edges_file = Path(output_path + "interactome_edges.tsv")
    with codecs.open(edges_file, 'w', "utf-8") as fh:
        for level in config.LEVELS:
            for e in interactomes[level].edges:
                index1 = 0
                index2 = 0
                if e[0].startswith("sm_"):
                    index1 = index(indexed_vertices, e[0], start_indexes["SimpleEntity"], end_indexes["SimpleEntity"])
                else:
                    index1 = index(indexed_vertices, e[0], start_indexes[level], end_indexes[level])

                if e[1].startswith("sm_"):
                    index2 = index(indexed_vertices, e[1], start_indexes["SimpleEntity"], end_indexes["SimpleEntity"])
                else:
                    index2 = index(indexed_vertices, e[1], start_indexes[level], end_indexes[level])
                fh.write(f"{index1}\t{index2}\n")
    print(f"Created edges file with indexed vertices")
    return


def save_ranges(start_indexes, end_indexes, output_path):
    ranges_file = Path(output_path + "interactome_ranges.tsv")
    with codecs.open(ranges_file, 'w', "utf-8") as fh:
        for level in config.LEVELS:
            fh.write(f"{start_indexes[level]}\t{end_indexes[level]}\n")
        fh.write(f"{start_indexes['SimpleEntity']}\t{end_indexes['SimpleEntity']}\n")


def save_interactomes_with_indexed_vertices(interactomes, graphs_path):
    """
    Save interactomes as three separate files with the list of edges using indexed vertices.

    :param interactomes: dict with a networkx graph for each level
    :param graphs_path:
    :return: void
    """
    # Gather the gene, protein, proteoform and small molecule vertices together and index all of them
    indexed_vertices, start_indexes, end_indexes = index_all_vertices(interactomes)
    assert end_indexes["SimpleEntity"] == len(indexed_vertices) - 1
    assert sorted(indexed_vertices[start_indexes["genes"]:end_indexes["genes"]])
    assert sorted(indexed_vertices[start_indexes["proteins"]:end_indexes["proteins"]])
    assert sorted(indexed_vertices[start_indexes["proteoforms"]:end_indexes["proteoforms"]])
    assert sorted(indexed_vertices[start_indexes["SimpleEntity"]:end_indexes["SimpleEntity"]])

    save_indexed_vertices(indexed_vertices, graphs_path)

    # Create three files with the interactomes
    save_edges_with_indexed_vertices(interactomes, graphs_path, indexed_vertices, start_indexes, end_indexes)

    save_ranges(start_indexes, end_indexes, graphs_path)

    return


if __name__ == '__main__':
    # print(f"Working directory: {os.getcwd()}")

    graphs_path = "../../resources/Reactome/"
    interactomes = {l: read_or_create_interactome(l, True, graphs_path) for l in config.LEVELS}
    print(f"Indexing vertices")
    save_interactomes_with_indexed_vertices(interactomes, graphs_path)
