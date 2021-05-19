import codecs
import json
import os
from bisect import bisect_left
from pathlib import Path

import networkx as nx
import pandas as pd

import config
from config import no_sm, with_sm, with_unique_sm, sm, LEVELS
from lib.graph_database import get_reaction_participants_by_pathway, get_complex_components_by_pathway, get_pathway_name


def print_interactome_details(g):
    print(f"Graph for {g.graph['level']} ")
    print(f"Graph edges: {g.size()}")
    print(f"Graph nodes: {len(g.nodes)}")
    print(f"Graph {g.graph['level']} nodes: {g.graph['num_entities']}")
    print(f"Graph small molecule nodes: {g.graph['num_small_molecules']}")
    print("\n***********************\n\n")


def get_json_filename(level, method, out_path):
    return Path(out_path + level + "_" + method + ".json")


def add_nodes(G, df, method=with_sm):
    """
    Add the nodes represented by a dataframe records to the graph G

    :param G: Networkx graph to add the nodes
    :param df: Dataframe with the records resulting from the cypher query for participants of reactions or components
    of complexes
    :param level: genes, proteins or proteoforms. To handle the renaming of some fields accordingly
    :return: The same graph G with the added nodes
    """
    level = G.graph["level"]
    sm_id_column = 'Id' if not method == with_unique_sm else 'unique_sm'

    for index, row in df.iterrows():

        # Each node has an id, label, type, entity_color, role, reactions, pathways, complexes.
        # id: unique identifier across all levels and entity types. For example, to avoid duplicated nodes when a gene
        # or a small molecule are called the same.
        # label: gene name, small molecule name, protein name, and for proteoforms the protein name with possible manual
        # annotation of modifications.
        # print(f"Adding node {index} of entity {entity}")
        unique_id = 'sm_' + row[sm_id_column] if row['Type'] == 'SimpleEntity' else row['Id']

        if unique_id not in G.nodes:
            G.add_node(unique_id,
                       label=row['Id'],
                       type=(row['Type'] if row['Type'] == 'SimpleEntity' else level),
                       entity_color=config.get_entity_color(row['Type'], level),
                       roles=set(),
                       reactions=set(),
                       pathways=set(),
                       complexes=set()
                       )
            if row['Type'] == 'SimpleEntity':
                G.graph['num_small_molecules'] += 1
            else:
                G.graph['num_entities'] += 1

        if 'Complex' in df.columns:
            G.nodes[unique_id]['complexes'].add(row['Complex'])
        elif 'Role' in df.columns:
            G.nodes[unique_id]['roles'].add(row['Role'])
            G.nodes[unique_id]['reactions'].add(row['Reaction'])
            G.nodes[unique_id]['pathways'].add(row['Pathway'])

        if (G.nodes[unique_id]['type'] == "SimpleEntity") or G.nodes[row['Id']]['type'] == "genes":
            G.nodes[unique_id]['prevId'] = row['Id']
        elif G.nodes[unique_id]['type'] == "proteins" or G.nodes[row['Id']]['type'] == "proteoforms":
            G.nodes[unique_id]['prevId'] = row['PrevId']

    print(f"{level} level - small molecules: {G.graph['num_small_molecules']}")
    print(f"{level}: {G.graph['num_entities']}")


def read_graph(json_file):
    """
    Read interaction network in the json file into a networkx graph instance

    :param json_file:
    :return: networkx Graph
    """
    G = nx.Graph()
    with open(json_file, "r") as f:
        data = json.load(f)
        G = nx.json_graph.node_link_graph(data)
    return G


def save_json_graph(G, filename):
    print(f"Storing network level: {G.graph['level']} as: {filename}")
    with open(filename, 'w+') as outfile:
        data = nx.json_graph.node_link_data(G)
        json.dump(data, outfile)
    print(f"Created json file.")


def create_graph(pathway, level, sm, graphs_path="", v=False, save=False):
    """
    Converts a set of reactions with its participants into a graph

    :param df_reactions: Pandas dataframe with reactions and its participants
    :param df_complexes: Pandas dataframe with complexes and its components
    :param sm: False shows only EntityWithAccessionSequence participants
    :return: networkx graph with the interaction network representation of the reactions
    """

    G = nx.Graph()
    name = get_pathway_name(pathway)
    if len(name) == 0:
        return G

    G.graph["stId"] = pathway
    G.graph["level"] = level
    G.graph["sm"] = sm
    G.graph['num_' + level] = 0
    G.graph['num_small_molecules'] = 0

    df_reactions = get_reaction_participants_by_pathway(pathway, level, sm, v)
    df_complexes = get_complex_components_by_pathway(pathway, level, sm, v)

    add_nodes(G, df_reactions, level)
    add_nodes(G, df_complexes, level)

    add_edges_reaction_participants(G, df_reactions, level, v)
    add_edges_complex_components(G, df_complexes, level, v)

    if save:
        save_interactome(G, level, graphs_path, label=pathway + ("_no_sm" if not sm else ""))

    print(f"Created graph {pathway} - {level} - {sm}", flush=True)
    return G


def add_edges_from_product(G, c1, c2, v=False):
    for i in c1:
        for j in c2:
            if i != j:
                G.add_edge(i, j)
                if v:
                    print(f"Added edge from: {i} to {j}")


def add_edges(G, inputs, outputs, catalysts, regulators, reaction, v=False):
    if v:
        print(
            f"\n\nReaction: {reaction}:\nInputs: {inputs}\nCatalysts: {catalysts}\nOutputs: {outputs}\nRegulators: {regulators}")
    add_edges_from_product(G, inputs, outputs, v)
    add_edges_from_product(G, catalysts, outputs, v)
    add_edges_from_product(G, regulators, outputs, v)


def add_edges_reaction_participants(G, df, method=with_sm, v=False):
    """
    Add the edges inferred from the reaction participants in the records of the dataframe.

    :param G: Networkx graph with the participants or reactions. The nodes should be already in the graph.
    :param df: The dataframe with all the participants, their role, the reaction and pathway where they participate.
    :param v: v; show extra information messages
    :return: The same graph G with the edges added
    """
    if len(df) == 0:
        return G
    pathway = df.iloc[0]['Pathway']
    reaction = df.iloc[0]['Reaction']
    inputs = set()
    outputs = set()
    catalysts = set()
    regulators = set()

    sm_id_column = 'Id' if not method == with_unique_sm else 'unique_sm'

    for index, participant in df.iterrows():
        unique_id = 'sm_' + participant[sm_id_column] if participant['Type'] == 'SimpleEntity' else participant['Id']
        G.nodes[unique_id]['roles'].add(participant['Role'])
        G.nodes[unique_id]['reactions'].add(participant['Reaction'])
        G.nodes[unique_id]['pathways'].add(participant['Pathway'])
        if reaction != participant['Reaction'] or pathway != participant['Pathway']:
            add_edges(G, inputs, outputs, catalysts, regulators, reaction, v)
            reaction = participant['Reaction']
            pathway = participant['Pathway']
            inputs = set()
            outputs = set()
            catalysts = set()
            regulators = set()
        if participant['Role'] == 'input':
            inputs.add(unique_id)
        elif participant['Role'] == 'output':
            outputs.add(unique_id)
        elif participant['Role'] == 'regulatedBy':
            regulators.add(unique_id)
        elif participant['Role'] == 'catalystActivity':
            catalysts.add(unique_id)
    reaction = df.iloc[-1]['Reaction']
    pathway = df.iloc[-1]['Pathway']
    add_edges(G, inputs, outputs, catalysts, regulators, reaction, v)

    # Convert the reaction set of each node into a list
    for n in G.nodes:
        G.nodes[n]['reactions'] = list(G.nodes[n]['reactions'])
        G.nodes[n]['pathways'] = list(G.nodes[n]['pathways'])
        G.nodes[n]['roles'] = list(G.nodes[n]['roles'])

    if (v):
        print(f"From reactions, added {len(G.nodes)} nodes to the graph.")
        print(f"From reactions, added {len(G.edges)} edges to the graph.")


def add_edges_complex_components(G, df, method=with_sm, v=False):
    """
    Add the edges inferred from the complex components in the records of the dataframe.

    :param G: Networkx graph with the components of the complexes. The nodes should be already in the graph.
    :param df: The dataframe with all the components and the column of the complex to which they belong.
    :param v: v; show extra information messages
    :return: The same graph G with the edges added
    """
    if len(df) == 0:
        return G

    components = set()

    sm_id_column = 'Id' if not method == with_unique_sm else 'unique_sm'

    complex = df.iloc[0]['Complex']
    for index, record in df.iterrows():
        unique_id = 'sm_' + record[sm_id_column] if record['Type'] == 'SimpleEntity' else record['Id']
        G.nodes[unique_id]['complexes'].add(record['Complex'])

        if complex != record['Complex']:
            add_edges_from_product(G, components, components, v=v)
            complex = record['Complex']
            components = set()
        components.add(unique_id)
    add_edges_from_product(G, components, components, v=v)

    # Convert the complex set of each node into a list
    for n in G.nodes:
        G.nodes[n]['complexes'] = list(G.nodes[n]['complexes'])

    if (v):
        print(f"From complexes, added {len(G.nodes)} nodes to the graph.")
        print(f"From complexes, added {len(G.edges)} edges to the graph.")


def save_interactome(G, graphs_path="", label=""):
    """
    Create the json file with all attributes.
    Creates a vertices file for the accessioned entities.
    Creates a vertices file for the simple entities.
    Creates an edges file for the interactions.
    If the level is proteins write the mapping from proteins to genes
    If the level is proteoforms write the mapping from proteins to proteoforms

    :param G: networkx file with the
    :param level: genes, proteins or protoeforms
    :param graphs_path: directory where to store the files
    :param label: string phrase for the file names
    :return: void
    """
    level = G.graph["level"]
    method = G.graph["method"]
    print(f"Storing network level: {level} at directory: {graphs_path}")

    if len(label) > 0:
        label = label + "_"

    json_file = get_json_filename(level, method, graphs_path)
    accessioned_entities_file = Path(graphs_path + label + level + "_accessioned_entities.tsv")
    interactions_file = Path(graphs_path + label + level + "_interactions.tsv")
    small_molecules_file = Path(graphs_path + label + level + "_small_molecules.tsv")

    if len(str(graphs_path)) > 0:
        if not os.path.exists(graphs_path):
            Path(graphs_path).mkdir(parents=True, exist_ok=True)

    with open(json_file, 'w+') as outfile:
        data = nx.json_graph.node_link_data(G)
        json.dump(data, outfile)
    print(f"Created json file.")

    with open(interactions_file, 'wb') as fh:
        nx.write_edgelist(G, fh)
    print(f"Created edges file for {level}")

    with codecs.open(accessioned_entities_file, 'w', "utf-8") as fh:
        for n, t in G.nodes(data='type'):
            if t != "SimpleEntity":
                fh.write('%s\n' % n)
    print(f"Created entity vertices file for {level}")

    with codecs.open(small_molecules_file, 'w', "utf-8") as fh:
        for n, t in G.nodes(data='type'):
            if t == "SimpleEntity":
                fh.write('%s\n' % n)
                G.graph['num_small_molecules'] += 1
    print(f"Created small molecule vertices file for {level}")

    if level == "proteins":
        mapping_file = graphs_path + config.MAPPING_FILE.replace('level', 'genes')
        with codecs.open(mapping_file, 'w', "utf-8") as map_file:
            for n, t in G.nodes(data='type'):
                if t != "SimpleEntity":
                    map_file.write(f"{n}\t{G.nodes[n]['prevId']}\n")
        print(f"Created mapping proteins to genes file.")

    if level == config.proteoforms:
        mapping_file = graphs_path + config.MAPPING_FILE.replace('level', 'proteoforms')
        with codecs.open(mapping_file, 'w', "utf-8") as map_file:
            for n, t in G.nodes(data='type'):
                if t != "SimpleEntity":
                    map_file.write(f"{G.nodes[n]['prevId']}\t{n}\n")
        print(f"Created mapping proteins to proteoforms file.")


def create_interactome(level, method, participants, components, out_path=""):
    """
    Create interaction network with the participants and components provided as parameter

    Creates a networkx graph and stores it in a json file.

    :param level: genes, proteins or proteoforms
    :param method: "no_sm", "with_sm" or "with_unique_sm"
    :param participants: pandas dataframe with the reaction participants
    :param components: pandas dataframe with the complex components
    :param out_path: path to directory to store the json file
    :return: void
    """

    print("Creating interactome file...")
    json_file = get_json_filename(level, method, out_path)

    G = nx.Graph()
    G.graph["level"] = level
    G.graph["method"] = method
    G.graph['num_entities'] = 0
    G.graph['num_small_molecules'] = 0

    if method == no_sm:
        add_nodes(G, pd.concat([participants[level], components[level]]))
        add_edges_reaction_participants(G, participants[level])
        add_edges_complex_components(G, components[level])
    elif method == with_sm or method == with_unique_sm:
        add_nodes(G, pd.concat([participants[level], participants[sm], components[level], components[sm]]), method)
        add_edges_reaction_participants(G, participants[level])
        add_edges_complex_components(G, components[level])
    else:
        raise Exception("No such method to create the interactome")

    save_interactome(G, out_path, method + "_")
    print(f"Finished creating interactome file for {level}-{method}")


def get_interactome(level, method, participants, components, out_path=""):
    """
    Returns a networkx graph instance of the selected interaction network.
    Tries to read from a json file with the network. If the file does not exists it creates it.
    """

    filename = get_json_filename(level, method, out_path)

    if not Path(filename).exists():
        create_interactome(level, method, participants, components, out_path)
    g = read_graph(filename)

    return g


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


def get_nums(interactome_dict):
    num_interactions = pd.Series([interactome_dict[l].size() for l in LEVELS], index=LEVELS)
    num_entities = pd.Series([interactome_dict[l].graph['num_entities'] for l in LEVELS], index=LEVELS)
    num_small_molecules = pd.Series([interactome_dict[l].graph['num_small_molecules'] for l in LEVELS], index=LEVELS)

    return num_interactions, num_entities, num_small_molecules


if __name__ == '__main__':
    # print(f"Working directory: {os.getcwd()}")

    graphs_path = "../../resources/Reactome/"
    interactomes = {l: create_interactome(l, True, graphs_path) for l in config.LEVELS}
    print(f"Indexing vertices")
    save_interactomes_with_indexed_vertices(interactomes, graphs_path)
