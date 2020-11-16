import codecs
import json
import os
from collections import namedtuple
from pathlib import Path

import networkx as nx
import pandas as pd

import config
from config import LEVELS, get_entity_color
from queries import get_reaction_participants_by_pathway, get_complex_components_by_pathway, \
    get_pathway_name

Pathway_graphs = namedtuple("Graphs",
                            ['genes', "genes_no_small_molecules", "proteins", "proteins_no_small_molecules",
                             "proteoforms", "proteoforms_no_small_molecules"])


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


def add_edges_reaction_participants(G, df, v=False):
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
    for index, participant in df.iterrows():
        unique_id = 'sm_' + participant['Id'] if participant['Type'] == 'SimpleEntity' else participant['Id']
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


def add_edges_complex_components(G, df, v=False):
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
    complex = df.iloc[0]['Complex']
    for index, record in df.iterrows():
        unique_id = 'sm_' + record['Id'] if record['Type'] == 'SimpleEntity' else record['Id']
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


def add_nodes(G, df, level):
    """
    Add the nodes represented by a dataframe records to the graph G

    :param G: Networkx graph to add the nodes
    :param df: Dataframe with the records resulting from the cypher query for participants of reactions or components
    of complexes
    :param level: genes, proteins or proteoforms. To handle the renaming of some fields accordingly
    :return: The same graph G with the added nodes
    """
    for index, entity in df.iterrows():

        # Each node has an id, label, type, entity_color, role, reactions, pathways, complexes.
        # id: unique identifier across all levels and entity types. For example, to avoid duplicated nodes when a gene
        # or a small molecule are called the same.
        # label: gene name, small molecule name, protein name, and for proteoforms the protein name with possible manual
        # annotation of modifications.
        # print(f"Adding node {index} of entity {entity}")
        unique_id = 'sm_' + entity['Id'] if entity['Type'] == 'SimpleEntity' else entity['Id']
        if unique_id not in G.nodes:
            G.add_node(unique_id,
                       label=entity['Id'],
                       type=(entity['Type'] if entity['Type'] == 'SimpleEntity' else level),
                       entity_color=get_entity_color(entity['Type'], level),
                       roles=set(),
                       reactions=set(),
                       pathways=set(),
                       complexes=set()
                       )
            if entity['Type'] == 'SimpleEntity':
                G.graph['num_small_molecules'] += 1
            else:
                G.graph['num_' + level] += 1

        if 'Complex' in df.columns:
            G.nodes[unique_id]['complexes'].add(entity['Complex'])
        elif 'Role' in df.columns:
            G.nodes[unique_id]['roles'].add(entity['Role'])
            G.nodes[unique_id]['reactions'].add(entity['Reaction'])
            G.nodes[unique_id]['pathways'].add(entity['Pathway'])

        if (G.nodes[unique_id]['type'] == "SimpleEntity") or G.nodes[entity['Id']]['type'] == "genes":
            G.nodes[unique_id]['prevId'] = entity['Id']
        elif G.nodes[unique_id]['type'] == "proteins" or G.nodes[entity['Id']]['type'] == "proteoforms":
            G.nodes[unique_id]['prevId'] = entity['PrevId']

    print(f"{level} level - small molecules: {G.graph['num_small_molecules']}")
    print(f"{level}: {G.graph['num_' + level]}")


def read_graph(json_file):
    """
    Read interaction network in the json file into a networkx graph instance

    :param json_file:
    :return: networkx Graph
    """
    G = nx.Graph()
    with open(json_file) as f:
        data = json.load(f)
        G = nx.json_graph.node_link_graph(data)

    return G


def save_graph(G, level, graphs_path="", label=""):
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
    print(f"Storing network level: {level} at directory: {graphs_path}")

    if len(label) > 0:
        label = label + "_"

    json_file = Path(graphs_path + label + level + ".json")
    vertices_file = Path(graphs_path + label + level + "_vertices.tsv")
    edges_file = Path(graphs_path + label + level + "_interactions.tsv")
    small_molecules_file = Path(graphs_path + label + level + "_small_molecules.tsv")

    if len(str(graphs_path)) > 0:
        if not os.path.exists(graphs_path):
            Path(graphs_path).mkdir(parents=True, exist_ok=True)

    with open(json_file, 'w') as outfile:
        data = nx.json_graph.node_link_data(G)
        json.dump(data, outfile)
    print(f"Created json file.")

    with open(edges_file, 'wb') as fh:
        nx.write_edgelist(G, fh)
    print(f"Created edges file for {level}")

    with codecs.open(vertices_file, 'w', "utf-8") as fh:
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

    if level == "proteoforms":
        mapping_file = graphs_path + config.MAPPING_FILE.replace('level', 'proteoforms')
        with codecs.open(mapping_file, 'w', "utf-8") as map_file:
            for n, t in G.nodes(data='type'):
                if t != "SimpleEntity":
                    map_file.write(f"{G.nodes[n]['prevId']}\t{n}\n")
        print(f"Created mapping proteins to proteoforms file.")


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
        save_graph(G, level, graphs_path, label=pathway + ("_no_sm" if not sm else ""))

    print(f"Created graph {pathway} - {level} - {sm}", flush=True)
    return G


def create_pathway_graphs(pathway, graphs_path="../../reports/pathways/", v=False):
    """
    Creates interactomes for the pathway in all three levels, with and wihout small molecules

    If pathway does not exists, then returns empty interactomes.
    :param pathway: Pathway stId string
    :param graphs_path:
    :return: Get two lists of 3 pathways.
    """

    name = get_pathway_name(pathway)
    if len(name) == 0:
        return [nx.Graph(), nx.Graph(), nx.Graph()], [nx.Graph(), nx.Graph(), nx.Graph()]
    else:
        graphs = [create_graph(pathway, level, True, graphs_path, v, True) for level in LEVELS]
        graphs_no_small_molecules = [create_graph(pathway, level, False, graphs_path, v, True) for level in LEVELS]
        return graphs, graphs_no_small_molecules


def create_graphs(pathways, graphs_path="../../reports/pathways/"):
    """
    Create the interaction networks of each pathway of the arguments
    :param pathways: pandas dataframe with column "stId"
    :param graphs_path: location for output files
    :return: List of namedtuples()
    """
    result = list()
    if type(pathways) == pd.DataFrame and len(pathways) > 0:
        for pathway in pathways['stId']:
            name = get_pathway_name(pathway)
            name = name.iloc[0]['Name']
            print(f"Creating networks for pathway {pathway}")
            graphs, graphs_no_small_molecules = create_pathway_graphs(pathway, graphs_path)

            gs = Pathway_graphs(graphs[0], graphs_no_small_molecules[0], graphs[1], graphs_no_small_molecules[1],
                                graphs[2],
                                graphs_no_small_molecules[2])
            result.append(gs)
    return result


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
