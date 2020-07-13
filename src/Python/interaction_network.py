import codecs
import os
from collections import namedtuple
from pathlib import Path

import networkx as nx
import pandas as pd

from config import LEVELS, get_entity_color
from network_topology_queries import get_reaction_participants_by_pathway, get_complex_components_by_pathway, \
    get_pathway_name

Pathway_graphs = namedtuple("Graphs",
                            ['genes', "genes_no_small_molecules", "proteins", "proteins_no_small_molecules",
                             "proteoforms", "proteoforms_no_small_molecules"])


def add_edges_from_product(G, c1, c2, verbose=False):
    for i in c1:
        for j in c2:
            if i != j:
                G.add_edge(i, j)
                if verbose:
                    print(f"Added edge from: {i} to {j}")


def add_edges(G, inputs, outputs, catalysts, regulators, reaction, verbose):
    if verbose:
        print(
            f"\n\nReaction: {reaction}:\nInputs: {inputs}\nCatalysts: {catalysts}\nOutputs: {outputs}\nRegulators: {regulators}")
    add_edges_from_product(G, inputs, outputs, verbose)
    add_edges_from_product(G, catalysts, outputs, verbose)
    add_edges_from_product(G, regulators, outputs, verbose)


def connect_reaction_participants(G, df, level, v):
    """
    G must contain all the nodes specified by the edges
    :param G:
    :param df:
    :param level:
    :param v:
    :return:
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
        G.nodes[participant['Id']]['roles'].add(participant['Role'])
        G.nodes[participant['Id']]['reactions'].add(participant['Reaction'])
        G.nodes[participant['Id']]['pathways'].add(participant['Pathway'])
        if reaction != participant['Reaction'] or pathway != participant['Pathway']:
            add_edges(G, inputs, outputs, catalysts, regulators, reaction, v)
            reaction = participant['Reaction']
            pathway = participant['Pathway']
            inputs = set()
            outputs = set()
            catalysts = set()
            regulators = set()
        if participant['Role'] == 'input':
            inputs.add(participant['Id'])
        elif participant['Role'] == 'output':
            outputs.add(participant['Id'])
        elif participant['Role'] == 'regulatedBy':
            regulators.add(participant['Id'])
        elif participant['Role'] == 'catalystActivity':
            catalysts.add(participant['Id'])
    reaction = df.iloc[-1]['Reaction']
    pathway = df.iloc[-1]['Pathway']
    add_edges(G, inputs, outputs, catalysts, regulators, reaction, v)

    if (v):
        print(f"From reactions, added {len(G.nodes)} nodes to the graph.")
        print(f"From reactions, added {len(G.edges)} edges to the graph.")


def connect_complex_components(G, df, level="proteins", v=True):
    if len(df) == 0:
        return G

    components = set()
    complex = df.iloc[0]['Complex']
    for index, component in df.iterrows():
        G.nodes[component['Id']]['complexes'].add(component['Complex'])

        if complex != component['Complex']:
            add_edges_from_product(G, components, components, verbose=v)
            complex = component['Complex']
            components = set()
        components.add(component['Id'])
    complex = df.iloc[-1]['Complex']
    add_edges_from_product(G, components, components, verbose=v)

    if (v):
        print(f"From complexes, added {len(G.nodes)} nodes to the graph.")
        print(f"From complexes, added {len(G.edges)} edges to the graph.")


def add_nodes(G, df, level):
    for index, entity in df.iterrows():
        # print(f"Adding node {index} of entity {entity}")
        if entity['Id'] not in G.nodes:
            G.add_node(entity['Id'],
                       id=entity['Id'],
                       type=(entity['Type'] if entity['Type'] == 'SimpleEntity' else level),
                       entity_color=get_entity_color(entity['Type'], level),
                       roles=set(),
                       reactions=set(),
                       pathways=set(),
                       complexes=set(),
                       )
            if entity['Type'] == 'SimpleEntity':
                G.graph['num_small_molecules'] += 1
            else:
                G.graph['num_' + level] += 1

        if 'Complex' in df.columns:
            G.nodes[entity['Id']]['complexes'].add(entity['Complex'])
        elif 'Role' in df.columns:
            G.nodes[entity['Id']]['roles'].add(entity['Role'])
            G.nodes[entity['Id']]['reactions'].add(entity['Reaction'])
            G.nodes[entity['Id']]['pathways'].add(entity['Pathway'])

        if level == 'proteins':
            if (G.nodes[entity['Id']]['type'] == "SimpleEntity"):
                G.nodes[entity['Id']]['prevId'] = entity['Id']
            else:
                G.nodes[entity['Id']]['prevId'] = entity['PrevId']
        elif level == 'proteoforms':
            if (G.nodes[entity['Id']]['type'] == "SimpleEntity"):
                G.nodes[entity['Id']]['prevId'] = entity['Id']
            else:
                G.nodes[entity['Id']]['prevId'] = entity['PrevId']


def read_graph(graphs_path, file_preffix, level):
    graphs_path = Path(graphs_path)

    vertices_file = graphs_path / ((file_preffix + "_" if len(file_preffix) > 0 else "") + level + "_vertices.tsv")
    edges_file = graphs_path / ((file_preffix + "_" if len(file_preffix) > 0 else "") + level + "_interactions.tsv")
    small_molecules_file = graphs_path / ((file_preffix + "_" if len(file_preffix) > 0 else "") + level + "_small_molecules_vertices.tsv")

    if not Path(vertices_file).exists() or not Path(edges_file).exists() or not Path(small_molecules_file).exists():
        print("Graph files not found:")
        print(f"-- {vertices_file}")
        print(f"-- {edges_file}")
        print(f"-- {small_molecules_file}")

    print(f"    Reading {edges_file}")
    G = nx.read_adjlist(edges_file.as_posix())
    G.graph["level"] = level
    G.graph['num_' + level] = 0
    G.graph['num_small_molecules'] = 0

    print(f"    Reading {vertices_file}")
    with vertices_file.open() as fh:
        entity = fh.readline()
        while entity:
            attrs = {entity.strip(): {'type': level}}
            nx.set_node_attributes(G, attrs)
            G.graph['num_' + level] += 1
            entity = fh.readline()

    print(f"    Reading {small_molecules_file}")
    with small_molecules_file.open() as fh:
        small_molecule = fh.readline()
        while small_molecule:
            attrs = {small_molecule.strip(): {'type': "SimpleEntity"}}
            nx.set_node_attributes(G, attrs)
            G.graph['num_small_molecules'] += 1
            small_molecule = fh.readline()

    print(f"        Graph edges: {G.size()}")
    print(f"        Graph nodes: {len(G.nodes)}")
    print(f"        Graph {level} nodes: {G.graph['num_' + level]}")
    print(f"        Graph small molecule nodes: {G.graph['num_small_molecules']}")
    return G


def save_graph(G, level, graphs_path="", label=""):
    print(f"Storing network level: {level}")

    graphs_path = Path(graphs_path)

    if len(label) > 0:
        label = label + "_"

    vertices_file = graphs_path / (label + level + "_vertices.tsv")
    edges_file = graphs_path / (label + level + "_interactions.tsv")
    small_molecules_file = graphs_path / (label + level + "_small_molecules_vertices.tsv")

    if len(str(graphs_path)) > 0:
        if not os.path.exists(graphs_path):
            Path(graphs_path).mkdir(parents=True, exist_ok=True)
    fh = open(edges_file, 'wb')
    nx.write_adjlist(G, fh, delimiter="\t")

    print(f"Created edges file for {level}")
    fh.close()

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

    connect_reaction_participants(G, df_reactions, level, v)
    connect_complex_components(G, df_complexes, level, v)

    if save:
        save_graph(G, level, graphs_path, label=pathway + ("_no_sm" if not sm else ""))

    print(f"Created graph {pathway} - {level} - {sm}", flush=True)
    return G


def create_pathway_graphs(pathway, graphs_path="../../reports/pathways/", v=False):
    """
    Creates graphs for the pathway in all three levels, with and wihout small molecules

    If pathway does not exists, then returns empty graphs.
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
    # - Vertices: Composition of nodes in all graphs
    #    - Merge: sets of Reactions, Pathways, Complexes
    #    - Copy value of: Id, Type, Entity Color
    # - Edges: Composition of edges in all graphs
    full_graph = nx.compose_all(graphs)  # Add all nodes setting  Id, Type, and entity_color

    for graph in graphs:
        for node in graph.nodes:
            full_graph.nodes[node]['reactions'].update(graph.nodes[node]['reactions'])
            full_graph.nodes[node]['pathways'].update(graph.nodes[node]['pathways'])
            full_graph.nodes[node]['roles'].update(graph.nodes[node]['roles'])
            full_graph.nodes[node]['complexes'].update(graph.nodes[node]['complexes'])


if __name__ == '__main__':
    create_graph("R-HSA-8981607", "proteins", True, save=True)

    G = read_graph("", "R-HSA-8981607", "proteins")
    print(G.nodes)

    G = read_graph("../../resources/Reactome/", "", "genes")
    print(list(G)[:10])