import os
from collections import namedtuple

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
    if len(df) == 0:
        return G
    pathway = df.iloc[0]['Pathway']
    reaction = df.iloc[0]['Reaction']
    inputs = set()
    outputs = set()
    catalysts = set()
    regulators = set()
    for index, participant in df.iterrows():
        if not G.has_node(participant['Id']):
            G.add_node(participant['Id'],
                       type=(participant['Type'] if participant['Type'] == 'SimpleEntity' else level),
                       entity_color=get_entity_color(participant['Type'], level),
                       roles=set(), reactions=set(), pathways=set(), complexes=set())
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
        if not G.has_node(component['Id']):
            G.add_node(component['Id'],
                       type=(component['Type'] if component['Type'] == 'SimpleEntity' else level),
                       entity_color=get_entity_color(component['Type'], level),
                       roles=set(), reactions=set(), pathways=set(), complexes=set())
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


def create_graph(pathway, level, sm, graphs_path="", v=False):
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

    df_reactions = get_reaction_participants_by_pathway(pathway, level, sm, v)
    df_complexes = get_complex_components_by_pathway(pathway, level, sm, v)

    connect_reaction_participants(G, df_reactions, level, v)
    connect_complex_components(G, df_complexes, level, v)

    if v:
        print(f"Storing network")

    out_file = pathway + "_" + level
    if not sm:
        out_file += "_no_small_molecules"
    out_file += "_edge_list"

    if len(str(graphs_path)) > 0:
        if not os.path.exists(graphs_path):
            os.mkdir(graphs_path)
    fh = open(os.path.join(str(graphs_path), out_file), 'wb')
    nx.write_edgelist(G, fh, data=True)

    print(f"Created graph {pathway} - {level} - {sm}", flush=True)
    fh.close()
    return G


def get_pathway_graphs(pathway, graphs_path):
    name = get_pathway_name(pathway)
    if len(name) == 0:
        return [nx.Graph(), nx.Graph(), nx.Graph()], [nx.Graph(), nx.Graph(), nx.Graph()]
    else:
        graphs = [create_graph(pathway, level, True, graphs_path, True) for level in LEVELS]
        graphs_no_small_molecules = [create_graph(pathway, level, False, graphs_path, True) for level in LEVELS]
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
            graphs, graphs_no_small_molecules = get_pathway_graphs(pathway, graphs_path)

            gs = Pathway_graphs(graphs[0], graphs_no_small_molecules[0], graphs[1], graphs_no_small_molecules[1],
                                graphs[2],
                                graphs_no_small_molecules[2])
            result.append(gs)
    return result


if __name__ == '__main__':
    create_graph("R-HSA-8981607", "proteins", True)
