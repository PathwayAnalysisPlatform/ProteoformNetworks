from collections import namedtuple

import networkx as nx
import pandas as pd
import os

from config import get_entity_color, LEVELS
from network_topology_queries import get_reaction_participants_by_pathway, get_complex_components_by_pathway, \
    get_pathway_name, get_reactions, get_complexes, get_reaction_participants_by_reaction, \
    get_complex_components_by_complex

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


def add_edges(g, inputs, outputs, catalysts, regulators, reaction, verbose):
    if verbose:
        print(
            f"\n\nReaction: {reaction}:\nInputs: {inputs}\nCatalysts: {catalysts}\nOutputs: {outputs}\nRegulators: {regulators}")
    add_edges_from_product(g, inputs, outputs, verbose)
    add_edges_from_product(g, catalysts, outputs, verbose)
    add_edges_from_product(g, regulators, outputs, verbose)


def connect_reaction_participants(g, df, level, verbose):
    if len(df) == 0:
        return g
    reaction = df.iloc[0]['Reaction']
    inputs = set()
    outputs = set()
    catalysts = set()
    regulators = set()
    for index, participant in df.iterrows():
        g.add_node(participant['Id'],
                   type=(participant['Type'] if participant["Type"] == "SimpleEntity" else level[:-1]),
                   role=participant['Role'],
                   reaction=participant['Reaction'],
                   color=get_entity_color(participant["Type"], level))
        if reaction != participant["Reaction"]:
            add_edges(g, inputs, outputs, catalysts, regulators, reaction, verbose)
            reaction = participant["Reaction"]
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
    add_edges(g, inputs, outputs, catalysts, regulators, reaction, verbose)

    if (verbose):
        print(f"From reactions, added {len(g.nodes)} nodes to the graph.")
        print(f"From reactions, added {len(g.edges)} edges to the graph.")


def connect_complex_components(g, df, level="proteins", verbose=True):
    if len(df) == 0:
        return g
    complex = df.iloc[0]['Complex']
    components = set()

    for index, record in df.iterrows():
        if not g.has_node(record['Id']):
            g.add_node(record['Id'],
                       type=(record['Type'] if record["Type"] == "SimpleEntity" else level[:-1]),
                       complex=record['Complex'],
                       color=get_entity_color(record["Type"], level))

        if complex != record['Complex']:
            add_edges_from_product(g, components, components, verbose=verbose)
            complex = record['Complex']
            components = set()
        components.add(record['Id'])
    complex = df.iloc[-1]['Complex']
    add_edges_from_product(g, components, components, verbose=verbose)

    if (verbose):
        print(f"From complexes, added {len(g.nodes)} nodes to the graph.")
        print(f"From complexes, added {len(g.edges)} edges to the graph.")



def create_graph(pathway, level, showSmallMolecules, graphs_path="", verbose=False):
    """
    Converts a set of reactions with its participants into a graph

    :param df_reactions: Pandas dataframe with reactions and its participants
    :param df_complexes: Pandas dataframe with complexes and its components
    :param showSmallMolecules: False shows only EntityWithAccessionSequence participants
    :return: networkx graph with the interaction network representation of the reactions
    """


    G = nx.Graph()
    name = get_pathway_name(pathway)
    if len(name) == 0:
        return G

    G.graph["stId"] = pathway
    G.graph["level"] = level
    G.graph["showSmallMolecules"] = showSmallMolecules

    df_reactions = get_reaction_participants_by_pathway(pathway, level, showSmallMolecules, verbose)
    df_complexes = get_complex_components_by_pathway(pathway, level, showSmallMolecules, verbose)

    connect_reaction_participants(G, df_reactions, level, verbose)
    connect_complex_components(G, df_complexes, level, verbose)

    if verbose:
        print(f"Storing network")

    out_file = pathway + "_" + level
    if not showSmallMolecules:
        out_file += "_no_small_molecules"
    out_file += "_edge_list"

    if len(str(graphs_path)) > 0:
        if not os.path.exists(graphs_path):
            os.mkdir(graphs_path)
    fh = open(os.path.join(str(graphs_path), out_file), 'wb')
    nx.write_edgelist(G, fh, data=True)

    print(f"Created graph {pathway} - {level} - {showSmallMolecules}", flush=True)
    fh.close()
    return G


def get_pathway_graphs(pathway, graphs_path):
    name = get_pathway_name(pathway)
    if len(name) == 0:
        return [nx.Graph(), nx.Graph(), nx.Graph()], [nx.Graph(), nx.Graph(), nx.Graph()]
    else:
        graphs = [create_graph(pathway, level, True, graphs_path) for level in LEVELS]
        graphs_no_small_molecules = [create_graph(pathway, level, False, graphs_path) for level in LEVELS]
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


def create_full_graph(level, showSmallMolecules, graphs_path="", verbose=False):
    G = nx.Graph()
    G.graph["level"] = level
    G.graph["showSmallMolecules"] = showSmallMolecules

    reactions = get_reactions()
    complexes = get_complexes()

    reactions = reactions[:5]
    complexes = complexes[:5]

    index = 0
    for reaction in reactions['stId']:
        df_participants = get_reaction_participants_by_reaction(reaction, level, showSmallMolecules, False)
        connect_reaction_participants(G, df_participants, level, True)
        index += 1
        if verbose:
            print(
                f"Added reaction {reaction}\t{index}/{len(reactions)}\t\tNodes: {len(G.nodes)}, \t\t Edges: {len(G.edges)}",
                flush=True)

    index = 0
    for complex in complexes['stId']:
        df_components = get_complex_components_by_complex(complex, level, showSmallMolecules, False)
        connect_complex_components(G, df_components, level, True)
        index += 1
        if verbose:
            print(
                f"Added complex {complex}\t{index}/{len(complexes)}\t\tNodes: {len(G.nodes)}, \t\t Edges: {len(G.edges)}",
                flush=True)

    if len(graphs_path) > 0:
        Path(graphs_path).mkdir(parents=True, exist_ok=True)
    out_file = graphs_path + level
    if not showSmallMolecules:
        out_file += "_no_small_molecules"
    out_file += "_edge_list.csv"

    if verbose:
        print(f"Storing network: {out_file}")

    fh = open(out_file, 'wb')
    nx.write_edgelist(G, out_file, data=True)


    return G

# create_graph("R-HSA-9634600", "proteins", True, "hello/")

#
# print("Started creating the graph:")
# create_full_graph("genes", True, "", True)
# print("Finished!!")
