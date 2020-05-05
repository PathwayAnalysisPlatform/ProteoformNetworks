import re

import networkx as nx

from Python.network_topology_queries import get_reactions_and_participants_by_pathway


def add_edges_from_product(G, c1, c2, color, reaction, verbose=False):
    for i in c1:
        for j in c2:
            if i != j:
                G.add_edge(i, j, edge_color=color, reaction=reaction)
                if verbose:
                    print(f"Added edge from: {i} to {j}")


def add_edges(G, inputs, outputs, catalysts, regulators, reaction, verbose=False):
    print(f"\n\nReaction: {reaction}:\nInputs: {inputs}\nCatalysts: {catalysts}\nOutputs: {outputs}\nRegulators: {regulators}")
    add_edges_from_product(G, inputs, outputs, "black", reaction, True)
    add_edges_from_product(G, catalysts, outputs, "green", reaction, True)
    add_edges_from_product(G, regulators, outputs, "red", reaction, True)


def create_graph(df, directed=False, showSmallMolecules=True):
    """
    Converts a set of reactions with its participants into a graph

    :param df: Pandas dataframe with reactions and its participants
    :param directed: True creates a directed graph, False an undirected graph
    :param showSmallMolecules: False shows only EntityWithAccessionSequence participants
    :return: networkx graph with the interaction network representation of the reactions
    """

    df['Name'] = df['Name'].apply(lambda x: re.sub("\s*\[[\s\w]*\]\s*", '', x))

    # TODO Connect complex neighbours

    # Add all nodes
    G = nx.Graph()

    reaction = df.iloc[1]['Reaction']
    inputs = set()
    outputs = set()
    catalysts = set()
    regulators = set()
    for index, participant in df.iterrows():
        G.add_node(participant['Id'],
                   type=participant['Type'],
                   role=participant['Role'],
                   reaction=participant['Reaction'])
        if reaction != participant["Reaction"]:
            add_edges(G, inputs, outputs, catalysts, regulators, reaction, True)
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
    add_edges(G, inputs, outputs, catalysts, regulators, reaction, True)

    print(f"Added {len(G.nodes)} nodes to the graph.")
    print(f"Added {len(G.edges)} edges to the graph.")
    return G


df = get_reactions_and_participants_by_pathway("R-HSA-9634600", level="proteoforms")
# print(df.columns)
print(df[['Id', 'Role']])
G = create_graph(df)
print(G.nodes)
print(G.edges)

