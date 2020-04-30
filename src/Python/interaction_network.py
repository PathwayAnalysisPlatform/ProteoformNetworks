import re

import networkx as nx

from Python.network_topology_queries import get_reactions_and_participants_by_pathway


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

    reaction = df.loc[1]['Reaction']
    print(f"First reaction is: {reaction}")
    inputs = set()
    outputs = set()
    catalysts = set()
    regulators = set()
    for index, participant in df.iterrows():
        G.add_node(participant['Name'],
                   type=participant['Type'],
                   role=participant['Role'],
                   reaction=participant['Reaction'])
        if reaction != participant["Reaction"]:
            # print(f"Reaction: {reaction}:\n")
            # print(f"Inputs: {inputs}")
            # print(f"Catalysts: {catalysts}")
            # print(f"Outputs: {outputs}")
            # print(f"Regulators: {regulators}")
            for i in inputs:
                for o in outputs:
                    G.add_edge(i, o, edge_color='black')
                    # print(f"Added edge from: {i} to {o}")
            for c in catalysts:
                for o in outputs:
                    G.add_edge(c, o, edge_color='green')
                    # print(f"Added edge from: {c} to {o}")
            for r in regulators:
                for o in outputs:
                    G.add_edge(r, o, edge_color='red')
                    # print(f"Added edge from: {r} to {o}")
            reaction = participant["Reaction"]
            inputs = set()
            outputs = set()
            catalysts = set()
            regulators = set()
        else:
            if participant['Role'] == 'input':
                inputs.add(participant['Name'])
            elif participant['Role'] == 'output':
                outputs.add(participant['Name'])
            elif participant['Role'] == 'regulatedBy':
                regulators.add(participant['Name'])
            elif participant['Role'] == 'catalystActivity':
                catalysts.add(participant['Name'])

    print(f"Added {len(G.nodes)} nodes to the graph.")
    print(f"Added {len(G.edges)} edges to the graph.")
    return G


# df = get_reactions_and_participants_by_pathway("R-HSA-70171")
# G = create_graph(df)
# print(G.nodes)
# print(G.edges)
