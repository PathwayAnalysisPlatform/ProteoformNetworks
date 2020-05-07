import networkx as nx

from config import LEVELS_COLOR, get_entity_color
from network_topology_queries import get_reaction_participants_by_pathway, get_complex_components_by_pathway


def add_edges_from_product(G, c1, c2, color, reaction, complex, verbose=False):
    for i in c1:
        for j in c2:
            if i != j:
                G.add_edge(i, j, edge_color=color, reaction=reaction, complex=complex)
                if verbose:
                    print(f"Added edge from: {i} to {j}")


def add_edges(g, inputs, outputs, catalysts, regulators, reaction, verbose=True):
    if verbose:
        print(
            f"\n\nReaction: {reaction}:\nInputs: {inputs}\nCatalysts: {catalysts}\nOutputs: {outputs}\nRegulators: {regulators}")
    add_edges_from_product(g, inputs, outputs, "black", reaction, "", True)
    add_edges_from_product(g, catalysts, outputs, "green", reaction, "", True)
    add_edges_from_product(g, regulators, outputs, "red", reaction, "", True)


def connect_reaction_participants(g, df, directed=False, level="proteins", verbose=True):
    reaction = df.iloc[1]['Reaction']
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
            add_edges(g, inputs, outputs, catalysts, regulators, reaction)
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
    add_edges(g, inputs, outputs, catalysts, regulators, reaction, True)

    if (verbose):
        print(f"From reactions, added {len(g.nodes)} nodes to the graph.")
        print(f"From reactions, added {len(g.edges)} edges to the graph.")


def connect_complex_components(g, df, directed=False, level="proteins", verbose=True):
    complex = df.iloc[1]['Complex']
    components = set()

    for index, record in df.iterrows():
        if not g.has_node(record['Id']):
            g.add_node(record['Id'],
                       type=(record['Type'] if record["Type"] == "SimpleEntity" else level[:-1]),
                       complex=record['Complex'],
                       color=get_entity_color(record["Type"], level))

        if complex != record['Complex']:
            add_edges_from_product(g, components, components, "blue", "", complex)
            components = set()
        components.add(record['Id'])
    complex = df.iloc[-1]['Complex']
    add_edges_from_product(g, components, components, "blue", "", complex)

    if (verbose):
        print(f"From complexes, added {len(g.nodes)} nodes to the graph.")
        print(f"From complexes, added {len(g.edges)} edges to the graph.")


def create_graph(pathway, level="proteins", directed=False, showSmallMolecules=True, verbose=True):
    """
    Converts a set of reactions with its participants into a graph

    :param df_reactions: Pandas dataframe with reactions and its participants
    :param df_complexes: Pandas dataframe with complexes and its components
    :param directed: True creates a directed graph, False an undirected graph
    :param showSmallMolecules: False shows only EntityWithAccessionSequence participants
    :return: networkx graph with the interaction network representation of the reactions
    """

    G = nx.Graph()
    G.graph["stId"] = pathway
    G.graph["level"] = level
    G.graph["showSmallMolecules"] = showSmallMolecules

    df_reactions = get_reaction_participants_by_pathway(pathway, showSmallMolecules=showSmallMolecules,
                                                        level=level, verbose=verbose)
    df_complexes = get_complex_components_by_pathway(pathway, showSmallMolecules=showSmallMolecules,
                                                     level=level, verbose=verbose)

    connect_reaction_participants(G, df_reactions, directed, level, verbose)
    connect_complex_components(G, df_complexes, directed, level, verbose)

    return G

# df = get_reaction_participants_by_pathway("R-HSA-9634600", level="proteoforms")
# print(df_reactions.columns)
# print(df[['Id', 'Role']])
# g = create_graph(df)
# print(g.nodes)
# print(g.edges)
