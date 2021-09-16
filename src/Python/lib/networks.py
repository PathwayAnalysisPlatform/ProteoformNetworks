import codecs
import json
import os
from bisect import bisect_left
from pathlib import Path

import networkx as nx
import pandas as pd

import config
from config import no_sm, with_sm, with_unique_sm, sm, LEVELS
from lib.graph_database_access import get_pathway_name, get_participants_by_pathway, get_components_by_pathway, \
    get_pathways, get_participants, get_components


def create_pathway_interaction_network(pathway, level, method, out_path=""):
    if level not in LEVELS:
        raise Exception(f"There is no {level} level.")

    filename = get_json_filename(level, method, out_path, pathway)
    if not Path(filename).exists():
        print(f"    * Creating network {filename}")
        participants = {level: get_participants_by_pathway(pathway, level, out_path) for level in [level, sm]}
        components = {level: get_components_by_pathway(pathway, level, out_path) for level in [level, sm]}
        create_interaction_network(level, method, participants, components, out_path, pathway)
    # else:
    #     print(f"    * Network {filename} already exists..")
    G = read_graph(filename)

    missing_bridges = False
    missing_articulations = False
    if any("Bridge" not in G.edges[edge] for edge in G.edges):
        set_bridges(G)
        missing_bridges = True

    if any("Articulation Point" not in G.nodes[node] for node in G.nodes):
        set_articulation_points(G)

    if missing_articulations or missing_bridges:
        update_json_file(G, level, method, out_path, pathway)

    return G


def create_pathway_interaction_networks(pathway, out_path=""):
    """
    Creates interaction networks for a pathway in all three levels, with the 3 contruction methods

    If pathway does not exists, then returns empty networks.
    :param pathway: Pathway stId string
    :param out_path:
    :return: Get dictionary {method: lists of 3 pathways}.
    """

    name = get_pathway_name(pathway)
    if len(name) == 0:
        print(f"Pathway {pathway} does not exist")
        return {m: {l: nx.Graph() for l in LEVELS} for m in config.METHODS}
    else:
        print(f"-- Creating interaction networks for pathway {pathway}")
        return {
            m: {level: create_pathway_interaction_network(pathway, level, m, out_path) for level in LEVELS}
            for m in config.METHODS
        }


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


def get_multiindex():
    arrays = [
        [*(["Not Included"] * 3), *(["Included"] * 3), *(["Reaction-Unique Included"] * 3)], [*(LEVELS * 3)]
    ]
    tuples = list(zip(*arrays))
    index = pd.MultiIndex.from_tuples(tuples, names=["Small Molecules", "Entity Level"])
    return index


def print_interactome_details(g):
    print(f"Graph for {g.graph['level']} ")
    print(f"Graph edges: {g.size()}")
    print(f"Graph nodes: {len(g.nodes)}")
    print(f"Graph {g.graph['level']} nodes: {g.graph['num_entities']}")
    print(f"Graph small molecule nodes: {g.graph['num_small_molecules']}")
    print("\n***********************\n\n")


def get_json_filename(level, method, out_path="", label=""):
    if not isinstance(out_path, str):
        out_path = out_path.dirname
    if len(out_path) > 0:
        if out_path[-1] != '/':
            out_path += "/"
    if len(label) > 0:
        return Path(out_path + label + "_" + level + "_" + method + ".json")
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
    sm_id_column = 'Id' if not method == with_unique_sm else 'UniqueId'

    for index, row in df.iterrows():

        # Each node has an id, label, type, entity_color, role, reactions, pathways, complexes.
        # id: unique identifier across all levels and entity types. For example, to avoid duplicated nodes when a gene
        # or a small molecule are called the same.
        # label: gene name, small molecule name, protein name, and for proteoforms the protein name with possible manual
        # annotation of modifications.
        # print(f"Adding node {index} of entity {entity}")
        unique_id = row[sm_id_column] if row['Type'] == 'SimpleEntity' else row['Id']

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
            if row['Type'].startswith("S"):
                G.graph['num_small_molecules'] += 1
            else:
                G.graph['num_entities'] += 1

        if 'Complex' in df.columns:
            G.nodes[unique_id]['complexes'].add(row['Complex'])
        elif 'Role' in df.columns:
            G.nodes[unique_id]['roles'].add(row['Role'])
            G.nodes[unique_id]['reactions'].add(row['Reaction'])
            G.nodes[unique_id]['pathways'].add(row['Pathway'])

        if G.nodes[unique_id]['type'].startswith("S"):
            G.nodes[unique_id]['prevId'] = row[sm_id_column]
        elif G.nodes[row['Id']]['type'].startswith("g"):
            G.nodes[unique_id]['prevId'] = row['Id']
        else:  # G.nodes[unique_id]['type'] == "proteins" or G.nodes[row['Id']]['type'] == "proteoforms":
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


def add_edges_from_product(G, c1, c2, v=False):
    for i in c1:
        for j in c2:
            if i != j:
                G.add_edge(i, j)
                if v:
                    print(f"Added edge from: {i} to {j}")


def add_edges_for_reaction(G, inputs, outputs, catalysts, regulators, reaction, v=False):
    if v:
        print(
            f"\n\nReaction: {reaction}:\nInputs: {inputs}\nCatalysts: {catalysts}\nOutputs: {outputs}\nRegulators: {regulators}")
    add_edges_from_product(G, inputs, outputs, v)
    add_edges_from_product(G, catalysts, outputs, v)
    add_edges_from_product(G, regulators, outputs, v)


def add_edges_all_reaction_participants(G, df, method=with_sm, v=False):
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

    sm_id_column = 'Id' if not method == with_unique_sm else 'UniqueId'

    for index, participant in df.iterrows():
        unique_id = participant[sm_id_column] if participant['Type'] == 'SimpleEntity' else participant['Id']
        G.nodes[unique_id]['roles'].add(participant['Role'])
        G.nodes[unique_id]['reactions'].add(participant['Reaction'])
        G.nodes[unique_id]['pathways'].add(participant['Pathway'])
        if reaction != participant['Reaction'] or pathway != participant['Pathway']:
            add_edges_for_reaction(G, inputs, outputs, catalysts, regulators, reaction, v)
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
    add_edges_for_reaction(G, inputs, outputs, catalysts, regulators, reaction, v)

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

    sm_id_column = 'Id' if not method == with_unique_sm else 'UniqueId'

    complex = df.iloc[0]['Complex']
    for index, record in df.iterrows():
        unique_id = record[sm_id_column] if record['Type'] == 'SimpleEntity' else record['Id']
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


def update_json_file(G, level, method, out_path="", label=""):

    name_start = ""
    if len(label) > 0:
        name_start += label + "_"

    json_file = get_json_filename(level, method, out_path, label)
    # print(f"Updating file {json_file}")

    with open(json_file, 'w+') as outfile:
        data = nx.json_graph.node_link_data(G)
        json.dump(data, outfile)
    # print(f"Updated json file for {label}")

def save_interaction_network(G, level, method, out_path="", label=""):
    """
    Create the json file with all attributes.
    Creates a vertices file for the accessioned entities.
    Creates a vertices file for the simple entities.
    Creates an edges file for the interactions.
    If the level is proteins write the mapping from proteins to genes
    If the level is proteoforms write the mapping from proteins to proteoforms

    :param G: networkx file with the
    :param level: genes, proteins or protoeforms
    :param out_path: directory where to store the files
    :param label: string phrase for the file names
    :return: void
    """
    print(f"Storing network level: {level} at directory: {out_path}")

    if isinstance(out_path, str) and len(out_path) > 0:
        if out_path[-1] != '/':
            out_path += "/"

    name_start = ""
    if len(label) > 0:
        name_start += label + "_"

    name_start += level + "_" + method + "_"
    json_file = get_json_filename(level, method, out_path, label)
    accessioned_entities_file = Path(out_path + name_start + "accessioned_entities.tsv")
    interactions_file = Path(out_path + name_start + "interactions.tsv")
    small_molecules_file = Path(out_path + name_start + "small_molecules.tsv")

    if len(str(out_path)) > 0:
        if not os.path.exists(out_path):
            Path(out_path).mkdir(parents=True, exist_ok=True)

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
        mapping_file = out_path + config.MAPPING_FILE.replace('level', 'genes')
        with codecs.open(mapping_file, 'w', "utf-8") as map_file:
            for n, t in G.nodes(data='type'):
                if t != "SimpleEntity":
                    map_file.write(f"{n}\t{G.nodes[n]['prevId']}\n")
        print(f"Created mapping proteins to genes file.")

    if level == config.proteoforms:
        mapping_file = out_path + config.MAPPING_FILE.replace('level', 'proteoforms')
        with codecs.open(mapping_file, 'w', "utf-8") as map_file:
            for n, t in G.nodes(data='type'):
                if t != "SimpleEntity":
                    map_file.write(f"{G.nodes[n]['prevId']}\t{n}\n")
        print(f"Created mapping proteins to proteoforms file.")


def create_interaction_network(level, method, participants, components, out_path="", label=""):
    """
    Create interaction network with the participants and components provided as parameter.
    It does not care which level the participants are, it simply connects reaction parcitipants and complex participants.
    Creates a networkx graph and stores it in a json file.

    :param level: genes, proteins or proteoforms. This attribute is just to set it as graph property.
    :param method: "no_sm", "with_sm" or "with_unique_sm". This is just to set is as graph property.
    :param participants: pandas dataframe with the reaction participants
    :param components: pandas dataframe with the complex components
    :param out_path: path to directory to store the json file
    :return: The networkx interaction network
    """

    print("Creating interaction network...")

    G = nx.Graph()
    G.graph["level"] = level
    G.graph["method"] = method
    G.graph['num_entities'] = 0
    G.graph['num_small_molecules'] = 0

    if method == no_sm:
        add_nodes(G, pd.concat([participants[level], components[level]]))
        add_edges_all_reaction_participants(G, participants[level])
        add_edges_complex_components(G, components[level])
    elif method == with_sm or method == with_unique_sm:
        df_both_participants = pd.concat([participants[level], participants[sm]]).sort_values(["Pathway", "Reaction"])
        df_both_components = pd.concat([components[level], components[sm]]).sort_values(["Complex"])
        add_nodes(G, pd.concat([df_both_participants, df_both_components]), method)
        add_edges_all_reaction_participants(G, df_both_participants, method)
        add_edges_complex_components(G, df_both_components, method)
    else:
        raise Exception("No such method to create the interactome")

    for node in G.nodes:
        G.nodes[node]['complexes'] = list(G.nodes[node]['complexes'])
        G.nodes[node]['roles'] = list(G.nodes[node]['roles'])
        G.nodes[node]['reactions'] = list(G.nodes[node]['reactions'])
        G.nodes[node]['pathways'] = list(G.nodes[node]['pathways'])

    set_bridges(G)
    set_articulation_points(G)

    save_interaction_network(G, level, method, out_path, label)
    print(f"Finished creating interactome file for {level}-{method}")
    return G


def get_or_create_interaction_network(level, method, participants, components, out_path="", label=""):
    """
    Returns a networkx graph instance of the selected interaction network.
    Tries to read from a json file with the network. If the file does not exists it creates it.

    :param level: genes, proteins or proteoforms. This attribute is just to set it as graph property.
    :param method: "no_sm", "with_sm" or "with_unique_sm". This is just to set is as graph property.
    :param participants: pandas dataframe with the reaction participants
    :param components: pandas dataframe with the complex components
    :param out_path: path to directory to store the json file
    :param label: Any text to distinguish the graph file, ex. Pathway name
    :return: The networkx interaction network
    """

    filename = get_json_filename(level, method, out_path, label)

    if not Path(filename).exists():
        create_interaction_network(level, method, participants, components, out_path, label)
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


def save_indexed_vertices(indexed_vertices, out_path):
    vertices_file = Path(out_path + "interactome_vertices.tsv")
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


def save_interactomes_with_indexed_vertices(interactomes, out_path):
    """
    Save interactomes as three separate files with the list of edges using indexed vertices.

    :param interactomes: dict with a networkx graph for each level
    :param out_path:
    :return: void
    """
    # Gather the gene, protein, proteoform and small molecule vertices together and index all of them
    indexed_vertices, start_indexes, end_indexes = index_all_vertices(interactomes)
    assert end_indexes["SimpleEntity"] == len(indexed_vertices) - 1
    assert sorted(indexed_vertices[start_indexes["genes"]:end_indexes["genes"]])
    assert sorted(indexed_vertices[start_indexes["proteins"]:end_indexes["proteins"]])
    assert sorted(indexed_vertices[start_indexes["proteoforms"]:end_indexes["proteoforms"]])
    assert sorted(indexed_vertices[start_indexes["SimpleEntity"]:end_indexes["SimpleEntity"]])

    save_indexed_vertices(indexed_vertices, out_path)

    # Create three files with the interactomes
    save_edges_with_indexed_vertices(interactomes, out_path, indexed_vertices, start_indexes, end_indexes)

    save_ranges(start_indexes, end_indexes, out_path)

    return


def get_sizes(interactome_dict):
    num_interactions = pd.Series([interactome_dict[l].size() for l in LEVELS], index=LEVELS)
    num_entities = pd.Series([interactome_dict[l].graph['num_entities'] for l in LEVELS], index=LEVELS)
    num_small_molecules = pd.Series([interactome_dict[l].graph['num_small_molecules'] for l in LEVELS], index=LEVELS)

    return num_interactions, num_entities, num_small_molecules


def get_nodes_and_adjacent(nodes, G):
    """
    Get the same input nodes along with their adjacent small molecule nodes. Returns only those nodes that actually
    exist in G. The non-existing nodes are ignored.

    :param nodes: node set
    :param G: the complete graph
    :return: a set of nodes with their adjacent small molecule nodes
    """

    result = set()

    for node in nodes:
        if node in G.nodes:
            result.add(node)
            for neighbor in G.neighbors(node):
                if G.nodes[neighbor]['type'].startswith("S"):
                    result.add(neighbor)

    return result


def set_articulation_points(G):
    if not any("Articulation Point" in G.nodes[node] for node in G.nodes):
        nx.set_node_attributes(G, False, "Articulation Point")

        # Calculate articulation points
        art_points = list(nx.articulation_points(G))
        for node in art_points:
            G.nodes[node]["Articulation Point"] = True

        G.graph['Articulation Points'] = len(art_points)

def set_num_articulation_points(G):
    G.graph['Articulation Points'] = len([x for x,y in G.nodes(data=True) if y['Articulation Point']])

def set_bridges(G):
    if not any("Bridge" in G.edges[edge] for edge in G.edges):
        nx.set_edge_attributes(G, False, "Bridge")

        # Calculate bridges
        for edge in list(nx.bridges(G)):
            nx.set_edge_attributes(G, {edge: {"Bridge": True}})

def set_num_bridges(G):
    G.graph["Bridges"] = len([True for node1, node2, data in G.edges(data=True) if data['Bridge']])

def get_interactomes(graphs_path):
    participant_records = {l: get_participants(l, graphs_path) for l in [*LEVELS, sm]}
    components_records = {l: get_components(l, graphs_path) for l in [*LEVELS, sm]}

    interactomes_no_sm = {
        l: get_or_create_interaction_network(l, no_sm, participant_records, components_records, graphs_path) for l in
        LEVELS}
    interactomes_with_sm = {
        l: get_or_create_interaction_network(l, with_sm, participant_records, components_records, graphs_path) for l in
        LEVELS}
    interactomes_with_unique_sm = {
        l: get_or_create_interaction_network(l, with_unique_sm, participant_records, components_records, graphs_path)
        for l in LEVELS}
    interactomes = [*interactomes_no_sm.values(), *interactomes_with_sm.values(), *interactomes_with_unique_sm.values()]
    return interactomes

if __name__ == '__main__':
    print(f"Working directory: {os.getcwd()}")
    # graphs = create_pathway_interaction_networks("R-HSA-9673163", "../../../resources/pathway_networks/")

    df_pathways = get_pathways()
    for pathway in df_pathways["stId"]:
        create_pathway_interaction_networks(pathway, "../../../resources/pathway_networks/")
