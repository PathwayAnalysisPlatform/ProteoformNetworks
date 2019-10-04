# %%
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

LEVELS = {"gene", "protein", "proteoform"}


# %%
# Creates an image of a visual representation of two modules
# trait1, trait2: string values of the desired trait modules
# level: gene, protein or proteoform
def plot_module_pair(trait1, trait2, level, path_to_root="../../../"):
    # Build your graph
    G = nx.Graph()
    print("Creating graph 1")
    g1 = get_graph(trait1, level, path_to_root=path_to_root)
    # g1 = nx.Graph()
    # g1.add_edges_from([(1, 2), (1, 3)])
    print("Creating graph 2")
    g2 = get_graph(trait2, level, path_to_root=path_to_root)
    # g2 = nx.Graph()
    # g2.add_edges_from([(2, 4), (3, 4)])

    print("Merging graphs")
    G.add_edges_from(g1.edges)
    G.add_edges_from(g2.edges)

    # And a data frame with characteristics for your nodes
    locations = []
    for node in G.nodes:
        if node in g1.nodes and node in g2.nodes:
            locations.append(3)
        elif node in g1.nodes:
            locations.append(1)
        else:
            locations.append(2)
    carac = pd.DataFrame(
        {'ID': list(G.nodes), 'location': locations})

    print("Number of nodes: ", G.number_of_nodes())
    print("Number of edges: ", G.number_of_edges())

    # The order of the node for networkX is the following order:
    print(G.nodes())
    # Thus, we cannot give directly the 'myvalue' column to netowrkX, we need to arrange the order!

    # Here is the tricky part: I need to reorder carac to assign the good color to each node
    carac = carac.set_index('ID')
    carac = carac.reindex(G.nodes())

    # And I need to transform my categorical column in a numerical value: group1->1, group2->2...
    carac['location'] = pd.Categorical(carac['location'])
    carac['location'].cat.codes

    print(f"Nodes in {trait1}:");
    print(g1.nodes)

    print(f"Nodes in {trait2}:")
    print(g2.nodes)

    print("Locations: ")
    print(locations)

    print("Plotting")
    options = {
        'node_size': 500,
        'width': 3,
        'font_size': 8,
        'with_labels': True,
        'node_color': carac['location'].cat.codes,
        'cmap': plt.cm.Set1
    }
    plt.subplot(221)
    nx.draw_spring(G, **options)

    options = {
        'node_size': 500,
        'font_size': 6,
        'width': 3,
        'with_labels': True,
        'node_color': carac['location'].cat.codes,
        'cmap': plt.cm.Set1
    }
    plt.subplot(222)
    nx.draw_random(G, **options)

    file_name = f"{path_to_root}figures/module_pairs/{trait1}_-{trait2}_{level}.png"
    file_name = file_name.replace("\"", "").replace(" ", "-")
    plt.savefig(file_name)
    plt.show()


import os

print("WD: ", os.getcwd())
plot_module_pair("\"Anemia, Sickle Cell\"", "Bilirubin", "gene", path_to_root="../../")
plot_module_pair("\"Anemia, Sickle Cell\"", "Bilirubin", "protein", path_to_root="../../")
plot_module_pair("\"Anemia, Sickle Cell\"", "Bilirubin", "proteoform", path_to_root="../../")


# %%
def get_graph(trait, level, path_to_root="../../../"):
    if level not in LEVELS:
        raise ValueError("level must be one of %r." % LEVELS)

    G = nx.Graph()

    # Traverse file with module members. Get set of members for the trait
    print("Reading members of module ", trait)
    path_file_modules = path_to_root + "reports/modules/" + level + '_modules.tsv'
    with open(path_file_modules) as file_modules:
        file_modules.readline()  # Read header
        line = file_modules.readline()
        while line:
            columns = line.split("\t")
            if columns[0] == trait:
                G.add_node(columns[1].strip())
            line = file_modules.readline()

    # Traverse file with interactions of the level. Get the set of edges for this trait
    print("Reading edges in module ", trait)
    path_file_interactions = path_to_root + "resources/Reactome/v70/"
    if level == "gene":
        path_file_interactions += "Genes/geneInternalEdges.tsv"
    elif level == "protein":
        path_file_interactions += "Proteins/proteinEdges.tsv"
    else:
        path_file_interactions += "Proteoforms/proteoformEdges.tsv"

    nodes = set(G.nodes)
    with open(path_file_interactions) as file_interactions:
        file_interactions.readline()  # Read header
        line = file_interactions.readline()
        while line:
            columns = line.split("\t")
            if columns[0] in nodes and columns[1] in nodes:
                G.add_edge(columns[0], columns[1])
            line = file_interactions.readline()

    return G
