# %%
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

LEVELS = {"gene", "protein", "proteoform"}


def get_graph(trait, level, path_to_root="../../../"):
    if level not in LEVELS:
        raise ValueError("level must be one of %r." % LEVELS)

    G = nx.Graph()

    # Traverse file with module members. Get set of members for the trait
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


# Creates an image of a visual representation of two modules
# trait1, trait2: string values of the desired trait modules
# level: gene, protein or proteoform
def plot_module_pair(trait1, trait2, level):
    # Build your graph
    G = nx.Graph()
    trait1_graph = get_graph(trait1, level)
    trait2_graph = get_graph(trait2, level)
    G.add_edges_from(trait1_graph)
    G.add_edges_from(trait2_graph)

    # And a data frame with characteristics for your nodes
    carac = pd.DataFrame(
        {'ID': ['A', 'B', 'C', 'D', 'E'], 'location': ['group1', 'group1', 'group2', 'overlap', 'overlap']})

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

    # Custom the nodes:
    plt.subplot(2, 1, 1)
    nx.draw(G, with_labels=True, node_color=carac['location'].cat.codes, cmap=plt.cm.Set1, node_size=1500)
    plt.savefig("../../figures/module_pairs/test.png")
    plt.show()

# plot_module_pair("\"Anemia, Sickle Cell\"", "Bilirubin", "gene")
