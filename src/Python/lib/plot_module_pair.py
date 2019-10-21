# %%
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

# %%
# Creates an image of a visual representation of two modules
# trait1, trait2: string values of the desired trait modules
# level: gene, protein or proteoform
def plot_module_pair_old(trait1, trait2, level, path_to_root="../../../"):
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