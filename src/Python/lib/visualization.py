# %%
import networkx as nx
import matplotlib.pyplot as plt

from bokeh.io import output_file, show
from bokeh.models import Plot, Range1d, MultiLine, Circle, HoverTool, BoxZoomTool, ResetTool
from bokeh.models.graphs import from_networkx


# Plot the frequency of each overlap score
# Only considers positive scores, where sets overlap.
from lib.data_read_write import get_graph


def plot(x, score_label='overlap'):
    # Consider only overlapping sets
    x = x[x != 0]
    prob = x.value_counts()
    plt.hist(x, bins=40)
    plt.title(f"Histogram of {score_label} score")
    plt.xlabel('Score')
    plt.ylabel('Frequency')
    plt.show()


# Create the label from the file name.
# All scoring files must follow the naming convention: scores_label.tsv, where 'label' is the name
# of the overlap scoring method
def getScoreLabel(file_name):
    index = file_name.rfind('/')
    score_label = file_name[(0 if index == -1 else index + 1):]
    score_label = score_label.replace("scores_", "").replace("_", " ").replace(".tsv", "")
    return score_label


def plot_module_pair(trait1, trait2, level, path_to_modules):
    """Visualize the pair of modules in an interactive graph."""
    # Prepare the graph
    G = nx.Graph()
    print("Creating graph 1")
    g1 = get_graph(trait1, level, path_to_modules)
    # g1 = nx.gnp_random_graph(5, 0, 5)
    print("Creating graph 2")
    g2 = get_graph(trait2, level, path_to_modules)
    # g2 = nx.gnp_random_graph(7, 0.5)

    G.add_edges_from(g1.edges)
    G.add_edges_from(g2.edges)

    region = {}
    for node in G.nodes:
        if node in g1.nodes and node in g2.nodes:
            region[node] = "OVERLAP"
        elif node in g1.nodes:
            region[node] = "1"
        else:
            region[node] = "2"
    nx.set_node_attributes(G, region, "region")
    print(region)

    # %%
    CROSSING_EDGES, IN_MODULE_EDGES, IN_OVERLAP_EDGES = "red", "black", "orange"

    edge_attrs = {}
    for start_node, end_node, _ in G.edges(data=True):
        if G.nodes[start_node]["region"] == G.nodes[end_node]["region"]:
            if G.nodes[start_node]["region"] == "OVERLAP":
                edge_color = IN_OVERLAP_EDGES
            else:
                edge_color = IN_MODULE_EDGES
        else:
            edge_color = CROSSING_EDGES
        edge_attrs[(start_node, end_node)] = edge_color
    print(edge_attrs)
    nx.set_edge_attributes(G, edge_attrs, "edge_color")
    # %%
    # Show with Bokeh
    plot = Plot(plot_width=1500, plot_height=1000,
                x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1))
    plot.title.text = "Modules for " + trait1 + " and " + trait2

    node_hover_tool = HoverTool(tooltips=[("name", "@index"), ("region", "@region")])
    plot.add_tools(node_hover_tool, BoxZoomTool(), ResetTool())

    graph_renderer = from_networkx(G, nx.spring_layout, scale=1, center=(0, 0))

    graph_renderer.node_renderer.glyph = Circle(size=15, fill_color="blue")
    graph_renderer.edge_renderer.glyph = MultiLine(line_color="edge_color", line_alpha=0.8, line_width=1)
    plot.renderers.append(graph_renderer)

    output_file("interactive_graphs.html")
    show(plot)
