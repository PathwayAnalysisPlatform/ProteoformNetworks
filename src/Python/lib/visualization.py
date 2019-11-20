# %%
import matplotlib.pyplot as plt
import networkx as nx
from bokeh.io import output_file, show
from bokeh.layouts import layout
from bokeh.models import Div
from bokeh.models import Plot, Range1d, MultiLine, Circle, HoverTool, BoxZoomTool, ResetTool, PointDrawTool
from bokeh.models.graphs import from_networkx

from config import levels
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


def plot_module_pair(trait1, trait2, level, path_to_modules, path_to_figures):
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
    colors = {}
    for node in G.nodes:
        if node in g1.nodes and node in g2.nodes:
            region[node] = "OVERLAP"
            colors[node] = "red"
        elif node in g1.nodes:
            region[node] = trait1
            colors[node] = "green"
        else:
            region[node] = trait2
            colors[node] = "yellow"
    nx.set_node_attributes(G, colors, "colors")
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
    plot = Plot(plot_width=600, plot_height=600,
                x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1))
    plot.title.text = "Modules for \"" + trait1 + "\" and \"" + trait2 + "\""
    plot.title.text_font_size = '24pt'

    TOOLTIPS = """
        <div>
            <span style="font-size: 18px; font-weight: bold;">Name: @index</span>
        </div>
        <div>
            <span style="font-size: 18px; font-weight: bold; color: @colors;">Region: @region</span>
        </div>
    """
    node_hover_tool = HoverTool(tooltips=TOOLTIPS)
    plot.add_tools(node_hover_tool, BoxZoomTool(), ResetTool())

    graph_renderer = from_networkx(G, nx.spring_layout, scale=1, center=(0, 0))

    graph_renderer.node_renderer.glyph = Circle(size=20, fill_color="colors")
    graph_renderer.edge_renderer.glyph = MultiLine(line_color="edge_color", line_alpha=0.8, line_width=2)
    plot.renderers.append(graph_renderer)

    output_file(f"{path_to_figures}{trait1}_to_{trait2}_{level}.html")
    show(plot)
    return plot


def create_graph(trait1, trait2, level, path_to_modules, only_interface=False):
    graph_complete = nx.Graph()
    graph_interface = nx.Graph()

    traits = [trait1, trait2]
    graphs = {trait: get_graph(trait, level, path_to_modules) for trait in traits}
    [graph_complete.add_edges_from(graph.edges) for graph in graphs.values()]  # Merge both graphs

    CROSSING_EDGES, IN_MODULE_1_EDGES, IN_MODULE_2_EDGES, IN_OVERLAP_EDGES = "orange", "darkgreen", "darkblue", "red"

    region = {}
    colors = {}
    for node in graph_complete.nodes:
        if node in graphs[trait1].nodes and node in graphs[trait2].nodes:
            region[node] = "OVERLAP"
            colors[node] = "red"
        elif node in graphs[trait1].nodes:
            region[node] = trait1
            colors[node] = "green"
        else:
            region[node] = trait2
            colors[node] = "blue"
    nx.set_node_attributes(graph_complete, colors, "colors")
    nx.set_node_attributes(graph_complete, region, "region")

    edge_attrs = {}
    interface_nodes = set()
    for start_node, end_node, _ in graph_complete.edges(data=True):
        if graph_complete.nodes[start_node]["region"] == graph_complete.nodes[end_node]["region"]:
            if graph_complete.nodes[start_node]["region"] == "OVERLAP":
                edge_color = IN_OVERLAP_EDGES
                interface_nodes.add(start_node)
                interface_nodes.add(end_node)
            elif graph_complete.nodes[start_node]["region"] == trait1:
                edge_color = IN_MODULE_1_EDGES
            else:
                edge_color = IN_MODULE_2_EDGES
        else:
            edge_color = CROSSING_EDGES
            interface_nodes.add(start_node)
            interface_nodes.add(end_node)
        edge_attrs[(start_node, end_node)] = edge_color
    nx.set_edge_attributes(graph_complete, edge_attrs, "edge_color")

    if only_interface:
        graph_interface = graph_complete.subgraph(interface_nodes)

    print("All nodes:")
    for node in graph_complete.nodes:
        print(node)

    print("Interface nodes: ")
    for node in graph_interface.nodes:
        print(node)

    return graph_interface if only_interface else graph_complete


def create_plot(level, graph):
    plot = Plot(plot_width=600, plot_height=450,
                x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1),
                toolbar_location="below")
    plot.title.text = f"{level.title()}"
    plot.title.align = 'center'
    plot.title.text_font_size = '16pt'

    TOOLTIPS = """
            <div>
                <span style="font-size: 18px; font-weight: bold;">Name: @index</span>
            </div>
            <div>
                <span style="font-size: 18px; font-weight: bold; color: @colors;">Region: @region</span>
            </div>
        """
    node_hover_tool = HoverTool(tooltips=TOOLTIPS)
    plot.add_tools(node_hover_tool, BoxZoomTool(), ResetTool(), PointDrawTool())

    graph_renderer = from_networkx(graph, nx.spring_layout, scale=1, center=(0, 0))

    if len(graph.nodes) > 0:
        graph_renderer.node_renderer.glyph = Circle(size=6, fill_color="colors")
        graph_renderer.edge_renderer.glyph = MultiLine(line_color="edge_color", line_alpha=0.8, line_width=1.5)
    plot.renderers.append(graph_renderer)
    return plot


def plot_module_pair(trait1, trait2, path_to_modules, path_to_figures):
    """Visualize the pair of modules in an interactive graph.

    trait1, trait2: names of the two sets as strings
    path_to_modules: directory where to find the input files
    path_to_figures: directory where to create the output files
    """

    graphs_complete = {level: create_graph(trait1, trait2, level, path_to_modules) for level in levels}
    graphs_interface = {level: create_graph(trait1, trait2, level, path_to_modules, only_interface=True)
                        for level in levels}

    figures_complete_modules = [create_plot(level, graph) for level, graph in graphs_complete.items()]
    figures_interfaces = [create_plot(level, graph) for level, graph in graphs_interface.items()]

    output_file(f"{path_to_figures}{trait1}_to_{trait2}.html")
    title = f"<p style=\"font-weight:bold;text-align:center;font-size:22px;width:1800px;\">" \
            f"<span style=\"color:green;\">{trait1}</span> with " \
            f"<span style=\"color:blue;\">{trait2}</span>" \
            f"</p>"

    show(layout([[Div(text=f"{title}")], figures_complete_modules, figures_interfaces]))
