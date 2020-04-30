# %%
import networkx as nx
from bokeh.io import output_file, show
from bokeh.models import Plot, Range1d, BoxZoomTool, ResetTool, Circle, MultiLine
from bokeh.models.graphs import from_networkx

from Python.interaction_network import create_graph
from Python.network_topology_queries import get_reactions_and_participants_by_pathway


def plot_graph(G, title=""):

    SMALL_MOLECULES = '#8DC7C5'
    PROTEINS = "teal"

    plot = Plot(plot_width=400, plot_height=400,
                x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1))
    plot.title.text = "Graph Interaction Demonstration"

    plot.add_tools(BoxZoomTool(), ResetTool())

    graph_renderer = from_networkx(G, nx.spring_layout, scale=1, center=(0, 0))

    spectral = ('#084594', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#deebf7', '#f7fbff')
    graph_renderer.node_renderer.glyph = Circle(size=15, fill_color=spectral[0])
    # graph_renderer.edge_renderer.glyph = MultiLine(line_color="edge_color", line_alpha=0.8, line_width=2)
    plot.renderers.append(graph_renderer)

    output_file("graph.html")
    show(plot)
    return plot

df = get_reactions_and_participants_by_pathway("R-HSA-70171")
G = create_graph(df)
plot_graph(G)
