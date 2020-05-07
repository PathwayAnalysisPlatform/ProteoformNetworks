# %%
import os
import networkx as nx
from bokeh.io import show, output_file
from bokeh.models import Plot, Range1d, BoxZoomTool, ResetTool, Circle, HoverTool, MultiLine
from bokeh.models.graphs import from_networkx
from interaction_network import create_graph


def plot_pathway(pathway, level="proteins", directed=False, showSmallMolecules=True, verbose=True):
    G = create_graph(pathway, level, directed, showSmallMolecules, verbose)
    P = plot_graph(G)
    return P


def plot_graph(G):
    """
    Plots the interaction network represented by Graph G
    :param G: networkx graph, it must have attributes "stId" with the Pathway id
    :return: the figure
    """

    plot = Plot(plot_width=250, plot_height=250,
                x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1),
                toolbar_location="below")
    plot.title.text = G.graph["level"].title()
    plot.title.text_font_size = '12pt'

    TOOLTIPS = """
            <div>
                <span style="font-size: 18px; font-weight: bold;">Name: @index</span>
            </div>
            <div>
                <span style="font-size: 18px; font-weight: bold; color: @color;">Type: @type</span>
            </div>
        """

    node_hover_tool = HoverTool(tooltips=TOOLTIPS)
    plot.add_tools(node_hover_tool, BoxZoomTool(), ResetTool())

    graph_renderer = from_networkx(G, nx.spring_layout, scale=1, center=(0, 0))

    graph_renderer.node_renderer.glyph = Circle(size=10, fill_color="color")
    graph_renderer.edge_renderer.glyph = MultiLine(line_color="edge_color", line_alpha=0.8, line_width=2)
    plot.renderers.append(graph_renderer)

    output_file("network.html")
    show(plot)
    return plot


def plot_pathway_all_levels():
    pass


# print(f"Located at: {os. getcwd()}")
# plot_pathway("R-HSA-9634600", level="proteins")
