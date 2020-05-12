# %%
import networkx as nx
from bokeh.io import output_file, save
from bokeh.layouts import layout
from bokeh.models import Plot, Range1d, BoxZoomTool, ResetTool, Circle, HoverTool, MultiLine, Div
from bokeh.models.graphs import from_networkx

from config import LEVELS, get_entity_color, COLOR_IO, COLOR_CO, COLOR_RO, COLOR_CC
from interaction_network import create_graph, get_pathway_graphs
from network_topology_queries import get_pathway_name, get_low_level_pathways


def plot_pathway(pathway, level="proteins", showSmallMolecules=True, verbose=True):
    G = create_graph(pathway, level, showSmallMolecules, verbose)
    P = plot_graph(G)
    return P


def plot_graph(G):
    """
    Plots the interaction network represented by Graph G
    :param G: networkx graph, it must have attributes "stId" with the Pathway id
    :return: the figure
    """

    plot = Plot(plot_width=600, plot_height=325,
                x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1),
                toolbar_location="below")
    plot.title.text = G.graph["level"].title()
    plot.title.text_font_size = '12pt'

    TOOLTIPS = """
            <div>
                <span style="font-size: 10px; font-weight: bold;">Name: @index</span>
            </div>
            <div>
                <span style="font-size: 10px; font-weight: bold; color: @color;">Type: @type</span>
            </div>
        """

    node_hover_tool = HoverTool(tooltips=TOOLTIPS)
    plot.add_tools(node_hover_tool, BoxZoomTool(), ResetTool())

    graph_renderer = from_networkx(G, nx.spring_layout, scale=1, center=(0, 0))

    graph_renderer.node_renderer.glyph = Circle(size=10, fill_color="color")
    graph_renderer.edge_renderer.glyph = MultiLine(line_color="edge_color", line_alpha=0.8, line_width=2)
    plot.renderers.append(graph_renderer)

    # output_file(f"{G.graph['level']}_network.html")
    # show(plot)
    return plot


def plot_pathway_all_levels(pathway, figures_path="../../figures/pathways/", graphs_path="../../reports/pathways/", verbose=False):
    name = get_pathway_name(pathway)
    if len(name) == 0:
        return

    graphs, graphs_no_small_molecules = get_pathway_graphs(pathway, graphs_path)

    plot_pathway_graphs(pathway, graphs, graphs_no_small_molecules, figures_path, verbose)


def plot_pathway_graphs(pathway, graphs, graphs_no_small_molecules, figures_path, verbose=False):
    name = get_pathway_name(pathway)
    if len(name) == 0:
        return
    name = name.iloc[0]['Name']

    figures = [plot_graph(graph) for graph in graphs]
    figures_no_small_molecules = [plot_graph(graph) for graph in graphs_no_small_molecules]

    if verbose:
        import os
        print("In function create_graph: ", os.getcwd())
        print(f"Output path is: {figures_path}{pathway}_network.html")
    output_file(f"{figures_path}{pathway}_network.html")

    title = f"<p style=\"font-weight:bold;text-align:left;font-size:22px;width:1800px;\">" \
            f"<span style=\"color:black;\">{name} ({pathway})</span>" \
            f"</p>"

    notes = f"""
    <style type="text/css">
    table {{
        border: 0px;
    }}

    .color_td {{
        border: inset 1px black;
        width: 50px
    }}

    table td#SimpleEntity {{
        background-color: {get_entity_color("SimpleEntity", "")};
    }}

    table td#EntityWithAccessionedSequence {{
        background-color: {get_entity_color("", "proteins")};
    }}

    table td#IO {{
        background-color: {COLOR_IO};
    }}

    table td#CO {{
        background-color: {COLOR_CO};
    }}

    table td#RO {{
        background-color: {COLOR_RO};
    }}

    table td#C {{
        background-color: {COLOR_CC};
    }}
</style>

<table>
    <tr>
        <td id="SimpleEntity" class="color_td"></td>
        <td>Small molecules</td>
    </tr>
    <tr>
        <td id="EntityWithAccessionedSequence" class="color_td"></td>
        <td>Genes, Proteins, Proteoforms</td>
    </tr>
    <tr>
        <td id="IO" class="color_td"></td>
        <td>Input -- Output</td>
    </tr>
    <tr>
        <td id="CO" class="color_td"></td>
        <td>Catalyst -- Output</td>
    </tr>
    <tr>
        <td id="RO" class="color_td"></td>
        <td>Regulator -- Output</td>
    </tr>
    <tr>
        <td id="C" class="color_td"></td>
        <td>Component -- Component</td>
    </tr>
</table>
    """

    l = layout([
        [Div(text=f"{title}")],
        [Div(text=notes)],
        figures,
        [Div(text=f"<p>Without small molecules:</p>")],
        figures_no_small_molecules
    ])

    save(l)


def plot_low_level_pathways(min_size=5, max_size=20, figures_path="figures/pathways/", graphs_path="reports/pathways/"):
    pathways = get_low_level_pathways()
    for pathway in pathways['stId']:
        name = get_pathway_name(pathway)
        name = name.iloc[0]['Name']
        print(f"Creating networks for pathway {pathway}")
        G = create_graph(pathway, "genes", True, graphs_path)

        if min_size <= len(G.nodes) <= max_size:
            plot_pathway_all_levels(pathway, figures_path, graphs_path)

    print(f"Plot {pathway}")
    plot_pathway_all_levels(pathway, figures_path, graphs_path)


plot_pathway_all_levels("R-HSA-165160")

# plot_low_level_pathways(figures_path="../../figures/pathways/", graphs_path="../../reports/pathways/")
# print("Finished")
