# %%
import os
import sys
from enum import Enum
from math import pi

import networkx as nx
import pandas as pd
from bokeh.io import output_file, save
from bokeh.io import show
from bokeh.layouts import layout
from bokeh.models import (ColumnDataSource, Legend, BoxSelectTool, PointDrawTool)
from bokeh.models import Div
from bokeh.palettes import Colorblind
# Prepare Data
from bokeh.plotting import figure
from bokeh.transform import cumsum

import config
from config import get_entity_color, COLOR_IO, COLOR_CO, COLOR_RO, COLOR_CC, no_sm
from lib.graph_database_access import get_pathway_name
from lib.networks import create_pathway_interaction_networks


class Coloring(Enum):
    ENTITY_TYPE = "Entity Type"
    REACTION = "Reaction"
    PATHWAY = "Pathway"
    BRIDGES = "Bridges"


def show_graph_with_lcc(g):
    """
    Plot a graph and highlight the nodes and edges forming the largest connected component
    :return: plot
    """
    pass


# Plot graph with multiple types of nodes

# Plot graph with multiple types of nodes and highlight the largest connected component


def getDefaultPalette(n):
    if n == 1:
        return [Colorblind[3][0]]
    elif n == 2:
        return [Colorblind[3][0], Colorblind[3][1]]
    else:
        return Colorblind[n]


def plot_interaction_network(G, coloring=Coloring.ENTITY_TYPE, **kwargs):
    """
    Plots the interaction network represented by Graph G
    :param G: networkx graph, it must have attributes "stId" with the Pathway id
    :param coloring: Type of coloring for the nodes and links, by entity type, reaction or pathway
    :return: the figure
    """

    plot_width = kwargs['plot_width'] if 'plot_width' in kwargs else 600
    plot_height = kwargs['plot_height'] if 'plot_height' in kwargs else 325
    toolbar_location = kwargs['toolbar_location'] if 'toolbar_location' in kwargs else None

    TOOLTIPS = [
        ("Id", "@id"),
        ("Type", "@type")
    ]



    f = figure(x_range=(-1.1, 1.1), y_range=(-1.1, 1.1), toolbar_location=toolbar_location, tooltips=TOOLTIPS)
    # node_hover_tool = HoverTool(tooltips=TOOLTIPS, names=['node'])

    f.title.text = kwargs['title'] if 'title' in kwargs else G.graph['level'].title()
    f.title.text_font_size = '12pt'
    f.title.align = 'center'
    f.axis.visible = False
    f.axis.axis_label = None
    f.grid.grid_line_color = None
    f.outline_line_width = kwargs['outline_line_width'] if 'outline_line_width' in kwargs else 0

    pos = kwargs['pos'] if 'pos' in kwargs else nx.spring_layout(G)
    if kwargs["v"] if "v" in kwargs else False:
        print(pos)

    # Display edges

    highlight_bridges = 'highlight_bridges' in kwargs and kwargs['highlight_bridges']

    for start, end in G.edges:
        if kwargs["v"] if "v" in kwargs else False:
            print(f"Edge from {start} -- {end}")
        x_start, y_start = pos[start]
        x_end, y_end = pos[end]
        line_color = "#FF0000" if highlight_bridges and G.edges[start, end]['Bridge'] else "black"
        line_width = 4 if highlight_bridges and G.edges[start, end]['Bridge'] else 2
        f.line([x_start, x_end], [y_start, y_end], line_color=line_color, line_width=line_width)

    if coloring == Coloring.ENTITY_TYPE:

        ids = []
        types = []
        x_values = []
        y_values = []
        entity_colors = []
        entity_border_colors = []
        entity_border_widths = []

        ids_sm = []
        types_sm = []
        x_values_sm = []
        y_values_sm = []
        entity_colors_sm = []
        entity_border_colors_sm = []
        entity_border_widths_sm = []

        for id in G.nodes:
            if G.nodes[id]['type'] == "SimpleEntity":
                if id not in pos.keys():
                    raise RuntimeError(f"Position not found for Simple Entity with key {id}")
                x_values_sm.append(pos[id][0])
                y_values_sm.append(pos[id][1])
                entity_colors_sm.append(G.nodes[id]['entity_color'])
                types_sm.append(G.nodes[id]['type'])
                ids_sm.append(G.nodes[id]['label'])
                if 'highlight_articulations' in kwargs and kwargs['highlight_articulations'] and G.nodes[id][
                    'Articulation Point']:
                    entity_border_colors_sm.append("#FF0000")
                    entity_border_widths_sm.append(3)
                else:
                    entity_border_colors_sm.append("black")
                    entity_border_widths_sm.append(1)


            else:
                if id not in pos.keys():
                    raise RuntimeError(f"Position not found for Accessioned Entity with key {id}")
                x_values.append(pos[id][0])
                y_values.append(pos[id][1])
                entity_colors.append(G.nodes[id]['entity_color'])
                types.append(G.nodes[id]['type'])
                ids.append(G.nodes[id]['label'])
                if 'highlight_articulations' in kwargs and kwargs['highlight_articulations'] and G.nodes[id][
                    'Articulation Point']:
                    entity_border_colors.append("#FF0000")
                    entity_border_widths.append(3)
                else:
                    entity_border_colors.append("black")
                    entity_border_widths.append(1)

        data = {
            'id': ids,
            'type': types,
            'x': x_values,
            'y': y_values,
            'entity_color': entity_colors,
            'entity_border_color': entity_border_colors,
            'line_width': entity_border_widths
        }
        data_sm = {
            'id': ids_sm,
            'type': types_sm,
            'x': x_values_sm,
            'y': y_values_sm,
            'entity_color': entity_colors_sm,
            'entity_border_color': entity_border_colors_sm,
            'line_width': entity_border_widths_sm

        }

        if len(data_sm['id']) > 0:
            f.circle(x='x', y='y',
                     size=20,
                     fill_color='entity_color', line_color='entity_border_color', line_width='line_width',
                     legend_label="Small Molecules",
                     source=ColumnDataSource(data_sm))

        f.square(x='x', y='y',
               size=20,
               fill_color='entity_color', line_color='entity_border_color', line_width='line_width',
               legend_label=G.graph["level"].title(),
               source=ColumnDataSource(data))
        f.legend.title = coloring.value

    else:
        if coloring == Coloring.REACTION:
            event = 'reaction'
        else:
            event = 'pathway'
        # Display nodes
        events = set()
        for k in G.nodes:
            for e in G.nodes[k][event + 's']:
                events.add(e)

        radius = kwargs['radius'] if 'radius' in kwargs else 0.1
        color_palette = kwargs['color_palette'] if 'color_palette' in kwargs else getDefaultPalette(len(events))

        ######## Create fake wedge to create full legend
        data = pd.Series(
            {e: 1 for e in events}).reset_index(name='value').rename(columns={'index': event})
        data['id'] = "-"
        data['type'] = "SimpleEntity"
        data['entity_color'] = get_entity_color("SimpleEntity", "genes")
        data['angle'] = data['value'] / data['value'].sum() * 2 * pi
        data[event + '_color'] = color_palette

        f.wedge(x=0, y=0, radius=0,
                start_angle=cumsum('angle', include_zero=True),
                end_angle=cumsum('angle'), legend_field=event, source=data, visible=False)

        ######## Add all nodes

        for k in G.nodes:
            x, y = pos[k]
            if kwargs["v"] if "v" in kwargs else False:
                print(f"Node {k} has position {x}, {y}")
                print(f"Belongs to {event}s {G.nodes[k][event + 's']}")
            data = pd.Series(
                {e: 1 if e in G.nodes[k][f"{event}s"] else 0 for e in events}).reset_index(
                name='value').rename(
                columns={'index': event})
            data['id'] = k
            data['type'] = G.nodes[k]['type']
            data['entity_color'] = G.nodes[k]['entity_color']
            data['angle'] = data['value'] / data['value'].sum() * 2 * pi
            data[event + '_color'] = color_palette

            i = 0
            count = 0
            index = 0
            for e in data['value']:
                if e == 1:
                    count += 1
                    index = i
                i += 1

            if count != 1:
                f.wedge(x=x, y=y, radius=radius,
                        start_angle=cumsum('angle', include_zero=True),
                        end_angle=cumsum('angle'),
                        line_color="black", fill_color=event + '_color', legend_field=event, source=data,
                        name='node')
            else:
                data = {
                    'id': [k], 'type': [G.nodes[k]['type']],
                    event + '_color': [data[event + '_color'][index]],
                    'entity_color': [G.nodes[k]['entity_color']]
                }
                f.circle(x=x, y=y, radius=radius,
                         line_color='black', fill_color=event + '_color',
                         source=ColumnDataSource(data=data), name='node')

        ###### Relocate the legend
        if 'legend_location' in kwargs:
            if kwargs['legend_location'] == 'right':
                legend = Legend(items=f.legend.items, location='top_right', title=coloring.value)
                f.legend.visible = True
                f.add_layout(legend, 'right')
            elif kwargs['legend_location']:
                f.legend.location = kwargs['legend_location']
            else:
                f.legend.visible = False

    output_file(f"{G.graph['level']}_network.html")
    # show(f)
    return f


def select_common_nodes(smaller_graph, larger_graph):
    """
    Create a dictionary with the mapping from one entity type to the other.
    Protein graph nodes, should know which gene they belong to.
    Proteoform graph nodes, should know which protein they belong to.

    :return:
    """
    common = {}
    if larger_graph.graph['level'] == 'proteins':
        for node in larger_graph.Nodes():
            if larger_graph.Nodes[node]['type'] == 'SimpleEntity':
                common[node] = node
            else:
                common[node] = larger_graph.Nodes[node]['gene']
    return common


def get_positions(graphs):
    """
    Get positions for the vertices of the 9 graphs, methods x levels combinations

    :param graphs: dictionary with methods to construct graphs as methods and a list of graphs (for each level) as values
    :return: dictionary with the same keys as graphs but with positions of vertices as values
    """
    # The graph that contains all other nodes and more is the proteoforms with unique small molecule ids in each reaction

    # Vertically:
    # -- With small molecules: Keep proteoforms in the same position and small molecules are recalculated
    # -- Without small molecules: Keep proteoforms in the same position.

    pos = {
        config.no_sm: [{} for g in graphs[config.no_sm]],
        config.with_sm: [{} for g in graphs[config.with_sm]],
        config.with_unique_sm: [{} for g in graphs[config.with_unique_sm]]
    }

    # For the proteoforms network With reaction-unique small molecules
    pos[config.with_unique_sm][2] = nx.spring_layout(graphs[config.with_unique_sm][2])

    # For the proteoforms network With not unique small molecules
    # Fix all proteoforms and set the position for the small molecules freely
    fixed_nodes = []
    fixed_positions = {}
    for node in graphs[config.with_unique_sm][2].nodes:  # For each node in the reaction-unique id
        if node.startswith("sm"):
            name_without_reaction = graphs[config.with_unique_sm][2].nodes[node]['label']
            fixed_positions[name_without_reaction] = pos[config.with_unique_sm][2][node]  # Get the name of the small molecule without the reaction
        else:
            fixed_nodes.append(node)
            fixed_positions[node] = pos[config.with_unique_sm][2][node]
    pos[config.with_sm][2] = nx.spring_layout(graphs[config.with_sm][2], pos=fixed_positions, fixed=fixed_positions.keys())

    # Horizontally:
    # With unique small molecules:
    # -- Proteins take the position of the first proteoform with the same accession. Small molecules are fixed.
    # -- Genes take the position of the first protein comming from the same gene. Small molecules are fixed.

    # For the proteoforms network With reaction-unique small molecules
    # 0 - genes, 1 - proteins, 2 proteoforms
    for i in reversed(range(2)):  # Select graphs of proteins and genes with indexes 1 and 0
        for node in graphs[config.with_unique_sm][i + 1].nodes:                              # For each node in the larger network
            prev_id = graphs[config.with_unique_sm][i + 1].nodes[node]['prevId']              # Get the predecesor
            pos[config.with_unique_sm][i][prev_id] = pos[config.with_unique_sm][i + 1][node]  # Set position in the smaller network
            print(f"Set position for {prev_id} using {node}")

    # For the proteoforms network With not unique small molecules
    # Leave small molecules fixed, recalculate the others
    # 0 - genes, 1 - proteins, 2 proteoforms
    for i in reversed(range(2)):  # Select graphs of proteins and genes with indexes 1 and 0
        fixed_positions = {}
        for node in graphs[config.with_sm][i + 1].nodes:                                    # For each node in the larger network
            if node.startswith("sm"):
                fixed_positions[node] = pos[config.with_sm][i+1][node]
        pos[config.with_sm][i] = nx.spring_layout(graphs[config.with_sm][i], pos=fixed_positions, fixed=fixed_positions.keys())

    # For the proteoforms network without small molecules
    # Copy directly the position of all proteoforms
    pos[config.no_sm] = pos[config.with_sm]

    return pos


def plot_pathway_all_levels(pathway, out_path="../../figures/pathways/", graphs_path="../../reports/pathways/",
                            coloring=Coloring.ENTITY_TYPE,
                            graphs=None, **kwargs):
    """
    Plot the interaction networks generated for a single pathway at the 3 levels; genes, proteins and proteoforms and
    using the 3 methods of constructing the interaction network.
    It locates the entities in the same coordinates as their translation products. Ex. Gene, protein and proteoform in
    the same relative location in its own plot.

    :param pathway: string id for the pathway
    :param out_path:
    :param graphs_path:
    :param coloring:
    :param graphs: [list of graphs one for each level]}
    :param v:
    :return:
    """
    name = get_pathway_name(pathway)['Name'][0]
    if len(name) == 0:
        return

    graphs = {
        config.no_sm: list(graphs[no_sm].values()),
        config.with_sm: list(graphs[config.with_sm].values()),
        config.with_unique_sm: list(graphs[config.with_unique_sm].values())
    }

    titles = {
        config.no_sm: ['A) Genes without SM', 'B) Proteins without SM', 'C) Proteoforms without SM'],
        config.with_sm: ['D) Proteins with SM', 'E) Proteins with SM', 'F) Proteoforms with SM'],
        config.with_unique_sm: ['G) Proteins with unique SM', 'H) Proteins with unique SM',
                                'I) Proteoforms with unique SM']
    }

    plot_size = 600
    legend_location_all = [None, None, 'right']
    plot_widths = [plot_size, plot_size, plot_size + 160]
    if coloring == Coloring.ENTITY_TYPE:
        legend_location_all = ['top_right', 'top_right', 'top_right']
        plot_widths = [plot_size, plot_size, plot_size]

    toolbar_location = kwargs['toolbar_location'] if 'toolbar_location' in kwargs else None

    pos = get_positions(graphs)

    figures = {
        method: [
            plot_interaction_network(graphs[method][i], coloring=coloring, pos=pos[method][i],
                                     plot_width=plot_widths[i],
                                     plot_height=plot_size,
                                     toolbar_location=toolbar_location, title=titles[method][i],
                                     legend_location=legend_location_all[i],
                                     highlight_articulations=(kwargs[
                                                                  'highlight_articulations'] if 'highlight_articulations' in kwargs else False),
                                     highlight_bridges=(
                                         kwargs['highlight_bridges'] if 'highlight_bridges' in kwargs else False))
            for i in range(3)
        ]
        for method in config.METHODS
    }

    if kwargs["v"] if "v" in kwargs else False:
        import os
        print("In function create_graph: ", os.getcwd())
        print(f"Output path is: {out_path}{pathway}_network.html")
    output_file(f"{out_path}{pathway}_{coloring.name}_network.html")

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
    """
    # <tr>
    #     <td id="IO" class="color_td"></td>
    #     <td>Input -- Output</td>
    # </tr>
    # <tr>
    #     <td id="CO" class="color_td"></td>
    #     <td>Catalyst -- Output</td>
    # </tr>
    # <tr>
    #     <td id="RO" class="color_td"></td>
    #     <td>Regulator -- Output</td>
    # </tr>
    # <tr>
    #     <td id="C" class="color_td"></td>
    #     <td>Component -- Component</td>
    # </tr>
    notes += """
    </table>
    """

    l = layout([
        [Div(text=f"{title}")],
        # [Div(text=notes)],
        figures[config.no_sm],
        figures[config.with_sm],
        figures[config.with_unique_sm]
    ])

    # show(l)
    save(l)
    return l


#
# def plot_low_level_pathways(min_size=5, max_size=20, figures_path="figures/pathways/", graphs_path="reports/pathways/"):
#     pathways = get_low_level_pathways()
#     for pathway in pathways['stId']:
#         name = get_pathway_name(pathway)
#         name = name.iloc[0]['Name']
#
#         G = create_graph(pathway, "genes", True, graphs_path, save=False)
#
#         if min_size <= len(G.nodes) <= max_size:
#             print(f"-- Plotting pathway \"{pathway}\" with size {len(G.nodes)}")
#             plot_pathway_all_levels(pathway, figures_path, graphs_path)
#         # else:
#         #     print(f"Skipping pathway: \"{pathway}\" with size {len(G.nodes)}")

def plot_pathways(pathways, level, sm, coloring, v=False):
    """
    Creates an interaction network composed by multiple pathways and plots the nodes according to the pathway
    each node belongs to

    :param pathways: Collection of pathway stId strings
    :param level: {trouble}
    :param sm:
    :param coloring: Coloring scheme for the nodes and edges: by entity type, reactions and pathways
    :param v:
    :return: bokeh figure with the plot
    """
    # graphs = [create_graph(pathway, level, sm) for pathway in pathways]
    # full_graph = merge_graphs(graphs)
    # f = plot_interaction_network(full_graph, Coloring.PATHWAY)
    # show(f)
    # return f
    pass


def main():
    pathway1 = "R-HSA-1236382"  # Constitutive Signaling by Ligand-Responsive EGFR Cancer Variants
    pathway2 = "R-HSA-6814122"  # Cooperation of PDCL (PhLP1) and TRiC/CCT in G-protein beta folding
    pathway3 = "R-HSA-9648002"  # Ras Processing

    os.chdir(os.path.dirname(os.path.abspath(sys.executable)) + "\\..\\..")
    print(os.getcwd())

    # g = create_pathway_interaction_network(pathway1, proteoforms, with_unique_sm, "resources/pathway_networks/" + pathway1)
    # p = plot_interaction_network(g, coloring=Coloring.ENTITY_TYPE, plot_width=600, plot_height=500,
    #                             title="Test title",
    #                             legend_location="right")
    graphs = create_pathway_interaction_networks(pathway1, "resources/pathway_networks/")
    p = plot_pathway_all_levels(pathway1, out_path="resources/pathway_networks/", graphs=graphs,
                                coloring=Coloring.ENTITY_TYPE,
                                highlight_articulations=True,
                                highlight_bridges=True,
                                toolbar_location="below")
    show(p)
    # p = plot_pathway_all_levels(pathway, graphs=graphs, coloring=Coloring.ENTITY_TYPE)
    # # show(p)
    # p = plot_pathway_all_levels(pathway, graphs=graphs, coloring=Coloring.REACTION)
    # # show(p)
    # p = plot_pathway_all_levels(pathway, graphs=graphs, coloring=Coloring.PATHWAY)
    # show(p)

    # plot_low_level_pathways(10, 20, "../../figures/pathways/", "../../reports/pathways/")
    print("Finished")


if __name__ == '__main__':
    main()
