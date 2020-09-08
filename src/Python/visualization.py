# %%
from enum import Enum
from math import pi

import networkx as nx
import pandas as pd
from bokeh.io import output_file, save
from bokeh.io import show
from bokeh.layouts import layout
from bokeh.models import (BoxZoomTool, HoverTool,
                          Range1d, ResetTool, ColumnDataSource, Legend)
from bokeh.models import Div
from bokeh.palettes import Colorblind
# Prepare Data
from bokeh.plotting import figure
from bokeh.transform import cumsum

from config import get_entity_color, COLOR_IO, COLOR_CO, COLOR_RO, COLOR_CC
from interaction_network import create_graph, merge_graphs
from interaction_network import create_pathway_graphs
from network_topology_queries import get_pathway_name, get_low_level_pathways


class Coloring(Enum):
    ENTITY_TYPE = "Entity Type"
    REACTION = "Reaction"
    PATHWAY = "Pathway"


def plot_pathway(pathway, level="proteins", sm=True, coloring=Coloring.ENTITY_TYPE, v=True):
    """
    Basic plotting function for a pathway. Shows only one granularity level and with/without small molecules.

    :param pathway: string stId of the pathway
    :param level: string {"genes", "proteins", "proteoforms"}
    :param sm: bool, show small molecules in the network
    :param coloring: Coloring instance
    :param v: bool v
    :return: bokeh plot
    """
    G = create_graph(pathway, level, sm, "", v)
    P = plot_interaction_network(G, coloring)
    return P


def getDefaultPalette(n):
    if n == 1:
        return [Colorblind[3][0]]
    elif n == 2:
        return [Colorblind[3][0], Colorblind[3][1]]
    else:
        return Colorblind[n]


def plot_interaction_network(G, coloring=Coloring.ENTITY_TYPE, v=False, **kwargs):
    """
    Plots the interaction network represented by Graph G
    :param G: networkx graph, it must have attributes "stId" with the Pathway id
    :param coloring: Type of coloring for the nodes and links, by entity type, reaction or pathway
    :param v: Verbose
    :return: the figure
    """

    plot_width = kwargs['plot_width'] if 'plot_width' in kwargs else 600
    plot_height = kwargs['plot_height'] if 'plot_height' in kwargs else 325
    toolbar_location = kwargs['toolbar_location'] if 'toolbar_location' in kwargs else "below"

    f = figure(plot_width=plot_width, plot_height=plot_height,
               x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1),
               toolbar_location=toolbar_location)
    f.title.text = kwargs['title'] if 'title' in kwargs else G.graph['level'].title()
    f.title.text_font_size = '12pt'
    f.title.align = 'center'
    f.axis.visible = False
    f.axis.axis_label = None
    f.grid.grid_line_color = None
    f.outline_line_width = kwargs['outline_line_width'] if 'outline_line_width' in kwargs else 0

    pos = kwargs['pos'] if 'pos' in kwargs else nx.spring_layout(G)
    if v:
        print(pos)

    # Display edges
    for start, end in G.edges:
        if v:
            print(f"Edge from {start} -- {end}")
        x_start, y_start = pos[start]
        x_end, y_end = pos[end]
        f.line([x_start, x_end], [y_start, y_end], line_width=2, line_color="black")

    if coloring == Coloring.ENTITY_TYPE:

        x_values = []
        y_values = []
        entity_colors = []
        types = []
        ids = []

        x_values_se = []
        y_values_se = []
        entity_colors_se = []
        types_se = []
        ids_se = []

        for id in G.nodes:
            if G.nodes[id]['type'] == "SimpleEntity":
                x_values_se.append(pos[id][0])
                y_values_se.append(pos[id][1])
                entity_colors_se.append(G.nodes[id]['entity_color'])
                types_se.append(G.nodes[id]['type'])
                ids_se.append(G.nodes[id]['id'])
            else:
                x_values.append(pos[id][0])
                y_values.append(pos[id][1])
                entity_colors.append(G.nodes[id]['entity_color'])
                types.append(G.nodes[id]['type'])
                ids.append(G.nodes[id]['id'])

        data = {'id': ids, 'type': types, 'x_values': x_values, 'y_values': y_values, 'entity_color': entity_colors}

        data_se = {'id': ids_se, 'type': types_se, 'x_values': x_values_se, 'y_values': y_values_se,
                   'entity_color': entity_colors_se}
        if len(data_se['id']) > 0:
            f.oval(x='x_values', y='y_values', color='entity_color', height=0.1, width=0.17, line_color='black',
                   source=ColumnDataSource(data=data_se), name='node', legend_label="Small Molecules")
        f.rect(x='x_values', y='y_values', color='entity_color', height=0.08, width=0.11, line_color='black',
               source=ColumnDataSource(data=data), name='node', legend_label=G.graph["level"].title())
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
            if v:
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
                f.legend.visible = False
                f.add_layout(legend, 'right')
            elif kwargs['legend_location']:
                f.legend.location = kwargs['legend_location']
            else:
                f.legend.visible = False

    TOOLTIPS = """
            <div>
                <span style="font-size: 12px; font-weight: bold;">Id: @id</span>
            </div>
            <div>
                <span style="font-size: 12px; font-weight: bold; color: @entity_color;">Type: @type</span>
            </div>
        """
    node_hover_tool = HoverTool(tooltips=TOOLTIPS, names=['node'])
    f.add_tools(node_hover_tool, BoxZoomTool(), ResetTool())

    # output_file(f"{G.graph['level']}_network.html")
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


def plot_pathway_all_levels(pathway, figures_path="../../figures/pathways/", graphs_path="../../reports/pathways/",
                            coloring=Coloring.ENTITY_TYPE,
                            graphs=None,
                            v=False):
    name = get_pathway_name(pathway)['Name'][0]
    if len(name) == 0:
        return

    graphs_sm, graphs_no_sm = graphs or create_pathway_graphs(pathway, graphs_path)

    pos_all = [nx.spring_layout(g) for g in graphs_sm]
    node_set_unasigned = [set(g.nodes) for g in graphs_sm]
    fixed_all = [set(), set(), set()]
    titles_sm = ['A', 'B', 'C']
    titles_no_sm = ['D', 'E', 'F']
    plot_size = 400
    legend_location_all = [None, None, 'right']
    plot_widths = [plot_size, plot_size, plot_size + 160]
    if coloring == Coloring.ENTITY_TYPE:
        legend_location_all = ['top_right', 'top_right', 'top_right']
        plot_widths = [plot_size, plot_size, plot_size]

    # For each node in the gene graph, assign it's position to the first node with that id in the protein network
    for i in reversed(range(2)):
        for node in graphs_sm[i + 1].nodes:  # For each node in the larger network
            prevId = graphs_sm[i + 1].nodes[node]['prevId']  # Get the predecesor
            if prevId in node_set_unasigned[i]:  # If predecesor has no position yet
                pos_all[i][prevId] = pos_all[i + 1][node]
                node_set_unasigned[i].remove(prevId)  # Mark predecesor as assigned to a position
                fixed_all[i].add(prevId)
        pos_all[i] = nx.spring_layout(graphs_sm[i], pos=pos_all[i], fixed=fixed_all[i])

    figures_sm = [
        plot_interaction_network(graphs_sm[i], coloring=coloring, pos=pos_all[i], plot_width=plot_widths[i],
                                 plot_height=plot_size,
                                 toolbar_location=None, title=titles_sm[i], legend_location=legend_location_all[i]) for
        i in range(3)
    ]
    if coloring != Coloring.ENTITY_TYPE:
        legend_location_all = [None, None, None]
        plot_widths = [plot_size, plot_size, plot_size]
    figures_no_sm = [
        plot_interaction_network(graphs_no_sm[i], coloring=coloring, pos=pos_all[i], plot_width=plot_widths[i],
                                 plot_height=plot_size,
                                 toolbar_location=None, title=titles_no_sm[i], legend_location=legend_location_all[i])
        for i in range(3)
    ]

    if v:
        import os
        print("In function create_graph: ", os.getcwd())
        print(f"Output path is: {figures_path}{pathway}_network.html")
    output_file(f"{figures_path}{pathway}_{coloring.name}_network.html")

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
        figures_sm,
        figures_no_sm
    ])

    # show(l)
    save(l)
    return l


def plot_low_level_pathways(min_size=5, max_size=20, figures_path="figures/pathways/", graphs_path="reports/pathways/"):
    pathways = get_low_level_pathways()
    for pathway in pathways['stId']:
        name = get_pathway_name(pathway)
        name = name.iloc[0]['Name']

        G = create_graph(pathway, "genes", True, graphs_path, save=False)

        if min_size <= len(G.nodes) <= max_size:
            print(f"-- Plotting pathway \"{pathway}\" with size {len(G.nodes)}")
            plot_pathway_all_levels(pathway, figures_path, graphs_path)
        # else:
        #     print(f"Skipping pathway: \"{pathway}\" with size {len(G.nodes)}")


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
    graphs = [create_graph(pathway, level, sm) for pathway in pathways]
    full_graph = merge_graphs(graphs)
    f = plot_interaction_network(full_graph, Coloring.PATHWAY)
    show(f)
    return f


def main():
    pathway = "R-HSA-9634600"
    graphs = create_pathway_graphs(pathway)
    p = plot_pathway_all_levels(pathway, graphs=graphs, coloring=Coloring.ENTITY_TYPE)
    # show(p)
    p = plot_pathway_all_levels(pathway, graphs=graphs, coloring=Coloring.REACTION)
    # show(p)
    p = plot_pathway_all_levels(pathway, graphs=graphs, coloring=Coloring.PATHWAY)
    # show(p)

    # plot_low_level_pathways(10, 20, "../../figures/pathways/", "../../reports/pathways/")
    print("Finished")


if __name__ == '__main__':
    main()
