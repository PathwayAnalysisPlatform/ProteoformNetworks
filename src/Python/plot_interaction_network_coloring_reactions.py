from math import pi

import networkx as nx
import pandas as pd
from bokeh.io import show
from bokeh.models import (BoxZoomTool, HoverTool,
                          Range1d, ResetTool)
from bokeh.models.graphs import from_networkx
from bokeh.palettes import cividis
# Prepare Data
from bokeh.plotting import figure
from bokeh.transform import cumsum

from interaction_network import create_graph

v = False

G = create_graph("R-HSA-70171", "genes", True)

# SAME_CLUB_COLOR, DIFFERENT_CLUB_COLOR = "black", "red"
# edge_attrs = {}
#
# for start_node, end_node, _ in G.edges(data=True):
#     edge_color = SAME_CLUB_COLOR if G.nodes[start_node]["club"] == G.nodes[end_node]["club"] else DIFFERENT_CLUB_COLOR
#     edge_attrs[(start_node, end_node)] = edge_color
#
# nx.set_edge_attributes(G, edge_attrs, "edge_color")

f = figure(plot_width=1000, plot_height=800,
           x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1))
f.title.text = "Graph Interaction Demonstration"
f.axis.visible = False
f.axis.axis_label = None
f.grid.grid_line_color = None
f.legend.title = 'Reactions'

graph_renderer = from_networkx(G, nx.spring_layout, scale=1, center=(0, 0))

if v:
    print(graph_renderer._property_values['layout_provider']._property_values['graph_layout'])
pos = graph_renderer._property_values['layout_provider']._property_values['graph_layout']

# Display edges
for start, end in G.edges:
    if v:
        print(f"Edge from {start} -- {end}")
    x_start, y_start = pos[start]
    x_end, y_end = pos[end]
    f.line([x_start, x_end], [y_start, y_end], line_width=2, line_color="black")

# Display nodes
reactions = set()
for k in G.nodes:
    for reaction in G.nodes[k]['reactions']:
        reactions.add(reaction)

r = 0.05
wedges = []
for k in G.nodes:
    x, y = pos[k]
    if v:
        print(f"Node {k} has position {x}, {y}")
        print(f"Belongs to reactions {G.nodes[k]['reactions']}")
    data = pd.Series({reaction: 1 if reaction in G.nodes[k]['reactions'] else 0 for reaction in reactions}).reset_index(
        name='value').rename(
        columns={'index': 'reaction'})
    data['id'] = k
    data['type'] = G.nodes[k]['type']
    data['entity_color'] = G.nodes[k]['entity_color']
    data['angle'] = data['value'] / data['value'].sum() * 2 * pi
    data['reaction_color'] = cividis(len(reactions))
    w = f.wedge(x=x, y=y, radius=r,
                start_angle=cumsum('angle', include_zero=True),
                end_angle=cumsum('angle'),
                line_color="white", fill_color='reaction_color', legend_field='reaction', source=data)
    wedges.append(w)

TOOLTIPS = """
            <div>
                <span style="font-size: 10px; font-weight: bold;">Name: @id</span>
            </div>
            <div>
                <span style="font-size: 10px; font-weight: bold; color: @entity_color;">Type: @type</span>
            </div>
        """
node_hover_tool = HoverTool(tooltips=TOOLTIPS)
f.add_tools(node_hover_tool, BoxZoomTool(), ResetTool())

show(f)
