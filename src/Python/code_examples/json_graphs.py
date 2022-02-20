from networkx.readwrite import json_graph
import json
import networkx as nx

G = nx.Graph([('A', 'B')])
G.nodes['A']['attr1'] = "valueReally"
G.nodes['B']['attr1'] = "value2"

data = json_graph.node_link_data(G)

print(G.nodes['B']['attr1'])
print(G.nodes['A']['attr1'])

# print(data)

with open("json_file.json", 'w' ) as outfile:
    json.dump(data, outfile)


H = json_graph.node_link_graph(data)

print(H.nodes['A']['attr1'])

