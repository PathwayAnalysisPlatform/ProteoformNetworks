import networkx as nx

G = nx.Graph()

G.add_edge("A", "B", edge_color="blue", reaction="r1", complex="TheComplex", weight=1)
G.nodes["A"]["IO"] = True


print(G.nodes(data=True))