import networkx as nx

G = nx.Graph()

G.add_edge("A", "B", edge_color="blue", reaction="r1", complex="TheComplex", weight=1)
G["A"]["B"]["IO"] = True


print(G.edges(data=True))