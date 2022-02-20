import matplotlib.pyplot as plt
import networkx as nx

G = nx.complete_graph(4)
pos = nx.spring_layout(G, seed=1)
print(pos)
fixed = [1]
pos = nx.spring_layout(G, pos=pos, fixed=fixed, seed=3)
print(pos)

index = 1
for node in G.nodes:
    G.nodes[node]['pathways'] = set()
    for i in range(index):
        G.nodes[node]['pathways'].add(chr(ord('A') + i))
    index += 1

fig = plt.figure(figsize=(50, 50))
ax = plt.axes([0, 0, 1, 1])
ax.set_aspect('equal', anchor='C')
nx.draw_networkx_edges(G, pos, ax=ax, width=5, style='solid')

plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)

trans = ax.transData.transform
trans2 = fig.transFigure.inverted().transform

piesize = 0.1
p2 = piesize / 2.0
for n in G:
    print("********************")
    print(f"Node {n}")
    print(f"Position: {pos[n]}")
    xx, yy = trans(pos[n])  # figure coordinates
    print(f"Figure coordinates: {(xx, yy)}")

    xa, ya = trans2((xx, yy))  # axes coordinates
    print(f"Axes coordinates: {(xa, ya)}")

    # a = plt.axes([xa - p2, ya - p2, piesize, piesize])
    # r = [0.5, 0.5, piesize, piesize]
    r = ([xa - p2, ya - p2, piesize, piesize])
    print(f"Rectangle: {r}")
    a = plt.axes(r)
    a.set_aspect('equal')

    ratio = 100 / len(G.nodes[n]['pathways'])
    fracs = [ratio] * len(G.nodes[n]['pathways'])
    a.pie(fracs)

plt.show()
