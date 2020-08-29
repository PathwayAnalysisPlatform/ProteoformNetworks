# Library with functions for percolation analysis
import copy
import itertools
import random

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


def getCompletenessValues(n):
    """
    Returns a list of n float values in range [0, 1]

    :param n: Number of values for completeness
    :return: list of completeness values
    """
    step = 1 / (n - 1)
    return np.arange(0, 1 + step, step)


def sampleWithNodePercolation(P, G, r):
    """
        Calculates relative size of the largest connected component at different network completeness values with r
        replicates to average the value.

        :param P: list of network completeness values
        :param G: networkx graph
        :param r: number of replicates for each completeness to average the value of the module size
        :return: pandas dataframe with columns ['Completeness', 'Relative Size']
        """

    P = sorted(P, reverse=True)

    completenesses = list(itertools.chain.from_iterable(itertools.repeat(P, r)))
    sizes = []
    # replicates = list(itertools.chain.from_iterable(itertools.repeat(x, len(P)) for x in range(r)))
    for replicate in range(r):
        g = copy.deepcopy(G)
        # print(f"Replicate {replicate} starting with {g.number_of_nodes()} nodes")

        original_num_nodes = g.number_of_nodes()
        for p in P:
            desired_size = int(p * original_num_nodes)
            print(f"\tCompleteness: {p} \t Graph size: {desired_size}")
            nodes_to_remove = g.number_of_nodes() - desired_size
            if nodes_to_remove > 0:
                g.remove_nodes_from(random.sample(list(g.nodes()), nodes_to_remove))
            size_of_lcc = len(max(nx.connected_components(g))) if g.number_of_nodes() > 0 else 0
            sizes.append(size_of_lcc)
    # print("Completeness:")
    # print(completenesses)
    # print("Sizes:")
    # print(sizes)
    df = pd.DataFrame(zip(completenesses, sizes), columns=['Completeness', 'Relative Size'])
    # print(df)
    df = df.groupby('Completeness').mean()
    df = df.reset_index()
    # print(df)
    return df

def sampleWithLinkPercolation(P, G, r):
    """
        Calculates relative size of the largest connected component at different network completeness values with r
        replicates to average the value.

        :param P: list of network completeness values
        :param G: networkx graph
        :param r: number of replicates for each completeness to average the value of the module size
        :return: pandas dataframe with columns ['Completeness', 'Relative Size']
        """

    P = sorted(P, reverse=True)

    completenesses = list(itertools.chain.from_iterable(itertools.repeat(P, r)))
    sizes = []
    for replicate in range(r):
        g = copy.deepcopy(G)
        # print(f"Replicate {replicate} starting with {g.number_of_nodes()} nodes")

        original_num_edges = g.number_of_edges()
        for p in P:
            desired_num_edges = int(p * original_num_edges)
            num_edges_to_remove = g.number_of_edges() - desired_num_edges
            if num_edges_to_remove > 0:
                g.remove_edges_from(random.sample(list(g.edges()), num_edges_to_remove))
            size_of_lcc = len(max(nx.connected_components(g)))
            sizes.append(size_of_lcc)
            print(f"\tCompleteness: {p} \t Graph edges: {desired_num_edges} \t Size of lcc: {size_of_lcc}")
    # print("Completeness:")
    # print(completenesses)
    # print("Sizes:")
    # print(sizes)
    df = pd.DataFrame(zip(completenesses, sizes), columns=['Completeness', 'Relative Size'])
    print(df)
    df = df.groupby('Completeness').mean()
    df = df.reset_index()
    print(df)
    return df

fig, (ax1, ax2) = plt.subplots(2, 1)
fig.subplots_adjust(hspace=0.8, wspace=0.4)

G = nx.fast_gnp_random_graph(100, 0.2, seed=1, directed=False)

P = getCompletenessValues(20)
df_node = sampleWithNodePercolation(P, G, 3)
print("Node percolation values:")
print(df_node)

df_link = sampleWithLinkPercolation(P, G, 3)
print("Link percolation values:")
print(df_link)

ax1.set_title("Node percolation")
ax1.plot(df_node['Completeness'], df_node['Relative Size'], 'o-')
ax1.set_xlabel('Completeness')
ax1.set_ylabel('Relative Size')
# df.plot(kind='scatter', x='Completeness', y='Relative Size')

ax2.set_title("Link percolation")
ax2.plot(df_link['Completeness'], df_link['Relative Size'], 'o-')
ax2.set_xlabel('Completeness')
ax2.set_ylabel('Relative Size')


plt.show()
