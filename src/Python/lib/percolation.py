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


def sampleWithLinkPercolation(G, r, P, v=False):
    """
        Calculates relative size of the largest connected component at P network completeness values with r replicates.

        :param G: networkx graph
        :param r: number of replicates for each completeness to average the value of the module size
        :param P: list of network completeness values
        :return: pandas dataframe with columns ['Completeness', 'Relative Size']
        """

    P = sorted(P, reverse=True)

    completenesses = list(itertools.chain.from_iterable(itertools.repeat(P, r)))
    sizes = []
    for replicate in range(r):
        g = copy.deepcopy(G)
        if v: print(f"Replicate {replicate} starting with {g.number_of_nodes()} nodes")

        original_num_nodes = g.number_of_nodes()
        original_num_edges = g.number_of_edges()
        for p in P:
            desired_num_edges = int(p * original_num_edges)
            num_edges_to_remove = g.number_of_edges() - desired_num_edges
            if num_edges_to_remove > 0:
                g.remove_edges_from(random.sample(list(g.edges()), num_edges_to_remove))
            rel_size_of_lcc = len(max(nx.connected_components(g))) * 1 / original_num_nodes
            sizes.append(rel_size_of_lcc)
            if v: print(f"\tCompleteness: {p} \t Graph edges: {desired_num_edges} \t Size of lcc: {rel_size_of_lcc}")

    df = pd.DataFrame(zip(completenesses, sizes), columns=['Completeness', 'Relative Size'])
    df = df.reset_index()
    return df


def sampleWithLinkPercolation(G, r, s=0.1, v=False):
    """
        Calculates relative size of the largest connected component at multiple network completeness values with
        r replicates.

        :param s: step size; percentage of edges (links) to remove at each iteration
        :param G: networkx graph
        :param r: number of replicates for each completeness to average the value of the module size
        :return: pandas dataframe with columns ['Completeness', 'Relative Size']
        """

    completenesses = list(itertools.chain.from_iterable(itertools.repeat(P, r)))
    sizes = []
    for replicate in range(r):
        g = copy.deepcopy(G)
        if v: print(f"Replicate {replicate} starting with {g.number_of_nodes()} nodes")

        original_num_nodes = g.number_of_nodes()
        original_num_edges = g.number_of_edges()
        for p in P:
            desired_num_edges = int(p * original_num_edges)
            num_edges_to_remove = g.number_of_edges() - desired_num_edges
            if num_edges_to_remove > 0:
                g.remove_edges_from(random.sample(list(g.edges()), num_edges_to_remove))
            rel_size_of_lcc = len(max(nx.connected_components(g))) * 1 / original_num_nodes
            sizes.append(rel_size_of_lcc)
            if v: print(f"\tCompleteness: {p} \t Graph edges: {desired_num_edges} \t Size of lcc: {rel_size_of_lcc}")

    df = pd.DataFrame(zip(completenesses, sizes), columns=['Completeness', 'Relative Size'])
    df = df.reset_index()
    return df

if __name__ == '__main__':

    G = nx.fast_gnp_random_graph(83, 0.2, seed=1, directed=False)

    P = getCompletenessValues(20)
    df_link = sampleWithLinkPercolation(P, G, 10, v=True)
    print("Link percolation values:")
    print(df_link)

    plt.plot(df_link['Completeness'], df_link['Relative Size'], 'o')
    plt.title("Link percolation")
    plt.xlabel("Completeness")
    plt.ylabel("Relative size of largest connected component")
    plt.show()

    print(df_link.shape)