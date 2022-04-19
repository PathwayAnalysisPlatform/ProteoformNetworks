# Library with functions for percolation analysis
import copy
import itertools
import random

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline


def getCompletenessValues(n):
    """
    Returns a list of n float values in range [0, 1]

    :param n: Number of values for completeness
    :return: list of completeness values
    """
    step = 1 / (n - 1)
    return np.arange(0, 1 + step, step)


def sampleLinkPercolationWithPoints(G, r, P, v=False):
    """
        Calculates relative size of the largest connected component at P network completeness values with r replicates.

        :param G: networkx graph
        :param r: number of replicates for each completeness to average the value of the module size
        :param P: list of network completeness values
        :return: pandas dataframe with columns ['Completeness', 'Relative Size']
        """

    P = sorted(P, reverse=True)
    replicates = list()
    completenesses = list(itertools.chain.from_iterable(itertools.repeat(P, r)))
    sizes = []
    for replicate in range(r):
        g = copy.deepcopy(G)
        if v: print(f"Replicate {replicate} starting with {g.number_of_nodes()} nodes")

        original_num_nodes = len(max(nx.connected_components(g)))
        original_num_edges = g.number_of_edges()

        for p in P:
            desired_num_edges = int(p * original_num_edges)
            num_edges_to_remove = g.number_of_edges() - desired_num_edges
            if num_edges_to_remove > 0:
                g.remove_edges_from(random.sample(list(g.edges()), num_edges_to_remove))
            rel_size_of_lcc = len(max(nx.connected_components(g))) / original_num_nodes
            sizes.append(rel_size_of_lcc)
            if v: print(f"\tCompleteness: {p} \t Graph edges: {desired_num_edges} \t Size of lcc: {rel_size_of_lcc}")
            replicates.append(replicate)

    df = pd.DataFrame(zip(completenesses, sizes, replicates), columns=['Completeness', 'Relative Size', 'Replicate'])
    df = df.reset_index()
    return df


def sampleLinkPercolationWithPercentages(G, r, s=0.1, v=False):
    """
        Calculates relative size of the largest connected component at multiple network completeness values with
        r replicates. This method makes sure to return the completeness value for Relative Size of Lcc at 0.5 for each
        replicate.

        :param s: step size; percentage of edges (links) to remove at each iteration. Float value in [0, 1]
        :param G: networkx graph
        :param r: number of replicates for each completeness to average the value of the module size
        :return: pandas dataframe with columns ['Completeness', 'Relative Size']
        """

    replicates = list()
    completenesses = list()
    sizes = []
    for replicate in range(r):
        g = copy.deepcopy(G)
        if v: print(f"Replicate {replicate} starting with {g.number_of_nodes()} nodes")

        original_size_lcc = len(max(nx.connected_components(g)))
        original_num_edges = g.number_of_edges()
        while g.number_of_edges() > 0:
            desired_num_edges = int((1 - s) * g.number_of_edges())  # How many edges should be after this iteration
            num_edges_to_remove = max(g.number_of_edges() - desired_num_edges, 0)
            g.remove_edges_from(random.sample(list(g.edges()), num_edges_to_remove))  # Remove the edges
            rel_size_lcc = len(max(nx.connected_components(g))) / original_size_lcc
            sizes.append(rel_size_lcc)
            p = g.number_of_edges() / original_num_edges
            completenesses.append(p)
            replicates.append(replicate)
            if v: print(f"\tCompleteness: {p} \t Graph edges: {desired_num_edges} \t Size of lcc: {rel_size_lcc}")
        print(f"Done replicate {replicate}")

    df = pd.DataFrame(zip(completenesses, sizes, replicates), columns=['Completeness', 'Relative Size', 'Replicate'])
    df = df.reset_index()
    return df


if __name__ == '__main__':
    G_1 = nx.fast_gnp_random_graph(50, 0.2, seed=1, directed=False)
    G_2 = nx.fast_gnp_random_graph(150, 0.2, seed=1, directed=False)

    # P = getCompletenessValues(5)
    df_1 = sampleLinkPercolationWithPercentages(G_1, 20, 0.1)
    df_2 = sampleLinkPercolationWithPercentages(G_2, 20, 0.1)
    print("Link percolation values 1:")
    print(df_1)
    print("Link percolation values 2:")
    print(df_1.groupby('Replicate'))

    plt.style.use('seaborn-darkgrid')
    for label, sub_df in df_1.groupby('Replicate'):
        print(label)
        plt.plot(sub_df['Completeness'], sub_df['Relative Size'], marker='', linewidth=0.5, alpha=0.3, color='green')

    for label, sub_df in df_2.groupby('Replicate'):
        print(label)
        plt.plot(sub_df['Completeness'], sub_df['Relative Size'], marker='', linewidth=0.5, alpha=0.3, color='blue')

    df_1_mean = df_1.groupby('Completeness').mean()
    df_1_mean = df_1_mean.reset_index()
    df_1_mean = df_1_mean.drop(columns=['Replicate', 'index'])
    print(df_1_mean.head())

    x_new = np.linspace(0, 1, 300)
    spl = make_interp_spline(df_1_mean['Completeness'], df_1_mean['Relative Size'])
    y_new = spl(x_new)
    plt.plot(x_new, y_new, marker='', linewidth=2, alpha=1, color='green', label='Genes')
    #
    df_2_mean = df_2.groupby('Completeness').mean()
    df_2_mean = df_2_mean.reset_index()
    df_2_mean = df_2_mean.drop(columns=['Replicate', 'index'])
    print(df_2_mean)
    # plt.plot(df_2_mean['Completeness'], df_2_mean['Relative Size'], marker='', linewidth=2, alpha=1, color='blue', label='Proteoforms')

    x_new = np.linspace(0, 1, 300)
    spl = make_interp_spline(df_2_mean['Completeness'], df_2_mean['Relative Size'])
    y_new = spl(x_new)
    plt.plot(x_new, y_new, marker='', linewidth=2, alpha=1, color='blue', label='Proteoforms')

    plt.legend()
    plt.title("Link percolation")
    plt.xlabel("Completeness")
    plt.ylabel("Relative size of largest connected component")
    plt.show()
