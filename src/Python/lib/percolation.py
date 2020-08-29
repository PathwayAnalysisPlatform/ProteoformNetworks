# Library with functions for percolation analysis
import copy
from enum import Enum
import random

import networkx as nx
import numpy as np
import pandas as pd
import itertools


class PercolationType(Enum):
    LINK = 1
    NODE = 2


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
        :return: pandas dataframe with columns ["Completeness", "Size"]
        """

    P = sorted(P, reverse=True)

    completenesses = list(itertools.chain.from_iterable(itertools.repeat(P, r)))
    sizes = []
    replicates = list(itertools.chain.from_iterable(itertools.repeat(x, len(P)) for x in range(r)))
    for replicate in range(r):
        g = copy.deepcopy(G)
        print(f"Replicate {replicate} starting with {g.number_of_nodes()} nodes")

        original_num_nodes = g.number_of_nodes()
        for p in P:
            desired_size = int(p * original_num_nodes)
            print(f"\tCompleteness: {p} \t Graph size: {desired_size}")
            nodes_to_remove = g.number_of_nodes() - desired_size
            if nodes_to_remove > 0:
                g.remove_nodes_from(random.sample(list(g.nodes()), nodes_to_remove))
            size_of_lcc = len(max(nx.connected_components(g))) if g.number_of_nodes() > 0 else 0
            sizes.append(size_of_lcc)
    print("Completeness:")
    print(completenesses)
    print("Sizes:")
    print(sizes)
    print("Replicates:")
    print(replicates)
    df = pd.DataFrame(zip(completenesses, sizes, replicates), columns=['Completeness', 'Relative Size', 'replicate'])
    print(df)
    return df.groupby('Completeness').mean()


def sampleModuleSizeAtCompletenessLevels(P, G, t, r):
    """
    Calculates relative size of the largest connected component at different network completeness values with r
    replicates to average the value.

    :param P: list of network completeness values
    :param G: networkx graph
    :param t: percolation type, link or node
    :param r: number of replicates for each completeness to average the value of the module size
    :return: pandas dataframe with columns ["Completeness", "Size"]
    """

    if t == PercolationType.NODE:
        return sampleWithNodePercolation(P, G, r)

G = nx.complete_graph(9)
P = getCompletenessValues(5)
df = sampleModuleSizeAtCompletenessLevels(P, G, PercolationType.NODE, 3)
print(df)