import networkx as nx
import pytest

from config import genes
from lib.networks import get_nodes_and_adjacent


@pytest.fixture(scope="session")
def path_graph():
    G = nx.path_graph(4)

    G.nodes[0]['type'] = genes
    G.nodes[1]['type'] = genes
    G.nodes[2]['type'] = genes
    G.nodes[3]['type'] = "SimpleEntity"
    return G


def test_get_nodes_and_adjacent(path_graph):
    nodes = [1, 2]
    node_set = get_nodes_and_adjacent(nodes, path_graph)
    assert 0 not in node_set
    assert 1 in node_set
    assert 2 in node_set
    assert 3 in node_set


def test_get_nodes_and_adjacent_ignores_non_existent(path_graph):
    nodes = [2, 5]
    node_set = get_nodes_and_adjacent(nodes, path_graph)
    assert 0 not in node_set
    assert 1 not in node_set
    assert 2 in node_set
    assert 3 in node_set
    assert 5 not in node_set
