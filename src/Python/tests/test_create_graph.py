import os
from collections import namedtuple

import pandas as pd
import pytest

from interaction_network import create_graph, create_graphs, Pathway_graphs
from network_topology_queries import get_pathways


@pytest.fixture
def glycolysis_graph(tmpdir):
    return create_graph("R-HSA-70171", level="proteins", sm=True, graphs_path=tmpdir)


@pytest.fixture(scope="session")
def some_graphs(tmpdir_factory):
    pathways = get_pathways()
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    Result = namedtuple("Result", "graphs graphs_path pathways")
    return Result(create_graphs(pathways[:3], graphs_path), graphs_path, pathways)


# def test_create_graph_edges_have_color_attribute_existing_pathway(glycolysis_graph):
#     G = glycolysis_graph
#
#     print(f"Edges data: {G.edges.data()}")
#
#     for edge, attrs in G.edges.items():
#         print(f"Checking edge: {edge}")
#         assert "edge_color" in attrs.keys()
#         print(f" it has color {G.edges[edge]['edge_color']}")
#
#     if os.path.exists("R-HSA-70171_proteins_edge_list"):
#         os.remove("R-HSA-70171_proteins_edge_list")


def test_create_graph_edges_have_color_attribute_non_existing_pathway():
    """
    Its enough that the function does not crash and returns an empty graph
    :return:
    """
    G = create_graph("apple", "genes", False)

    print(f"Edges data: {G.edges.data()}")

    for edge, attrs in G.edges.items():
        print(f"Checking edge: {edge}")
        assert "edge_color" in attrs.keys()
        print(f" it has color {G.edges[edge]['edge_color']}")


def test_create_graphs_returns_empty_result_when_no_pathway_arguments():
    result = create_graphs(pd.DataFrame)
    assert type(result) == list
    assert len(result) == 0


# def test_create_graphs_returns_the_correct_number_of_results(some_graphs):
#
#     assert type(some_graphs.graphs) == list
#     assert len(some_graphs.graphs) == 3
#
#     for entry in some_graphs.graphs:
#         assert type(entry) == Pathway_graphs
#         assert "genes" in entry._fields
#         assert "genes_no_small_molecules" in entry._fields
#         assert "proteins" in entry._fields
#         assert "proteins_no_small_molecules" in entry._fields
#         assert "proteoforms" in entry._fields
#         assert "proteoforms_no_small_molecules" in entry._fields


# def test_create_graphs_stores_the_graphs(some_graphs):
#     # Check graphs exist
#     assert os.path.exists(some_graphs.graphs_path)
#     for pathway in some_graphs.pathways:
#         assert os.path.exists(os.path.join(str(some_graphs.graphs_path), pathway + "_genes_edge_list"))
#         assert os.path.exists(os.path.join(str(some_graphs.graphs_path), pathway + "_proteins_edge_list"))
#         assert os.path.exists(os.path.join(str(some_graphs.graphs_path), pathway + "_proteoforms_edge_list"))
#         assert os.path.exists(os.path.join(str(some_graphs.graphs_path), pathway + "_genes_no_small_molecules_edge_list"))
#         assert os.path.exists(os.path.join(str(some_graphs.graphs_path), pathway + "_proteins_no_small_molecules_edge_list"))
#         assert os.path.exists(os.path.join(str(some_graphs.graphs_path), pathway + "_proteoforms_no_small_molecules_edge_list"))
