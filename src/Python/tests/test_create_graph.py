import os
import shutil
from collections import namedtuple
from pathlib import Path

import pandas as pd
import pytest

from interaction_network import create_graph, create_graphs, Pathway_graphs
from network_topology_queries import get_pathways


@pytest.fixture
def delete_graphs_path(scope="function"):
    graphs_path = __name__ + "/"

    yield
    # Delete graphs if created
    print(f"Deleting graphs at {graphs_path}")
    dirpath = Path(graphs_path)
    if dirpath.exists() and dirpath.is_dir():
        shutil.rmtree(dirpath)


class TestCreate_graph():

    def test_create_graph_edges_have_color_attribute_existing_pathway(self):
        G = create_graph("R-HSA-70171", level="proteins", showSmallMolecules=True)

        print(f"Edges data: {G.edges.data()}")

        for edge, attrs in G.edges.items():
            print(f"Checking edge: {edge}")
            assert "edge_color" in attrs.keys()
            print(f" it has color {G.edges[edge]['edge_color']}")

        if os.path.exists("R-HSA-70171_proteins_edge_list"):
            os.remove("R-HSA-70171_proteins_edge_list")

    def test_create_graph_edges_have_color_attribute_non_existing_pathway(self):
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


    def test_create_graphs_returns_empty_result_when_no_pathway_arguments(self):
        result = create_graphs(pd.DataFrame)
        assert type(result) == list
        assert len(result) == 0



    def test_create_graphs_returns_the_correct_number_of_results(self, delete_graphs_path):
        pathways = get_pathways()
        graphs_path = __name__ + "/"
        result = create_graphs(pathways[:3], graphs_path)
        assert type(result) == list
        assert len(result) == 3

        for entry in result:
            assert type(entry) == Pathway_graphs
            assert "genes" in entry._fields
            assert "genes_no_small_molecules" in entry._fields
            assert "proteins" in entry._fields
            assert "proteins_no_small_molecules" in entry._fields
            assert "proteoforms" in entry._fields
            assert "proteoforms_no_small_molecules" in entry._fields


    def test_create_graphs_stores_the_graphs(self, delete_graphs_path):
        pathways = get_pathways()
        pathways = pathways[:3]
        graphs_path = __name__ + "/"

        # Check graphs do not exist
        assert not os.path.exists(graphs_path)
        for pathway in pathways['stId']:
            assert not os.path.exists(graphs_path + pathway)

        create_graphs(pathways, graphs_path)

        # Check graphs exist
        assert os.path.exists(graphs_path)
        for pathway in pathways['stId']:
            assert os.path.exists(graphs_path + pathway+ "_genes_edge_list")
            assert os.path.exists(graphs_path + pathway+ "_proteins_edge_list")
            assert os.path.exists(graphs_path + pathway+ "_proteoforms_edge_list")
            assert os.path.exists(graphs_path + pathway+ "_genes_no_small_molecules_edge_list")
            assert os.path.exists(graphs_path + pathway+ "_proteins_no_small_molecules_edge_list")
            assert os.path.exists(graphs_path + pathway+ "_proteoforms_no_small_molecules_edge_list")