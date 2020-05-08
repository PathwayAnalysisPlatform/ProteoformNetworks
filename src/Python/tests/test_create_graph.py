from unittest import TestCase

from interaction_network import create_graph


class TestCreate_graph(TestCase):

    def test_create_graph_edges_have_color_attribute_existing_pathway(self):
        G = create_graph("R-HSA-70171", level="proteins", showSmallMolecules=True, verbose=True)

        print(f"Edges data: {G.edges.data()}")

        for edge, attrs in G.edges.items():
            print(f"Checking edge: {edge}")
            assert "edge_color" in attrs.keys()
            print(f" it has color {G.edges[edge]['edge_color']}")


    def test_create_graph_edges_have_color_attribute_non_existing_pathway(self):
        """
        Its enough that the function does not crash and returns an empty graph
        :return:
        """
        G = create_graph("apple", "genes", False, False)

        print(f"Edges data: {G.edges.data()}")

        for edge, attrs in G.edges.items():
            print(f"Checking edge: {edge}")
            assert "edge_color" in attrs.keys()
            print(f" it has color {G.edges[edge]['edge_color']}")

