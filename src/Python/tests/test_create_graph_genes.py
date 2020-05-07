import pytest

from interaction_network import create_graph
from network_topology_queries import get_reaction_participants_by_pathway


class TestGeneNetworkClass:
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    G_glycolysis = create_graph("R-HSA-9634600", level="genes")

    # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
    G_traffic = create_graph("R-HSA-5627083", level="genes")

    def test_create_graph_wrong_level_raises_exception(self):
        with pytest.raises(Exception):
            get_reaction_participants_by_pathway("R-HSA-9634600", level="banana")

    def test_create_graph_num_edges(self):
        G = TestGeneNetworkClass.G_glycolysis
        print(G.edges)
        assert len(G.edges) == 46

    def test_create_graph_num_vertices(self):
        G = TestGeneNetworkClass.G_glycolysis
        print(G.nodes)
        assert len(G.nodes) == 18

    # Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
    def test_connects_inputs_with_outputs(self):
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
        G = TestGeneNetworkClass.G_glycolysis

        # Input to output interactions for reaction R-HSA-163773:
        assert ("PFKFB1", "ADP") in G.edges
        assert ("ATP", "ADP") in G.edges
        assert ("ATP", "PFKFB1") in G.edges

        # Input to output interactions for reaction R-HSA-163750
        assert ("PFKFB1", "Pi") in G.edges
        assert ("H2O", "Pi") in G.edges
        assert ("H2O", "PFKFB1") in G.edges

        # Input to output interaction for reaction R-HSA-71802
        assert ("Fru(6)P", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("Fru(6)P", "ADP") in G.edges
        assert ("ATP", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("ATP", "ADP") in G.edges

    def test_connects_input_genes_with_small_outputs_not_when_is_same_molecule(self):
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
        G = TestGeneNetworkClass.G_glycolysis

        assert not ("PFKFB1", "PFKFB1") in G.edges

    def test_connects_catalysts_with_outputs(self):
        G = TestGeneNetworkClass.G_glycolysis
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600

        # Catalysts to output interactions for reaction R-HSA-163773:
        assert ("PRKACA", "PFKFB1") in G.edges
        assert ("PRKACA", "ADP") in G.edges
        assert ("PRKACB", "PFKFB1") in G.edges
        assert ("PRKACB", "ADP") in G.edges

        # Catalysts to output interactions for reaction R-HSA-163750
        assert ("PPP2R1A", "PFKFB1") in G.edges
        assert ("PPP2R1A", "Pi") in G.edges
        assert ("PPP2R1B", "PFKFB1") in G.edges
        assert ("PPP2R1B", "Pi") in G.edges
        assert ("PPP2CA", "PFKFB1") in G.edges
        assert ("PPP2CA", "Pi") in G.edges
        assert ("PPP2CB", "PFKFB1") in G.edges
        assert ("PPP2CB", "Pi") in G.edges
        assert ("PPP2R5D", "PFKFB1") in G.edges
        assert ("PPP2R5D", "Pi") in G.edges

        # Catalysts to output interactions for reaction R-HSA-70262
        assert ("PFKFB1", "Fru(6)P") in G.edges
        assert ("PFKFB1", "Pi") in G.edges

        # Catalysts to output interactions for reaction R-HSA-71802
        assert ("PFKFB3", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("PFKFB3", "ADP") in G.edges
        assert ("PFKFB4", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("PFKFB4", "ADP") in G.edges
        assert ("PFKFB1", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("PFKFB1", "ADP") in G.edges
        assert ("PFKFB2", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("PFKFB2", "ADP") in G.edges

    def test_connects_regulators_with_outputs(self):
        # Pathway "RHO GTPases regulate CFTR trafficking" R-HSA-5627083
        G = TestGeneNetworkClass.G_traffic

        # Regulator to output interactions for reaction R-HSA-5627071:
        assert ("RHOQ", "CFTR") in G.edges
        assert ("GOPC", "CFTR") in G.edges
        assert not ("CFTR", "CFTR") in G.edges
        assert ("GTP", "CFTR") in G.edges

        assert ("RHOQ", "GOPC") in G.edges
        assert not ("GOPC", "GOPC") in G.edges
        assert ("CFTR", "GOPC") in G.edges
        assert ("GTP", "GOPC") in G.edges

        assert not ("RHOQ", "RHOQ") in G.edges
        assert ("GOPC", "RHOQ") in G.edges
        assert ("CFTR", "RHOQ") in G.edges
        assert ("GTP", "RHOQ") in G.edges

        assert ("RHOQ", "GTP") in G.edges
        assert ("GOPC", "GTP") in G.edges
        assert ("CFTR", "GTP") in G.edges
        assert not ("GTP", "GTP") in G.edges

    def test_connects_components_of_same_complex(self):
        G = TestGeneNetworkClass.G_glycolysis

        # As part of PP2A-ABdeltaC complex R-HSA-165961
        assert not ("PPP2R1A", "PPP2R1A") in G.edges
        assert ("PPP2R1A", "PPP2R1B") in G.edges
        assert ("PPP2R1B", "PPP2R1A") in G.edges
        assert ("PPP2R1A", "PPP2CA") in G.edges
        assert ("PPP2R1A", "PPP2CB") in G.edges
        assert ("PPP2R1A", "PPP2R5D") in G.edges
        assert ("PPP2R5D", "PPP2R1B") in G.edges
        assert ("PPP2R5D", "PPP2R1A") in G.edges
        assert ("PPP2R5D", "PPP2CB") in G.edges
