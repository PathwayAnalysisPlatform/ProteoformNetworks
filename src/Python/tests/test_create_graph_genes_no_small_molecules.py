from interaction_network import create_graph
from network_topology_queries import get_reaction_participants_by_pathway


# Test graph creation with genes and small molecules disabled
class TestGeneNetworkNoSmallMoleculesClass:
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    G_glycolysis = create_graph("R-HSA-9634600", level="genes", showSmallMolecules=False)

    # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
    G_traffic = create_graph("R-HSA-5627083", level="genes", showSmallMolecules=False)

    # Pathway "FRS-mediated FGFR3 signaling"
    G_signal = create_graph("R-HSA-5654706", level="genes", showSmallMolecules=False)

    def test_create_graph_num_edges(self):
        G = TestGeneNetworkNoSmallMoleculesClass.G_glycolysis
        print(G.edges)
        assert len(G.edges) == 18

    def test_create_graph_num_vertices(self):
        G = TestGeneNetworkNoSmallMoleculesClass.G_glycolysis
        print(G.nodes)
        assert len(G.nodes) == 12

    # Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
    def test_connects_inputs_with_outputs(self):
        G = TestGeneNetworkNoSmallMoleculesClass.G_glycolysis
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600

        # Input to output interactions for reaction R-HSA-163750
        assert not ("PFKFB1", "Pi") in G.edges
        assert not ("PFKFB1", "PFKFB1") in G.edges
        assert not ("H2O", "Pi") in G.edges
        assert not ("H2O", "PFKFB1") in G.edges

        # Input to output interactions for reaction R-HSA-163773:
        assert not ("PFKFB1", "ADP") in G.edges
        assert not ("ATP", "ADP") in G.edges
        assert not ("ATP", "PFKFB1") in G.edges

        # Input to output interaction for reaction R-HSA-71802
        assert not ("Fru(6)P", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("Fru(6)P", "ADP") in G.edges
        assert not ("ATP", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("ATP", "ADP") in G.edges

    def test_connects_catalysts_with_outputs(self):
        G = TestGeneNetworkNoSmallMoleculesClass.G_glycolysis
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600

        # Catalysts to output interactions for reaction R-HSA-163773:
        assert ("PRKACA", "PFKFB1") in G.edges
        assert ("PRKACG", "PFKFB1") in G.edges
        assert ("PRKACB", "PFKFB1") in G.edges
        assert not ("PRKACA", "ADP") in G.edges
        assert not ("PRKACB", "ADP") in G.edges

        # Catalysts to output interactions for reaction R-HSA-163750
        assert ("PPP2R1A", "PFKFB1") in G.edges
        assert not ("PPP2R1A", "Pi") in G.edges
        assert ("PPP2R1B", "PFKFB1") in G.edges
        assert not ("PPP2R1B", "Pi") in G.edges
        assert ("PPP2CA", "PFKFB1") in G.edges
        assert not ("PPP2CA", "Pi") in G.edges
        assert ("PPP2CB", "PFKFB1") in G.edges
        assert not ("PPP2CB", "Pi") in G.edges
        assert ("PPP2R5D", "PFKFB1") in G.edges
        assert not ("PPP2R5D", "Pi") in G.edges

        # Catalysts to output interactions for reaction R-HSA-70262
        assert not ("PFKFB1", "Fru(6)P") in G.edges
        assert not ("PFKFB1", "Pi") in G.edges

        # Catalysts to output interactions for reaction R-HSA-71802
        assert not ("PFKFB3", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("PFKFB3", "ADP") in G.edges
        assert not ("PFKFB4", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("PFKFB4", "ADP") in G.edges
        assert not ("PFKFB1", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("PFKFB1", "ADP") in G.edges
        assert not ("PFKFB2", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("PFKFB2", "ADP") in G.edges

    def test_connects_regulators_with_outputs(self):
        # Pathway "RHO GTPases regulate CFTR trafficking" R-HSA-5627083
        G = TestGeneNetworkNoSmallMoleculesClass.G_traffic

        # Regulator to output interactions for reaction R-HSA-5627071:
        assert ("RHOQ", "CFTR") in G.edges
        assert ("GOPC", "CFTR") in G.edges
        assert not ("CFTR", "CFTR") in G.edges
        assert not ("GTP", "CFTR") in G.edges

        assert ("RHOQ", "GOPC") in G.edges
        assert not ("GOPC", "GOPC") in G.edges
        assert ("CFTR", "GOPC") in G.edges
        assert not ("GTP", "GOPC") in G.edges

        assert not ("RHOQ", "RHOQ") in G.edges
        assert ("GOPC", "RHOQ") in G.edges
        assert ("CFTR", "RHOQ") in G.edges
        assert not ("GTP", "RHOQ") in G.edges

        assert not ("RHOQ", "GTP") in G.edges
        assert not ("GOPC", "GTP") in G.edges
        assert not ("CFTR", "GTP") in G.edges
        assert not ("GTP", "GTP") in G.edges

    def test_connects_regulators_with_outputs_2(self):
        # Pathway "FRS-mediated FGFR3 signaling" R-HSA-5654706
        G = TestGeneNetworkNoSmallMoleculesClass.G_signal

        # In reaction R-HSA-5654408
        assert ("FGF20", "FGF16") or ("FGF16", "FGF20") in G.edges
        assert ("FGF17", "FGF20") in G.edges
        assert ("FGFR3", "FGF20") in G.edges
        assert not ("FGF23", "HS") in G.edges
        assert not ("FGF23", "FGF23") in G.edges
        assert ("FGF23", "FRS2") in G.edges
        assert not ("HS", "FGF17") in G.edges

    def test_connects_components_of_same_complex(self):
        G = TestGeneNetworkNoSmallMoleculesClass.G_glycolysis

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
