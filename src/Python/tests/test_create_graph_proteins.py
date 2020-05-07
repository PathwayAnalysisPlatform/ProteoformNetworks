from interaction_network import create_graph


# Test graph creation for proteins
class TestProteinNetworkClass:
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    G_glycolysis = create_graph("R-HSA-9634600", level="proteins")

    # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
    G_traffic = create_graph("R-HSA-5627083", level="proteins")

    def test_create_graph_num_edges(self):
        G = TestProteinNetworkClass.G_glycolysis
        print(G.edges)
        assert len(G.edges) == 46

    def test_create_graph_num_vertices(self):
        G = TestProteinNetworkClass.G_glycolysis
        print(G.nodes)
        assert len(G.nodes) == 18

    # Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
    def test_connects_inputs_with_outputs(self):
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
        G = TestProteinNetworkClass.G_glycolysis
        # Input to output interactions for reaction R-HSA-163773:
        assert ("P16118", "ADP") in G.edges
        assert ("ATP", "ADP") in G.edges
        assert ("ATP", "P16118") in G.edges

        # Input to output interactions for reaction R-HSA-163750
        assert ("P16118", "Pi") in G.edges
        assert ("H2O", "Pi") in G.edges
        assert ("H2O", "P16118") in G.edges

        # Input to output interaction for reaction R-HSA-71802
        assert ("Fru(6)P", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("Fru(6)P", "ADP") in G.edges
        assert ("ATP", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("ATP", "ADP") in G.edges

    def test_connects_input_genes_with_small_outputs_not_when_is_same_molecule(self):
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
        G = TestProteinNetworkClass.G_glycolysis
        assert not ("P16118", "P16118") in G.edges

    def test_connects_catalysts_with_outputs(self):
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
        G = TestProteinNetworkClass.G_glycolysis

        # Catalysts to output interactions for reaction R-HSA-163773:
        assert ("P17612", "P16118") in G.edges
        assert ("P17612", "ADP") in G.edges
        assert ("P22694", "P16118") in G.edges
        assert ("P22694", "ADP") in G.edges

        # Catalysts to output interactions for reaction R-HSA-163750
        assert ("P30153", "P16118") in G.edges
        assert ("P30153", "Pi") in G.edges
        assert ("P30154", "P16118") in G.edges
        assert ("P30154", "Pi") in G.edges
        assert ("P67775", "P16118") in G.edges
        assert ("P67775", "Pi") in G.edges
        assert ("P62714", "P16118") in G.edges
        assert ("P62714", "Pi") in G.edges
        assert ("Q14738", "P16118") in G.edges
        assert ("Q14738", "Pi") in G.edges

        # Catalysts to output interactions for reaction R-HSA-70262
        assert ("P16118", "Fru(6)P") in G.edges
        assert ("P16118", "Pi") in G.edges

        # Catalysts to output interactions for reaction R-HSA-71802
        assert ("Q16875", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("Q16875", "ADP") in G.edges
        assert ("Q16877", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("Q16877", "ADP") in G.edges
        assert ("P16118", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("P16118", "ADP") in G.edges
        assert ("O60825", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("O60825", "ADP") in G.edges

    def test_connects_regulators_with_outputs(self):
        # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
        G = TestProteinNetworkClass.G_traffic

        # Regulator to output interactions for reaction R-HSA-5627071:
        assert ("P17081", "P13569") in G.edges
        assert ("Q9HD26", "P13569") in G.edges
        assert not ("P13569", "P13569") in G.edges
        assert ("GTP", "P13569") in G.edges

        assert ("P17081", "Q9HD26") in G.edges
        assert not ("Q9HD26", "Q9HD26") in G.edges
        assert ("P13569", "Q9HD26") in G.edges
        assert ("GTP", "Q9HD26") in G.edges

        assert not ("P17081", "P17081") in G.edges
        assert ("Q9HD26", "P17081") in G.edges
        assert ("P13569", "P17081") in G.edges
        assert ("GTP", "P17081") in G.edges

        assert ("P17081", "GTP") in G.edges
        assert ("Q9HD26", "GTP") in G.edges
        assert ("P13569", "GTP") in G.edges
        assert not ("GTP", "GTP") in G.edges

    def test_connects_components_of_same_complex(self):
        G = TestProteinNetworkClass.G_glycolysis

        # As part of PP2A-ABdeltaC complex R-HSA-165961
        assert not ("P30153", "P30153") in G.edges
        assert ("P30153", "P30154") in G.edges
        assert ("P30154", "P30153") in G.edges
        assert ("P30153", "P67775") in G.edges
        assert ("P30153", "P62714") in G.edges
        assert ("P30153", "Q14738") in G.edges
        assert ("Q14738", "P30154") in G.edges
        assert ("Q14738", "P30153") in G.edges
        assert ("Q14738", "P62714") in G.edges
