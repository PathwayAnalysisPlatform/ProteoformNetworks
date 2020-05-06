from Python.interaction_network import create_graph


# Test graph creation for proteins
class TestProteoformNetworkClass:
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    G_glycolysis = create_graph("R-HSA-9634600", level="proteoforms")

    # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
    G_traffic = create_graph("R-HSA-5627083", level="proteoforms")

    # Pathway "FRS-mediated FGFR3 signaling" R-HSA-5654706
    G_signal = create_graph("R-HSA-5654706", level="proteoforms")

    def test_create_graph_num_edges(self):
        G = TestProteoformNetworkClass.G_glycolysis
        print(G.edges)
        assert len(G.edges) == 52

    def test_create_graph_num_vertices(self):
        G = TestProteoformNetworkClass.G_glycolysis
        print(G.nodes)
        assert len(G.nodes) == 19

    # Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
    def test_connects_inputs_with_outputs(self):
        G = TestProteoformNetworkClass.G_glycolysis
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600

        # Input to output interactions for reaction R-HSA-163773:
        assert ("P16118", "ADP") in G.edges
        assert ("P16118", "P16118;00046:33") in G.edges
        assert ("ATP", "ADP") in G.edges
        assert not ("ATP", "P16118") in G.edges
        assert ("ATP", "P16118;00046:33") in G.edges

        # Input to output interactions for reaction R-HSA-163750
        assert not ("P16118", "Pi") in G.edges
        assert ("P16118;00046:33", "Pi") in G.edges
        assert ("P16118;00046:33", "P16118") in G.edges
        assert ("H2O", "Pi") in G.edges
        assert ("H2O", "P16118") in G.edges

        # Input to output interaction for reaction R-HSA-71802
        assert ("Fru(6)P", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("Fru(6)P", "ADP") in G.edges
        assert ("ATP", "D-Fructose 2,6-bisphosphate") in G.edges
        assert ("ATP", "ADP") in G.edges

    def test_does_not_connect_input_to_output_when_same_molecule(self):
        # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
        G = TestProteoformNetworkClass.G_traffic

        # In reaction R-HSA-5627071 there are same set of proteoforms as input and output, since they are in
        # different sub cellular locations. Nevertheless, there are no self interactions in our network
        assert not ("P13569", "P13569") in G.edges
        assert not ("P17081", "P17081") in G.edges
        assert not ("Q9HD26", "Q9HD26") in G.edges
        assert not ("GTP", "GTP") in G.edges

    def test_connects_catalysts_with_outputs(self):
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
        G = TestProteoformNetworkClass.G_glycolysis

        # Catalysts to output interactions for reaction R-HSA-163773:
        assert not ("P17612", "P16118") in G.edges
        assert ("P17612", "P16118;00046:33") in G.edges
        assert ("P17612", "ADP") in G.edges
        assert not ("P22694", "P16118") in G.edges
        assert ("P22694", "P16118;00046:33") in G.edges
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
        assert not ("P16118", "Fru(6)P") in G.edges
        assert ("P16118;00046:33", "Fru(6)P") in G.edges
        assert not ("P16118", "Pi") in G.edges
        assert ("P16118;00046:33", "Pi") in G.edges

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
        G = TestProteoformNetworkClass.G_traffic

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

    def test_connects_regulators_with_outputs_2(self):
        # Pathway "FRS-mediated FGFR3 signaling" R-HSA-5654706
        G = TestProteoformNetworkClass.G_signal

        # In reaction R-HSA-5654408
        assert ("Q9NP95", "O43320") or ("O43320", "Q9NP95") in G.edges
        assert ("O60258-1", "Q9NP95") in G.edges
        assert ("P22607-1;00048:577,00048:647,00048:648,00048:724,00048:760,00048:770", "Q9NP95") in G.edges
        assert ("Q9GZV9;00164:178", "HS") in G.edges
        assert not ("Q9GZV9;00164:178", "Q9GZV9;00164:178") in G.edges
        assert ("Q9GZV9;00164:178", "Q8WU20;00048:196,00048:306,00048:349,00048:392,00048:436,00048:471") in G.edges
        assert ("HS", "O60258-1") in G.edges

    def test_connects_components_of_same_complex(self):
        G = TestProteoformNetworkClass.G_glycolysis

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

    def test_connect_components_of_same_complex_2(self):
        G = TestProteoformNetworkClass.G_signal

        assert ("P01116-1;00115:179,01116:186", "GTP") in G.edges
        assert ("P01111;00115:181,01116:186", "GTP") in G.edges
        assert ("P01112;00115:181,00115:184,01116:186", "GTP") in G.edges
        assert ("P01116-2;01116:186", "GTP") in G.edges

        assert ("P01112;00115:181,00115:184,01116:186", "P01111;00115:181,01116:186") in G.edges
        assert ("P01112;00115:181,00115:184,01116:186", "P01116-2;01116:186") in G.edges

        assert ("P01116-2;01116:186", "P01111;00115:181,01116:186") in G.edges
        assert ("P01116-2;01116:186", "P01112;00115:181,00115:184,01116:186") in G.edges

        assert ("P01111;00115:181,01116:186", "P01116-2;01116:186") in G.edges
        assert ("P01111;00115:181,01116:186", "P01112;00115:181,00115:184,01116:186") in G.edges

        assert ("Q07889", "P62993-1") in G.edges

        assert ("HS", "P22607-2;00048:579,00048:649,00048:650,00048:726,00048:762,00048:772") in G.edges
        assert ("P22607-2;00048:579,00048:649,00048:650,00048:726,00048:762,00048:772", "P05230") in G.edges

        assert ("P55075-1", "P09038") in G.edges
        assert ("O60258-1", "Q9GZV9;00164:178") in G.edges

