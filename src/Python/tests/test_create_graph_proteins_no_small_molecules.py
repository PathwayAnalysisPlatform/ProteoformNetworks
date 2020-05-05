import pytest

from Python.interaction_network import create_graph
from Python.network_topology_queries import get_reactions_and_participants_by_pathway


# Test graph creation with genes and small molecules disabled
class TestGeneNetworkNoSmallMoleculesClass:

    @pytest.fixture()
    def setup_gene_participants(self, scope='class'):
        print("\nSetup Class")
        df = get_reactions_and_participants_by_pathway("R-HSA-9634600", level="proteins", showSmallMolecules=False)
        return create_graph(df)

    def test_create_graph_num_edges(self, setup_gene_participants):
        G = setup_gene_participants
        print(G.edges)
        assert len(G.edges) == 8

    def test_create_graph_num_vertices(self, setup_gene_participants):
        G = setup_gene_participants
        print(G.nodes)
        assert len(G.nodes) == 12

    # Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
    def test_connects_inputs_with_outputs(self, setup_gene_participants):
        G = setup_gene_participants
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600

        # Input to output interactions for reaction R-HSA-163750
        assert not ("P16118", "Pi") in G.edges
        assert not ("P16118", "P16118") in G.edges
        assert not ("H2O", "Pi") in G.edges
        assert not ("H2O", "P16118") in G.edges

        # Input to output interactions for reaction R-HSA-163773:
        assert not ("P16118", "ADP") in G.edges
        assert not ("ATP", "ADP") in G.edges
        assert not ("ATP", "P16118") in G.edges

        # Input to output interaction for reaction R-HSA-71802
        assert not ("Fru(6)P", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("Fru(6)P", "ADP") in G.edges
        assert not ("ATP", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("ATP", "ADP") in G.edges

    def test_connects_catalysts_with_outputs(self, setup_gene_participants):
        G = setup_gene_participants
        # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600

        # Catalysts to output interactions for reaction R-HSA-163773:
        assert ("P17612", "P16118") in G.edges
        assert ("P22612", "P16118") in G.edges
        assert ("P22694", "P16118") in G.edges
        assert not ("P17612", "ADP") in G.edges
        assert not ("P22694", "ADP") in G.edges

        # Catalysts to output interactions for reaction R-HSA-163750
        assert ("P30153", "P16118") in G.edges
        assert not ("P30153", "Pi") in G.edges
        assert ("P30154", "P16118") in G.edges
        assert not ("P30154", "Pi") in G.edges
        assert ("P67775", "P16118") in G.edges
        assert not ("P67775", "Pi") in G.edges
        assert ("P62714", "P16118") in G.edges
        assert not ("P62714", "Pi") in G.edges
        assert ("Q14738", "P16118") in G.edges
        assert not ("Q14738", "Pi") in G.edges

        # Catalysts to output interactions for reaction R-HSA-70262
        assert not ("P16118", "Fru(6)P") in G.edges
        assert not ("P16118", "Pi") in G.edges

        # Catalysts to output interactions for reaction R-HSA-71802
        assert not ("Q16875", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("Q16875", "ADP") in G.edges
        assert not ("Q16877", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("Q16877", "ADP") in G.edges
        assert not ("P16118", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("P16118", "ADP") in G.edges
        assert not ("O60825", "D-Fructose 2,6-bisphosphate") in G.edges
        assert not ("O60825", "ADP") in G.edges

    def test_connects_regulators_with_outputs(self):
        df = get_reactions_and_participants_by_pathway("R-HSA-5627083", level="proteins", showSmallMolecules=False)
        G = create_graph(df)
        # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083

        # Regulator to output interactions for reaction R-HSA-5627071:
        assert ("P17081", "P13569") in G.edges
        assert ("Q9HD26", "P13569") in G.edges
        assert not ("P13569", "P13569") in G.edges
        assert not ("GTP", "P13569") in G.edges

        assert ("P17081", "Q9HD26") in G.edges
        assert not ("Q9HD26", "Q9HD26") in G.edges
        assert ("P13569", "Q9HD26") in G.edges
        assert not ("GTP", "Q9HD26") in G.edges

        assert not ("P17081", "P17081") in G.edges
        assert ("Q9HD26", "P17081") in G.edges
        assert ("P13569", "P17081") in G.edges
        assert not ("GTP", "P17081") in G.edges

        assert not ("P17081", "GTP") in G.edges
        assert not ("Q9HD26", "GTP") in G.edges
        assert not ("P13569", "GTP") in G.edges
        assert not ("GTP", "GTP") in G.edges
