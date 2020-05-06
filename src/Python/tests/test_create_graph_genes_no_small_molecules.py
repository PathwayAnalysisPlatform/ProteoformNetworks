import pytest

from Python.interaction_network import create_graph
from Python.network_topology_queries import get_reactions_and_participants_by_pathway


# Test graph creation with genes and small molecules disabled
class TestGeneNetworkNoSmallMoleculesClass:

    @pytest.fixture()
    def setup_gene_participants(self, scope='class'):
        print("\nSetup Class")
        df = get_reactions_and_participants_by_pathway("R-HSA-9634600", level="genes", showSmallMolecules=False)
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

    def test_connects_catalysts_with_outputs(self, setup_gene_participants):
        G = setup_gene_participants
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
        df = get_reactions_and_participants_by_pathway("R-HSA-5627083", level="genes", showSmallMolecules=False)
        G = create_graph(df)
        # Pathway "RHO GTPases regulate CFTR trafficking" R-HSA-5627083

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

# Test graph creation with proteins and small molecules disabled

# Test graph creation for proteins

# Test graph creation for proteoforms

# Test: Receive a Reaction with a direct participant SimpleEntity input and SimpleEntity output --> connects them
# Test: Receive a Reaction with a direct participant SimpleEntity catalyst and SimpleEntity output --> connects them
# Test: Receive a Reaction with a direct participant SimpleEntity regulator and SimpleEntity output --> connects them
# Test: Receive a Reaction with a complex input and EWAS output --> connects each complex component to the output
# Test: Receive a Reaction with a complex input --> does not connect complex components among them
# Test: Receive a Reaction with a complex made of complexes as input --> connects each component to the output
# Test: Receive a Reaction with a set as input --> connects each member to the output
# Test: Receive a Reaction with a set as input --> does not connect each member with each other
# Test: Receive a Reaction with a complex made of sets as input --> connects each member with each other
# Test: If small molecules disabled then connect only gene participants from input to output
# Test: If small molecules disabled then connect only gene participants from catalyst to output
# Test: If small molecules disabled then connect only gene participants from regulator to output
# Test: Given a pathway, the graph merges participants of all reactions in the pathway

# Test: Participant nodes have the attribute "Type" correctly as "Gene", "Protein", "SimpleEntity", etc.
# Test: Get pathway participants as proteoforms

# For protein network:
# All the tests of gene network

# For proteoform networks

# Run here