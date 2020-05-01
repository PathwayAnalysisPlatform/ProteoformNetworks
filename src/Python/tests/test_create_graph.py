from pytest import fail
import pandas

from Python.interaction_network import create_graph
from Python.network_topology_queries import get_reactions_and_participants_by_pathway, get_pathways, \
    get_reactions_by_pathway


def test_create_graph_num_edges():
    df = get_reactions_and_participants_by_pathway("R-HSA-70171")
    G = create_graph(df)
    assert len(G.edges) == 107


def test_create_graph_num_vertices():
    df = get_reactions_and_participants_by_pathway("R-HSA-70171")
    G = create_graph(df)
    assert len(G.nodes) == 61


# For gene network:
class TestQueriesClass:
    # Test: Cypher query to get all pathways gets all of them correctly
    def test_get_pathways(self):
        pathways = get_pathways()
        assert len(pathways) == 2362
        assert type(pathways) == pandas.DataFrame
        assert ((pathways['stId'] == 'R-HSA-9612973') & (pathways['displayName'] == 'Autophagy')).any()
        assert ((pathways['stId'] == 'R-HSA-1640170') & (pathways['displayName'] == 'Cell Cycle')).any()
        assert ((pathways['stId'] == 'R-HSA-70171') & (pathways['displayName'] == 'Glycolysis')).any()


    # Test: Query to get participants of a pathway gets all participant reactions.
    def test_get_reactions_of_pathway(self):
        reactions = get_reactions_by_pathway("R-HSA-70171")
        assert len(reactions) == 15
        assert (reactions['reaction'] == "R-HSA-8955794").any()
        assert (reactions['reaction'] == "R-HSA-6799604").any()
        assert (reactions['reaction'] == "R-HSA-70467").any()


    def test_query_for_pathway_participants_has_all_reactions(self):
        df = get_reactions_and_participants_by_pathway("R-HSA-70171")
        assert len(df['Reaction'].unique()) == 15
        assert (df['Reaction'] == "R-HSA-8955794").any()
        assert (df['Reaction'] == "R-HSA-6799604").any()
        assert (df['Reaction'] == "R-HSA-70467").any()

# Test: Query to get participants of a pathway returns all gene and small molecule participants
# Test: Query to get participants of a pathway returns all the gene and small molecule participants decomposing complexes
# Test: Query to get participants of a pathway returns all the gene and small molecule participants decomposing sets
# Test: If small molecules disabled, cypher query returns only the gene direct participants
# Test: If small molecules disabled, cypher query returns only the gene complex decomposed participants
# Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
# Test: Receive a Reaction with a direct participant EWAS input and EWAS output --> connects them
# Test: Receive a Reaction with a direct participant EWAS catalyst and EWAS output --> connects them
# Test: Receive a Reaction with a direct participant EWAS regulator and EWAS output --> connects them
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

# For protein network:
# All the tests of gene network

# For proteoform networks

# Run here
