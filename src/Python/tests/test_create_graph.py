import pandas as pd
import pytest

from Python.interaction_network import create_graph
from Python.network_topology_queries import get_reactions_and_participants_by_pathway, get_pathways, \
    get_reactions_by_pathway


@pytest.fixture()
def setup_participants():
    print("Setting up glycolysis participants")
    return get_reactions_and_participants_by_pathway("R-HSA-70171")


# Test: Cypher query to get all pathways gets all of them correctly
def test_get_pathways():
    pathways = get_pathways()
    assert len(pathways) == 2362
    assert type(pathways) == pd.DataFrame
    assert ((pathways['stId'] == 'R-HSA-9612973') & (pathways['displayName'] == 'Autophagy')).any()
    assert ((pathways['stId'] == 'R-HSA-1640170') & (pathways['displayName'] == 'Cell Cycle')).any()
    assert ((pathways['stId'] == 'R-HSA-70171') & (pathways['displayName'] == 'Glycolysis')).any()


# Test: Query to get reactions of pathway returns the correct list of reactions
def test_get_reactions_of_pathway():
    reactions = get_reactions_by_pathway("R-HSA-70171")
    assert len(reactions) == 15
    assert (reactions['reaction'] == "R-HSA-8955794").any()
    assert (reactions['reaction'] == "R-HSA-6799604").any()
    assert (reactions['reaction'] == "R-HSA-70467").any()


# Test: Query to get participants of a pathway gets all participant reactions.
def test_query_for_pathway_participants_has_all_reactions(setup_participants):
    df = setup_participants
    assert len(df['Reaction'].unique()) == 15
    assert (df['Reaction'] == "R-HSA-8955794").any()
    assert (df['Reaction'] == "R-HSA-6799604").any()
    assert (df['Reaction'] == "R-HSA-70467").any()


# Test: Query to get participants of a pathway returns all small molecule participants
def test_query_for_pathway_participans_returns_all_simple_molecules(setup_participants):
    df = setup_participants
    print(df.loc[df['Type'] == 'SimpleEntity'])
    assert len(df.loc[df['Type'] == 'SimpleEntity']['Entity'].unique()) == 32


# Test: Query to get participants of a pathway returns all gene participants
def test_query_for_pathway_participans_returns_all_ewas(setup_participants):
    df = setup_participants
    print(df.loc[df['Type'] == 'EntityWithAccessionedSequence'])
    assert len(df.loc[df['Type'] == 'EntityWithAccessionedSequence']['Entity'].unique()) == 31


# Test: Query for participants of a pathway returns all the gene and small molecule participants decomposing complexes
# Pathway R-HSA-1237044 hast complex participants in its reactions, some complexes also have complex components.
# Checks if the component molecules of the complex components are in the result
def test_query_for_pathway_participants_decomposes_complexes():
    df = get_reactions_and_participants_by_pathway("R-HSA-1237044")
    # Reaction R-HSA-1237325 in the pathway has complex participants: "R-HSA-1237320"
    # Complex R-HSA-1237320 has 6 participant molecules:
    assert (df['Entity'] == "R-ALL-71185").any()
    assert (df['Entity'] == "R-HSA-1008268").any()
    assert (df['Entity'] == "R-ALL-29368").any()
    assert (df['Entity'] == "R-ALL-71185").any()
    assert (df['Entity'] == "R-ALL-29368").any()
    assert (df['Entity'] == "R-HSA-1008196").any()


# Test: Query for participants of a pathway returns all the gene and small molecule participants decomposing sets
# Pathway R-HSA-70171 has reactions with EntitySets as participants, like reaction R-HSA-70420 with set R-HSA-450097
def test_query_for_pathway_participants_decomposes_sets():
    df = get_reactions_and_participants_by_pathway("R-HSA-70171")
    # DefinedSet R-HSA-450097 has 4 members
    assert (df['Entity'] == "R-HSA-450094").any()
    assert (df['Entity'] == "R-HSA-70395").any()
    assert (df['Entity'] == "R-HSA-70378").any()
    assert (df['Entity'] == "R-HSA-70412").any()

# Test: If small molecules disabled, cypher query returns only the gene direct participants
# Test: If small molecules disabled, cypher query returns only the gene complex decomposed participants

# Graph creation

def test_create_graph_num_edges(setup_participants):
    df = setup_participants
    G = create_graph(df)
    assert len(G.edges) == 107


def test_create_graph_num_vertices(setup_participants):
    df = setup_participants
    G = create_graph(df)
    assert len(G.nodes) == 61


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
