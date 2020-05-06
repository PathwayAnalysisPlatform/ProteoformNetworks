import pandas as pd
import pytest

from Python.network_topology_queries import get_reaction_participants_by_pathway, get_pathways, \
    get_reactions_by_pathway


@pytest.fixture()
def setup_participants():
    print("Setting up glycolysis participants")
    return get_reaction_participants_by_pathway("R-HSA-70171")


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
    df = get_reaction_participants_by_pathway("R-HSA-1237044")
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
    df = get_reaction_participants_by_pathway("R-HSA-70171")
    # DefinedSet R-HSA-450097 has 4 members
    assert (df['Entity'] == "R-HSA-450094").any()
    assert (df['Entity'] == "R-HSA-70395").any()
    assert (df['Entity'] == "R-HSA-70378").any()
    assert (df['Entity'] == "R-HSA-70412").any()


# Test: If small molecules disabled, cypher query returns only the gene direct participants
def test_query_for_pathway_participants_disable_small_molecules():
    df = get_reaction_participants_by_pathway("R-HSA-70171", showSmallMolecules=False)
    assert len(df) == 31  # Only EWAS
    # Has EWAS participants
    assert (df['Entity'] == "R-HSA-5696062").any()
    assert (df['Entity'] == "R-HSA-6798334").any()
    assert (df['Entity'] == "R-HSA-6799597").any()
    # Does not have Small molecule participants
    assert (df['Entity'] != "R-ALL-449701").all()
    assert (df['Entity'] != "R-ALL-70113").all()
    assert (df['Entity'] != "R-ALL-217314").all()


# Test: If small molecules disabled, cypher query returns only the gene complex decomposed participants
def test_query_for_pathway_participants_complexes_show_only_ewas():
    df = get_reaction_participants_by_pathway("R-HSA-70171", showSmallMolecules=False)
    # The pathway has the reaction R-HSA-5696021 which has the complex R-HSA-5696043 as participant
    # THe components of the complex are 2, one EWAS and one small molecule:
    assert (df['Entity'] == 'R-HSA-5696062').any()  # EWAS ADPGK
    assert (df['Entity'] != 'R-ALL-217314').all()  # Small molecule Mg2+


# Test: Get pathway participants as genes
def test_query_for_pathway_participants_as_genes():
    df = get_reaction_participants_by_pathway("R-HSA-70171", level="genes")
    assert len(df) == 102

    # Participant protein Q9BRR6 should be in the result as a gene name: ADPGK not as UniProt accession
    assert (df['Id'] == 'ADPGK').any()
    assert (df['Id'] != 'Q9BRR6').all()

    # Simmilar for the next proteins
    assert (df['Id'] == 'BPGM').any()
    assert (df['Id'] != 'P07738').all()

    assert (df['Id'] == 'BPGM').any()
    assert (df['Id'] != 'P07738').all()

    assert (df['Id'] == 'PKLR').any()
    assert (df['Id'] != 'P07738').all()


def test_query_for_pathway_participants_as_genes_trims_gene_id():
    df = get_reaction_participants_by_pathway("R-HSA-70171", level="genes")
    assert ((df['Entity'] == 'R-HSA-70097') & (df['Id'] == "PKLR")).any()
    assert ((df['Entity'] == 'R-HSA-450658') & (df['Id'] == "PKM")).any()
    assert not ((df['Entity'] == 'R-HSA-450658') & (df['Id'] == "PKM-2 [cytosol]")).any()
    assert not ((df['Entity'] == 'R-HSA-450658') & (df['Id'] == "P62993")).any()
    print(df.loc[(df['Entity'] == 'R-HSA-70097') & (df['Id'] == "PKLR")])
    assert not ((df['Entity'] == 'R-HSA-70097') & (df['Id'] == "PKLR-1 [cytosol]")).any()
    assert not ((df['Entity'] == 'R-HSA-211388') & (df['Id'] == "PKLR-2")).any()


def test_query_for_pathway_participants_replaces_small_molecule_names():
    df = get_reaction_participants_by_pathway("R-HSA-70171", level="genes")
    assert not ((df['Entity'] == 'R-ALL-29370') & (df['Id'] == '456216')).any()
    assert ((df['Entity'] == 'R-ALL-29370') & (df['Id'] == 'ADP')).any()
    assert not ((df['Entity'] == 'R-ALL-29926') & (df['Id'] == '18420')).any()
    assert not ((df['Entity'] == 'R-ALL-29926') & (df['Id'] == 'Mg2+ [cytosol]')).any()
    assert ((df['Entity'] == 'R-ALL-29926') & (df['Id'] == 'Mg2+')).any()


# Test: Get pathway participants as proteins
def test_query_for_pathway_participants_as_proteins(setup_participants):
    df = get_reaction_participants_by_pathway("R-HSA-70171", level="proteins")
    assert not ((df['Entity'] == 'R-HSA-5696062') & (df['Id'] == 'ADPGK')).any()
    assert ((df['Entity'] == 'R-HSA-5696062') & (df['Id'] == 'Q9BRR6')).any()

    assert not ((df['Entity'] == "R-HSA-6798334") & (df['Id'] == 'BPGM')).any()
    assert ((df['Entity'] == "R-HSA-6798334") & (df['Id'] == 'P07738')).any()

    assert not ((df['Entity'] == 'R-HSA-70412') & (df['Id'] == 'HK3')).any()
    assert ((df['Entity'] == 'R-HSA-70412') & (df['Id'] == 'P52790')).any()


def test_query_for_pathway_participants_as_proteins_implicit_parameter(setup_participants):
    df = setup_participants
    assert not ((df['Entity'] == 'R-HSA-5696062') & (df['Id'] == 'ADPGK')).any()
    assert ((df['Entity'] == 'R-HSA-5696062') & (df['Id'] == 'Q9BRR6')).any()

    assert not ((df['Entity'] == "R-HSA-6798334") & (df['Id'] == 'BPGM')).any()
    assert ((df['Entity'] == "R-HSA-6798334") & (df['Id'] == 'P07738')).any()

    assert not ((df['Entity'] == 'R-HSA-70412') & (df['Id'] == 'HK3')).any()
    assert ((df['Entity'] == 'R-HSA-70412') & (df['Id'] == 'P52790')).any()


def test_query_for_pathway_participants_as_proteins_complex_should_not_be_in_records(setup_participants):
    df = setup_participants
    assert not (df["Entity"] == "R-HSA-5696043").any()
    assert not (df["Name"] == "BPGM dimer [cytosol]").any()
    assert not (df["Entity"] == "R-HSA-6799598").any()
