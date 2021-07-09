import pytest

from config import genes, with_sm, no_sm, with_unique_sm
from lib.networks import create_pathway_interaction_network


@pytest.fixture(scope="session")
def glycolysis_graph_with_sm(tmpdir_factory):
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-9634600", genes, with_sm, graphs_path)

@pytest.fixture(scope="session")
def glycolysis_graph_with_unique_sm(tmpdir_factory):
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-9634600", genes, with_unique_sm, graphs_path)


@pytest.fixture(scope="session")
def traffic_graph_with_sm(tmpdir_factory):
    # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-5627083", genes, with_sm, graphs_path)


@pytest.fixture(scope="session")
def traffic_graph_with_unique_sm(tmpdir_factory):
    # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-5627083", genes, with_unique_sm, graphs_path)


@pytest.fixture(scope="session")
def pathway_free_fatty_acid_receptors(tmpdir_factory):
    # Pathway "Free fatty acid receptors" R-HSA-444209
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-444209", genes, with_sm, graphs_path)


def test_create_graph_wrong_level_raises_exception():
    with pytest.raises(Exception):
        create_pathway_interaction_network("R-HSA-9634600", "banana", no_sm)


def test_create_graph_num_edges(glycolysis_graph_with_sm):
    g = glycolysis_graph_with_sm
    print(g.edges)
    assert len(g.edges) == 46



def test_connects_input_genes_with_small_outputs_not_when_is_same_molecule(glycolysis_graph_with_sm):
    g = glycolysis_graph_with_sm
    assert not ("PFKFB1", "PFKFB1") in g.edges

    # Catalysts to output interactions for reaction R-HSA-163773:
    assert ("PRKACA", "PFKFB1") in g.edges
    assert ("PRKACA", "sm_ADP") in g.edges
    assert ("PRKACB", "PFKFB1") in g.edges
    assert ("PRKACB", "sm_ADP") in g.edges

    # Catalysts to output interactions for reaction R-HSA-163750
    assert ("PPP2R1A", "PFKFB1") in g.edges
    assert ("PPP2R1A", "sm_Pi") in g.edges
    assert ("PPP2R1B", "PFKFB1") in g.edges
    assert ("PPP2R1B", "sm_Pi") in g.edges
    assert ("PPP2CA", "PFKFB1") in g.edges
    assert ("PPP2CA", "sm_Pi") in g.edges
    assert ("PPP2CB", "PFKFB1") in g.edges
    assert ("PPP2CB", "sm_Pi") in g.edges
    assert ("PPP2R5D", "PFKFB1") in g.edges
    assert ("PPP2R5D", "sm_Pi") in g.edges

    # Catalysts to output interactions for reaction R-HSA-70262
    assert ("PFKFB1", "sm_Fru(6)P") in g.edges
    assert ("PFKFB1", "sm_Pi") in g.edges

    # Catalysts to output interactions for reaction R-HSA-71802
    assert ("PFKFB3", "sm_D-Fructose_2,6-bisphosphate") in g.edges
    assert ("PFKFB3", "sm_ADP") in g.edges
    assert ("PFKFB4", "sm_D-Fructose_2,6-bisphosphate") in g.edges
    assert ("PFKFB4", "sm_ADP") in g.edges
    assert ("PFKFB1", "sm_D-Fructose_2,6-bisphosphate") in g.edges
    assert ("PFKFB1", "sm_ADP") in g.edges
    assert ("PFKFB2", "sm_D-Fructose_2,6-bisphosphate") in g.edges
    assert ("PFKFB2", "sm_ADP") in g.edges


def test_connects_regulators_with_outputs(traffic_graph_with_sm):
    tg = traffic_graph_with_sm

    # Regulator to output interactions for reaction R-HSA-5627071:
    assert ("RHOQ", "CFTR") in tg.edges
    assert ("GOPC", "CFTR") in tg.edges
    assert not ("CFTR", "CFTR") in tg.edges
    assert ("sm_GTP", "CFTR") in tg.edges

    assert ("RHOQ", "GOPC") in tg.edges
    assert not ("GOPC", "GOPC") in tg.edges
    assert ("CFTR", "GOPC") in tg.edges
    assert ("sm_GTP", "GOPC") in tg.edges

    assert not ("RHOQ", "RHOQ") in tg.edges
    assert ("GOPC", "RHOQ") in tg.edges
    assert ("CFTR", "RHOQ") in tg.edges
    assert ("sm_GTP", "RHOQ") in tg.edges

    assert ("RHOQ", "sm_GTP") in tg.edges
    assert ("GOPC", "sm_GTP") in tg.edges
    assert ("CFTR", "sm_GTP") in tg.edges
    assert not ("sm_GTP", "sm_GTP") in tg.edges


def test_connects_components_of_same_complex(glycolysis_graph_with_sm):

    g = glycolysis_graph_with_sm
    # As part of PP2A-ABdeltaC complex R-HSA-165961
    assert not ("PPP2R1A", "PPP2R1A") in g.edges
    assert ("PPP2R1A", "PPP2R1B") in g.edges
    assert ("PPP2R1B", "PPP2R1A") in g.edges
    assert ("PPP2R1A", "PPP2CA") in g.edges
    assert ("PPP2R1A", "PPP2CB") in g.edges
    assert ("PPP2R1A", "PPP2R5D") in g.edges
    assert ("PPP2R5D", "PPP2R1B") in g.edges
    assert ("PPP2R5D", "PPP2R1A") in g.edges
    assert ("PPP2R5D", "PPP2CB") in g.edges

# Test it connectects correctly the reactions in pathway: R-HSA-444209
# Some reactions have the same gene as input and output

def test_connectes_participants_of_reaction_R_HSA_444191(pathway_free_fatty_acid_receptors):
    g = pathway_free_fatty_acid_receptors
    assert g.has_node("FFAR4")
    assert ("FFAR4", "FFAR4")
    assert ("sm_DTTA", "FFAR4")
    assert ("sm_OLEA", "FFAR4")

def test_connectes_participants_of_reaction_R_HSA_444191(pathway_free_fatty_acid_receptors):
    g = pathway_free_fatty_acid_receptors
    assert g.has_node("sm_12(S)-HETE")
    assert g.has_node("GPR31")
    assert ("sm_12(S)-HETE", "sm_12(S)-HETE")
    assert ("sm_12(S)-HETE", "GPR31")
    assert ("GPR31", "sm_12(S)-HETE")
