import pytest

from config import no_sm, genes
from lib.networks import create_pathway_interaction_network


@pytest.fixture(scope="session")
def glycolysis_graph(tmpdir_factory):
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-9634600", genes, no_sm, graphs_path)


@pytest.fixture(scope="session")
def traffic_graph(tmpdir_factory):
    # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-5627083", genes, no_sm, graphs_path)

@pytest.fixture(scope="session")
def signal_graph(tmpdir_factory):
    # Pathway "FRS-mediated FGFR3 signaling"
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-5654706", genes, no_sm, graphs_path)

def test_create_graph_num_edges(glycolysis_graph):
    print(glycolysis_graph.edges)
    assert len(glycolysis_graph.edges) == 18

def test_create_graph_num_vertices(glycolysis_graph):
    print(glycolysis_graph.nodes)
    assert len(glycolysis_graph.nodes) == 12

# Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
def test_connects_inputs_with_outputs(glycolysis_graph):

    # Input to output interactions for reaction R-HSA-163750
    assert not ("PFKFB1", "sm_Pi") in glycolysis_graph.edges
    assert not ("PFKFB1", "PFKFB1") in glycolysis_graph.edges
    assert not ("sm_H2O", "sm_Pi") in glycolysis_graph.edges
    assert not ("sm_H2O", "PFKFB1") in glycolysis_graph.edges

    # Input to output interactions for reaction R-HSA-163773:
    assert not ("PFKFB1", "sm_ADP") in glycolysis_graph.edges
    assert not ("sm_ATP", "sm_ADP") in glycolysis_graph.edges
    assert not ("sm_ATP", "PFKFB1") in glycolysis_graph.edges

    # Input to output interaction for reaction R-HSA-71802
    assert not ("sm_Fru(6)P", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("sm_Fru(6)P", "sm_ADP") in glycolysis_graph.edges
    assert not ("sm_ATP", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("sm_ATP", "sm_ADP") in glycolysis_graph.edges

def test_connects_catalysts_with_outputs(glycolysis_graph):

    # Catalysts to output interactions for reaction R-HSA-163773:
    assert ("PRKACA", "PFKFB1") in glycolysis_graph.edges
    assert ("PRKACG", "PFKFB1") in glycolysis_graph.edges
    assert ("PRKACB", "PFKFB1") in glycolysis_graph.edges
    assert not ("PRKACA", "sm_ADP") in glycolysis_graph.edges
    assert not ("PRKACB", "sm_ADP") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-163750
    assert ("PPP2R1A", "PFKFB1") in glycolysis_graph.edges
    assert not ("PPP2R1A", "sm_Pi") in glycolysis_graph.edges
    assert ("PPP2R1B", "PFKFB1") in glycolysis_graph.edges
    assert not ("PPP2R1B", "sm_Pi") in glycolysis_graph.edges
    assert ("PPP2CA", "PFKFB1") in glycolysis_graph.edges
    assert not ("PPP2CA", "sm_Pi") in glycolysis_graph.edges
    assert ("PPP2CB", "PFKFB1") in glycolysis_graph.edges
    assert not ("PPP2CB", "sm_Pi") in glycolysis_graph.edges
    assert ("PPP2R5D", "PFKFB1") in glycolysis_graph.edges
    assert not ("PPP2R5D", "sm_Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-70262
    assert not ("PFKFB1", "sm_Fru(6)P") in glycolysis_graph.edges
    assert not ("PFKFB1", "sm_Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-71802
    assert not ("PFKFB3", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("PFKFB3", "sm_ADP") in glycolysis_graph.edges
    assert not ("PFKFB4", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("PFKFB4", "sm_ADP") in glycolysis_graph.edges
    assert not ("PFKFB1", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("PFKFB1", "sm_ADP") in glycolysis_graph.edges
    assert not ("PFKFB2", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("PFKFB2", "sm_ADP") in glycolysis_graph.edges

def test_connects_regulators_with_outputs(traffic_graph):

    # Regulator to output interactions for reaction R-HSA-5627071:
    assert ("RHOQ", "CFTR") in traffic_graph.edges
    assert ("GOPC", "CFTR") in traffic_graph.edges
    assert not ("CFTR", "CFTR") in traffic_graph.edges
    assert not ("sm_GTP", "CFTR") in traffic_graph.edges

    assert ("RHOQ", "GOPC") in traffic_graph.edges
    assert not ("GOPC", "GOPC") in traffic_graph.edges
    assert ("CFTR", "GOPC") in traffic_graph.edges
    assert not ("sm_GTP", "GOPC") in traffic_graph.edges

    assert not ("RHOQ", "RHOQ") in traffic_graph.edges
    assert ("GOPC", "RHOQ") in traffic_graph.edges
    assert ("CFTR", "RHOQ") in traffic_graph.edges
    assert not ("GTP", "RHOQ") in traffic_graph.edges

    assert not ("RHOQ", "GTP") in traffic_graph.edges
    assert not ("GOPC", "GTP") in traffic_graph.edges
    assert not ("CFTR", "GTP") in traffic_graph.edges
    assert not ("GTP", "GTP") in traffic_graph.edges

def test_connects_regulators_with_outputs_2(signal_graph):

    # In reaction R-HSA-5654408
    assert ("FGF20", "FGF16") or ("FGF16", "FGF20") in signal_graph.edges
    assert ("FGF17", "FGF20") in signal_graph.edges
    assert ("FGFR3", "FGF20") in signal_graph.edges
    assert not ("FGF23", "sm_HS") in signal_graph.edges
    assert not ("FGF23", "FGF23") in signal_graph.edges
    assert ("FGF23", "FRS2") in signal_graph.edges
    assert not ("sm_HS", "FGF17") in signal_graph.edges

def test_connects_components_of_same_complex(glycolysis_graph):

    # As part of PP2A-ABdeltaC complex R-HSA-165961
    assert not ("PPP2R1A", "PPP2R1A") in glycolysis_graph.edges
    assert ("PPP2R1A", "PPP2R1B") in glycolysis_graph.edges
    assert ("PPP2R1B", "PPP2R1A") in glycolysis_graph.edges
    assert ("PPP2R1A", "PPP2CA") in glycolysis_graph.edges
    assert ("PPP2R1A", "PPP2CB") in glycolysis_graph.edges
    assert ("PPP2R1A", "PPP2R5D") in glycolysis_graph.edges
    assert ("PPP2R5D", "PPP2R1B") in glycolysis_graph.edges
    assert ("PPP2R5D", "PPP2R1A") in glycolysis_graph.edges
    assert ("PPP2R5D", "PPP2CB") in glycolysis_graph.edges
