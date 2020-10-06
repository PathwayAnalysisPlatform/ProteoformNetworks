from networks import create_graph
import pytest


@pytest.fixture(scope="session")
def glycolysis_graph(tmpdir_factory):
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_graph("R-HSA-9634600", "genes", False, graphs_path)


@pytest.fixture(scope="session")
def traffic_graph(tmpdir_factory):
    # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_graph("R-HSA-5627083", "genes", False, graphs_path)

@pytest.fixture(scope="session")
def signal_graph(tmpdir_factory):
    # Pathway "FRS-mediated FGFR3 signaling"
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_graph("R-HSA-5654706", "genes", False, graphs_path)

def test_create_graph_num_edges(glycolysis_graph):
    print(glycolysis_graph.edges)
    assert len(glycolysis_graph.edges) == 18

def test_create_graph_num_vertices(glycolysis_graph):
    print(glycolysis_graph.nodes)
    assert len(glycolysis_graph.nodes) == 12

# Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
def test_connects_inputs_with_outputs(glycolysis_graph):

    # Input to output interactions for reaction R-HSA-163750
    assert not ("PFKFB1", "Pi") in glycolysis_graph.edges
    assert not ("PFKFB1", "PFKFB1") in glycolysis_graph.edges
    assert not ("H2O", "Pi") in glycolysis_graph.edges
    assert not ("H2O", "PFKFB1") in glycolysis_graph.edges

    # Input to output interactions for reaction R-HSA-163773:
    assert not ("PFKFB1", "ADP") in glycolysis_graph.edges
    assert not ("ATP", "ADP") in glycolysis_graph.edges
    assert not ("ATP", "PFKFB1") in glycolysis_graph.edges

    # Input to output interaction for reaction R-HSA-71802
    assert not ("Fru(6)P", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("Fru(6)P", "ADP") in glycolysis_graph.edges
    assert not ("ATP", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("ATP", "ADP") in glycolysis_graph.edges

def test_connects_catalysts_with_outputs(glycolysis_graph):

    # Catalysts to output interactions for reaction R-HSA-163773:
    assert ("PRKACA", "PFKFB1") in glycolysis_graph.edges
    assert ("PRKACG", "PFKFB1") in glycolysis_graph.edges
    assert ("PRKACB", "PFKFB1") in glycolysis_graph.edges
    assert not ("PRKACA", "ADP") in glycolysis_graph.edges
    assert not ("PRKACB", "ADP") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-163750
    assert ("PPP2R1A", "PFKFB1") in glycolysis_graph.edges
    assert not ("PPP2R1A", "Pi") in glycolysis_graph.edges
    assert ("PPP2R1B", "PFKFB1") in glycolysis_graph.edges
    assert not ("PPP2R1B", "Pi") in glycolysis_graph.edges
    assert ("PPP2CA", "PFKFB1") in glycolysis_graph.edges
    assert not ("PPP2CA", "Pi") in glycolysis_graph.edges
    assert ("PPP2CB", "PFKFB1") in glycolysis_graph.edges
    assert not ("PPP2CB", "Pi") in glycolysis_graph.edges
    assert ("PPP2R5D", "PFKFB1") in glycolysis_graph.edges
    assert not ("PPP2R5D", "Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-70262
    assert not ("PFKFB1", "Fru(6)P") in glycolysis_graph.edges
    assert not ("PFKFB1", "Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-71802
    assert not ("PFKFB3", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("PFKFB3", "ADP") in glycolysis_graph.edges
    assert not ("PFKFB4", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("PFKFB4", "ADP") in glycolysis_graph.edges
    assert not ("PFKFB1", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("PFKFB1", "ADP") in glycolysis_graph.edges
    assert not ("PFKFB2", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("PFKFB2", "ADP") in glycolysis_graph.edges

def test_connects_regulators_with_outputs(traffic_graph):

    # Regulator to output interactions for reaction R-HSA-5627071:
    assert ("RHOQ", "CFTR") in traffic_graph.edges
    assert ("GOPC", "CFTR") in traffic_graph.edges
    assert not ("CFTR", "CFTR") in traffic_graph.edges
    assert not ("GTP", "CFTR") in traffic_graph.edges

    assert ("RHOQ", "GOPC") in traffic_graph.edges
    assert not ("GOPC", "GOPC") in traffic_graph.edges
    assert ("CFTR", "GOPC") in traffic_graph.edges
    assert not ("GTP", "GOPC") in traffic_graph.edges

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
    assert not ("FGF23", "HS") in signal_graph.edges
    assert not ("FGF23", "FGF23") in signal_graph.edges
    assert ("FGF23", "FRS2") in signal_graph.edges
    assert not ("HS", "FGF17") in signal_graph.edges

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
