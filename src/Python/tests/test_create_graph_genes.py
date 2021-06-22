import pytest

from interactomes import create_graph

@pytest.fixture(scope="session")
def glycolysis_graph(tmpdir_factory):
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_graph("R-HSA-9634600", "genes", True, graphs_path)


@pytest.fixture(scope="session")
def traffic_graph(tmpdir_factory):
    # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_graph("R-HSA-5627083", "genes", True, graphs_path)


def test_create_graph_wrong_level_raises_exception():
    with pytest.raises(Exception):
        get_reaction_participants_by_pathway("R-HSA-9634600", "banana", True)


def test_create_graph_num_edges(glycolysis_graph):
    print(glycolysis_graph.edges)
    assert len(glycolysis_graph.edges) == 46


def test_connects_input_genes_with_small_outputs_not_when_is_same_molecule(glycolysis_graph):
    assert not ("PFKFB1", "PFKFB1") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-163773:
    assert ("PRKACA", "PFKFB1") in glycolysis_graph.edges
    assert ("PRKACA", "ADP") in glycolysis_graph.edges
    assert ("PRKACB", "PFKFB1") in glycolysis_graph.edges
    assert ("PRKACB", "ADP") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-163750
    assert ("PPP2R1A", "PFKFB1") in glycolysis_graph.edges
    assert ("PPP2R1A", "Pi") in glycolysis_graph.edges
    assert ("PPP2R1B", "PFKFB1") in glycolysis_graph.edges
    assert ("PPP2R1B", "Pi") in glycolysis_graph.edges
    assert ("PPP2CA", "PFKFB1") in glycolysis_graph.edges
    assert ("PPP2CA", "Pi") in glycolysis_graph.edges
    assert ("PPP2CB", "PFKFB1") in glycolysis_graph.edges
    assert ("PPP2CB", "Pi") in glycolysis_graph.edges
    assert ("PPP2R5D", "PFKFB1") in glycolysis_graph.edges
    assert ("PPP2R5D", "Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-70262
    assert ("PFKFB1", "Fru(6)P") in glycolysis_graph.edges
    assert ("PFKFB1", "Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-71802
    assert ("PFKFB3", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("PFKFB3", "ADP") in glycolysis_graph.edges
    assert ("PFKFB4", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("PFKFB4", "ADP") in glycolysis_graph.edges
    assert ("PFKFB1", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("PFKFB1", "ADP") in glycolysis_graph.edges
    assert ("PFKFB2", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("PFKFB2", "ADP") in glycolysis_graph.edges


def test_connects_regulators_with_outputs(traffic_graph):

    # Regulator to output interactions for reaction R-HSA-5627071:
    assert ("RHOQ", "CFTR") in traffic_graph.edges
    assert ("GOPC", "CFTR") in traffic_graph.edges
    assert not ("CFTR", "CFTR") in traffic_graph.edges
    assert ("GTP", "CFTR") in traffic_graph.edges

    assert ("RHOQ", "GOPC") in traffic_graph.edges
    assert not ("GOPC", "GOPC") in traffic_graph.edges
    assert ("CFTR", "GOPC") in traffic_graph.edges
    assert ("GTP", "GOPC") in traffic_graph.edges

    assert not ("RHOQ", "RHOQ") in traffic_graph.edges
    assert ("GOPC", "RHOQ") in traffic_graph.edges
    assert ("CFTR", "RHOQ") in traffic_graph.edges
    assert ("GTP", "RHOQ") in traffic_graph.edges

    assert ("RHOQ", "GTP") in traffic_graph.edges
    assert ("GOPC", "GTP") in traffic_graph.edges
    assert ("CFTR", "GTP") in traffic_graph.edges
    assert not ("GTP", "GTP") in traffic_graph.edges


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

#TODO: Test it connectects correctly the reaction: R-HSA-6805205