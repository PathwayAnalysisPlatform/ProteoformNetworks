import pytest

from interaction_network import create_graph


@pytest.fixture(scope="session")
def glycolysis_graph(tmpdir_factory):
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_graph("R-HSA-9634600", "proteins", True, graphs_path)


@pytest.fixture(scope="session")
def traffic_graph(tmpdir_factory):
    # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_graph("R-HSA-5627083", "proteins", True, graphs_path)


def test_create_graph_num_edges(glycolysis_graph):
    print(glycolysis_graph.edges)
    assert len(glycolysis_graph.edges) == 46


def test_create_graph_num_vertices(glycolysis_graph):
    print(glycolysis_graph.nodes)
    assert len(glycolysis_graph.nodes) == 18


# Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
def test_connects_inputs_with_outputs(glycolysis_graph):
    # Input to output interactions for reaction R-HSA-163773:
    assert ("P16118", "ADP") in glycolysis_graph.edges
    assert ("ATP", "ADP") in glycolysis_graph.edges
    assert ("ATP", "P16118") in glycolysis_graph.edges

    # Input to output interactions for reaction R-HSA-163750
    assert ("P16118", "Pi") in glycolysis_graph.edges
    assert ("H2O", "Pi") in glycolysis_graph.edges
    assert ("H2O", "P16118") in glycolysis_graph.edges

    # Input to output interaction for reaction R-HSA-71802
    assert ("Fru(6)P", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("Fru(6)P", "ADP") in glycolysis_graph.edges
    assert ("ATP", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("ATP", "ADP") in glycolysis_graph.edges


def test_connects_input_genes_with_small_outputs_not_when_is_same_molecule(glycolysis_graph):
    assert not ("P16118", "P16118") in glycolysis_graph.edges


def test_connects_catalysts_with_outputs(glycolysis_graph):

    # Catalysts to output interactions for reaction R-HSA-163773:
    assert ("P17612", "P16118") in glycolysis_graph.edges
    assert ("P17612", "ADP") in glycolysis_graph.edges
    assert ("P22694", "P16118") in glycolysis_graph.edges
    assert ("P22694", "ADP") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-163750
    assert ("P30153", "P16118") in glycolysis_graph.edges
    assert ("P30153", "Pi") in glycolysis_graph.edges
    assert ("P30154", "P16118") in glycolysis_graph.edges
    assert ("P30154", "Pi") in glycolysis_graph.edges
    assert ("P67775", "P16118") in glycolysis_graph.edges
    assert ("P67775", "Pi") in glycolysis_graph.edges
    assert ("P62714", "P16118") in glycolysis_graph.edges
    assert ("P62714", "Pi") in glycolysis_graph.edges
    assert ("Q14738", "P16118") in glycolysis_graph.edges
    assert ("Q14738", "Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-70262
    assert ("P16118", "Fru(6)P") in glycolysis_graph.edges
    assert ("P16118", "Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-71802
    assert ("Q16875", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("Q16875", "ADP") in glycolysis_graph.edges
    assert ("Q16877", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("Q16877", "ADP") in glycolysis_graph.edges
    assert ("P16118", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("P16118", "ADP") in glycolysis_graph.edges
    assert ("O60825", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("O60825", "ADP") in glycolysis_graph.edges


def test_connects_regulators_with_outputs(traffic_graph):
    # Regulator to output interactions for reaction R-HSA-5627071:
    assert ("P17081", "P13569") in traffic_graph.edges
    assert ("Q9HD26", "P13569") in traffic_graph.edges
    assert not ("P13569", "P13569") in traffic_graph.edges
    assert ("GTP", "P13569") in traffic_graph.edges

    assert ("P17081", "Q9HD26") in traffic_graph.edges
    assert not ("Q9HD26", "Q9HD26") in traffic_graph.edges
    assert ("P13569", "Q9HD26") in traffic_graph.edges
    assert ("GTP", "Q9HD26") in traffic_graph.edges

    assert not ("P17081", "P17081") in traffic_graph.edges
    assert ("Q9HD26", "P17081") in traffic_graph.edges
    assert ("P13569", "P17081") in traffic_graph.edges
    assert ("GTP", "P17081") in traffic_graph.edges

    assert ("P17081", "GTP") in traffic_graph.edges
    assert ("Q9HD26", "GTP") in traffic_graph.edges
    assert ("P13569", "GTP") in traffic_graph.edges
    assert not ("GTP", "GTP") in traffic_graph.edges


def test_connects_components_of_same_complex(glycolysis_graph):
    # As part of PP2A-ABdeltaC complex R-HSA-165961
    assert not ("P30153", "P30153") in glycolysis_graph.edges
    assert ("P30153", "P30154") in glycolysis_graph.edges
    assert ("P30154", "P30153") in glycolysis_graph.edges
    assert ("P30153", "P67775") in glycolysis_graph.edges
    assert ("P30153", "P62714") in glycolysis_graph.edges
    assert ("P30153", "Q14738") in glycolysis_graph.edges
    assert ("Q14738", "P30154") in glycolysis_graph.edges
    assert ("Q14738", "P30153") in glycolysis_graph.edges
    assert ("Q14738", "P62714") in glycolysis_graph.edges
