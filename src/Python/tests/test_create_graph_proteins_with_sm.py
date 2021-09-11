import pytest

from config import proteins, with_sm
from lib.networks import create_pathway_interaction_network


@pytest.fixture(scope="session")
def glycolysis_graph(tmpdir_factory):
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-9634600", proteins, with_sm, graphs_path)


@pytest.fixture(scope="session")
def traffic_graph(tmpdir_factory):
    # Pathway "RHO sm_GTPases regulate P13569 trafficking" R-HSA-5627083
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-5627083", proteins, with_sm, graphs_path)


def test_create_graph_num_edges(glycolysis_graph):
    print(glycolysis_graph.edges)
    assert len(glycolysis_graph.edges) == 46


def test_create_graph_num_vertices(glycolysis_graph):
    print(glycolysis_graph.nodes)
    assert len(glycolysis_graph.nodes) == 18


# Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
def test_connects_inputs_with_outputs(glycolysis_graph):
    # Input to output interactions for reaction R-HSA-163773:
    assert ("P16118", "sm_ADP") in glycolysis_graph.edges
    assert ("sm_ATP", "sm_ADP") in glycolysis_graph.edges
    assert ("sm_ATP", "P16118") in glycolysis_graph.edges

    # Input to output interactions for reaction R-HSA-163750
    assert ("P16118", "sm_Pi") in glycolysis_graph.edges
    assert ("sm_H2O", "sm_Pi") in glycolysis_graph.edges
    assert ("sm_H2O", "P16118") in glycolysis_graph.edges

    # Input to output interaction for reaction R-HSA-71802
    assert ("sm_Fru(6)P", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("sm_Fru(6)P", "sm_ADP") in glycolysis_graph.edges
    assert ("sm_ATP", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("sm_ATP", "sm_ADP") in glycolysis_graph.edges


def test_connects_input_genes_with_small_outputs_not_when_is_same_molecule(glycolysis_graph):
    assert not ("P16118", "P16118") in glycolysis_graph.edges


def test_connects_catalysts_with_outputs(glycolysis_graph):
    # Catalysts to output interactions for reaction R-HSA-163773:
    assert ("P17612", "P16118") in glycolysis_graph.edges
    assert ("P17612", "sm_ADP") in glycolysis_graph.edges
    assert ("P22694", "P16118") in glycolysis_graph.edges
    assert ("P22694", "sm_ADP") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-163750
    assert ("P30153", "P16118") in glycolysis_graph.edges
    assert ("P30153", "sm_Pi") in glycolysis_graph.edges
    assert ("P30154", "P16118") in glycolysis_graph.edges
    assert ("P30154", "sm_Pi") in glycolysis_graph.edges
    assert ("P67775", "P16118") in glycolysis_graph.edges
    assert ("P67775", "sm_Pi") in glycolysis_graph.edges
    assert ("P62714", "P16118") in glycolysis_graph.edges
    assert ("P62714", "sm_Pi") in glycolysis_graph.edges
    assert ("Q14738", "P16118") in glycolysis_graph.edges
    assert ("Q14738", "sm_Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-70262
    assert ("P16118", "sm_Fru(6)P") in glycolysis_graph.edges
    assert ("P16118", "sm_Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-71802
    assert ("Q16875", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("Q16875", "sm_ADP") in glycolysis_graph.edges
    assert ("Q16877", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("Q16877", "sm_ADP") in glycolysis_graph.edges
    assert ("P16118", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("P16118", "sm_ADP") in glycolysis_graph.edges
    assert ("O60825", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("O60825", "sm_ADP") in glycolysis_graph.edges


def test_connects_regulators_with_outputs(traffic_graph):
    # Regulator to output interactions for reaction R-HSA-5627071:
    assert ("P17081", "P13569") in traffic_graph.edges
    assert ("Q9HD26", "P13569") in traffic_graph.edges
    assert not ("P13569", "P13569") in traffic_graph.edges
    assert ("sm_GTP", "P13569") in traffic_graph.edges

    assert ("P17081", "Q9HD26") in traffic_graph.edges
    assert not ("Q9HD26", "Q9HD26") in traffic_graph.edges
    assert ("P13569", "Q9HD26") in traffic_graph.edges
    assert ("sm_GTP", "Q9HD26") in traffic_graph.edges

    assert not ("P17081", "P17081") in traffic_graph.edges
    assert ("Q9HD26", "P17081") in traffic_graph.edges
    assert ("P13569", "P17081") in traffic_graph.edges
    assert ("sm_GTP", "P17081") in traffic_graph.edges

    assert ("P17081", "sm_GTP") in traffic_graph.edges
    assert ("Q9HD26", "sm_GTP") in traffic_graph.edges
    assert ("P13569", "sm_GTP") in traffic_graph.edges
    assert not ("sm_GTP", "sm_GTP") in traffic_graph.edges


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


def test_nodes_have_gene_attribute(glycolysis_graph):
    for node in glycolysis_graph.nodes:
        assert 'prevId' in glycolysis_graph.nodes[node]

def test_nodes_have_correct_gene_attribute(glycolysis_graph):
    assert glycolysis_graph.nodes['P30153']['prevId'] == 'PPP2R1A'
    assert glycolysis_graph.nodes['sm_ADP']['prevId'] == 'sm_ADP'
    assert glycolysis_graph.nodes['P16118']['prevId'] == 'PFKFB1'
