import pytest

from config import proteoforms, with_sm
from lib.networks import create_pathway_interaction_network


@pytest.fixture(scope="session")
def glycolysis_graph(tmpdir_factory):
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-9634600", proteoforms, with_sm, graphs_path)


@pytest.fixture(scope="session")
def traffic_graph(tmpdir_factory):
    # Pathway "RHO sm_GTPases regulate P13569 trafficking" R-HSA-5627083
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-5627083", proteoforms, with_sm, graphs_path)


@pytest.fixture(scope="session")
def signal_graph(tmpdir_factory):
    # Pathway "FRS-mediated FGFR3 signaling" R-HSA-5654706
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-5654706", proteoforms, with_sm, graphs_path)


def test_create_graph_num_edges(glycolysis_graph):
    print(glycolysis_graph.edges)
    assert len(glycolysis_graph.edges) == 47


def test_create_graph_num_vertices(glycolysis_graph):
    print(glycolysis_graph.nodes)
    assert len(glycolysis_graph.nodes) == 19


# Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
def test_connects_inputs_with_outputs(glycolysis_graph):
    # Input to output interactions for reaction R-HSA-163773:
    assert ("P16118;", "sm_ADP") in glycolysis_graph.edges
    assert ("P16118;", "P16118;00046:33") in glycolysis_graph.edges
    assert ("sm_ATP", "sm_ADP") in glycolysis_graph.edges
    assert not ("sm_ATP", "P16118;") in glycolysis_graph.edges
    assert ("sm_ATP", "P16118;00046:33") in glycolysis_graph.edges

    # Input to output interactions for reaction R-HSA-163750
    assert not ("P16118;", "sm_Pi") in glycolysis_graph.edges
    assert ("P16118;00046:33", "sm_Pi") in glycolysis_graph.edges
    assert ("P16118;00046:33", "P16118;") in glycolysis_graph.edges
    assert ("sm_H2O", "sm_Pi") in glycolysis_graph.edges
    assert ("sm_H2O", "P16118;") in glycolysis_graph.edges

    # Input to output interaction for reaction R-HSA-71802
    assert ("sm_Fru(6)P", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("sm_Fru(6)P", "sm_ADP") in glycolysis_graph.edges
    assert ("sm_ATP", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("sm_ATP", "sm_ADP") in glycolysis_graph.edges


def test_does_not_connect_input_to_output_when_same_molecule(traffic_graph):

    # In reaction R-HSA-5627071 there are same set of proteoforms as input and output, since they are in
    # different sub cellular locations. Nevertheless, there are no self interactions in our network
    assert not ("P13569;", "P13569;") in traffic_graph.edges
    assert not ("P17081;", "P17081;") in traffic_graph.edges
    assert not ("Q9HD26;", "Q9HD26;") in traffic_graph.edges
    assert not ("sm_GTP", "sm_GTP") in traffic_graph.edges


def test_connects_catalysts_with_outputs(glycolysis_graph):

    # Catalysts to output interactions for reaction R-HSA-163773:
    assert not ("P17612;", "P16118;") in glycolysis_graph.edges
    assert ("P17612;", "P16118;00046:33") in glycolysis_graph.edges
    assert ("P17612;", "sm_ADP") in glycolysis_graph.edges
    assert not ("P22694;", "P16118;") in glycolysis_graph.edges
    assert ("P22694;", "P16118;00046:33") in glycolysis_graph.edges
    assert ("P22694;", "sm_ADP") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-163750
    assert ("P30153;", "P16118;") in glycolysis_graph.edges
    assert ("P30153;", "sm_Pi") in glycolysis_graph.edges
    assert ("P30154;", "P16118;") in glycolysis_graph.edges
    assert ("P30154;", "sm_Pi") in glycolysis_graph.edges
    assert ("P67775;", "P16118;") in glycolysis_graph.edges
    assert ("P67775;", "sm_Pi") in glycolysis_graph.edges
    assert ("P62714;", "P16118;") in glycolysis_graph.edges
    assert ("P62714;", "sm_Pi") in glycolysis_graph.edges
    assert ("Q14738;", "P16118;") in glycolysis_graph.edges
    assert ("Q14738;", "sm_Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-70262
    assert not ("P16118;", "sm_Fru(6)P") in glycolysis_graph.edges
    assert ("P16118;00046:33", "sm_Fru(6)P") in glycolysis_graph.edges
    assert not ("P16118;", "sm_Pi") in glycolysis_graph.edges
    assert ("P16118;00046:33", "sm_Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-71802
    assert ("Q16875;", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("Q16875;", "sm_ADP") in glycolysis_graph.edges
    assert ("Q16877;", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("Q16877;", "sm_ADP") in glycolysis_graph.edges
    assert ("P16118;", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("P16118;", "sm_ADP") in glycolysis_graph.edges
    assert ("O60825;", "sm_D-Fructose_2,6-bisphosphate") in glycolysis_graph.edges
    assert ("O60825;", "sm_ADP") in glycolysis_graph.edges


def test_connects_regulators_with_outputs(traffic_graph):

    # Regulator to output interactions for reaction R-HSA-5627071:
    assert ("P17081;", "P13569;") in traffic_graph.edges
    assert ("Q9HD26;", "P13569;") in traffic_graph.edges
    assert not ("P13569;", "P13569;") in traffic_graph.edges
    assert ("sm_GTP", "P13569;") in traffic_graph.edges

    assert ("P17081;", "Q9HD26;") in traffic_graph.edges
    assert not ("Q9HD26;", "Q9HD26;") in traffic_graph.edges
    assert ("P13569;", "Q9HD26;") in traffic_graph.edges
    assert ("sm_GTP", "Q9HD26;") in traffic_graph.edges

    assert not ("P17081;", "P17081;") in traffic_graph.edges
    assert ("Q9HD26;", "P17081;") in traffic_graph.edges
    assert ("P13569;", "P17081;") in traffic_graph.edges
    assert ("sm_GTP", "P17081;") in traffic_graph.edges

    assert ("P17081;", "sm_GTP") in traffic_graph.edges
    assert ("Q9HD26;", "sm_GTP") in traffic_graph.edges
    assert ("P13569;", "sm_GTP") in traffic_graph.edges
    assert not ("sm_GTP", "sm_GTP") in traffic_graph.edges


def test_connects_regulators_with_outputs_2(signal_graph):

    # In reaction R-HSA-5654408
    assert ("Q9NP95;", "O43320;") or ("O43320;", "Q9NP95;") in signal_graph.edges
    assert ("O60258-1;", "Q9NP95;") in signal_graph.edges
    assert ("P22607-1;00048:577,00048:647,00048:648,00048:724,00048:760,00048:770", "Q9NP95;") in signal_graph.edges
    assert ("Q9GZV9;00164:178", "sm_HS") in signal_graph.edges
    assert not ("Q9GZV9;00164:178", "Q9GZV9;00164:178") in signal_graph.edges
    assert ("Q9GZV9;00164:178", "Q8WU20;00048:196,00048:306,00048:349,00048:392,00048:436,00048:471") in signal_graph.edges
    assert ("sm_HS", "O60258-1;") in signal_graph.edges


def test_connects_components_of_same_complex(glycolysis_graph):

    # As part of PP2A-ABdeltaC complex R-HSA-165961
    assert not ("P30153;", "P30153;") in glycolysis_graph.edges
    assert ("P30153;", "P30154;") in glycolysis_graph.edges
    assert ("P30154;", "P30153;") in glycolysis_graph.edges
    assert ("P30153;", "P67775;") in glycolysis_graph.edges
    assert ("P30153;", "P62714;") in glycolysis_graph.edges
    assert ("P30153;", "Q14738;") in glycolysis_graph.edges
    assert ("Q14738;", "P30154;") in glycolysis_graph.edges
    assert ("Q14738;", "P30153;") in glycolysis_graph.edges
    assert ("Q14738;", "P62714;") in glycolysis_graph.edges


def test_connect_components_of_same_complex_2(signal_graph):
    assert ("P01116-1;00115:179,01116:186", "sm_GTP") in signal_graph.edges
    assert ("P01111;00115:181,01116:186", "sm_GTP") in signal_graph.edges
    assert ("P01112;00115:181,00115:184,01116:186", "sm_GTP") in signal_graph.edges
    assert ("P01116-2;01116:186", "sm_GTP") in signal_graph.edges

    assert ("P01112;00115:181,00115:184,01116:186", "P01111;00115:181,01116:186") in signal_graph.edges
    assert ("P01112;00115:181,00115:184,01116:186", "P01116-2;01116:186") in signal_graph.edges

    assert ("P01116-2;01116:186", "P01111;00115:181,01116:186") in signal_graph.edges
    assert ("P01116-2;01116:186", "P01112;00115:181,00115:184,01116:186") in signal_graph.edges

    assert ("P01111;00115:181,01116:186", "P01116-2;01116:186") in signal_graph.edges
    assert ("P01111;00115:181,01116:186", "P01112;00115:181,00115:184,01116:186") in signal_graph.edges

    assert ("Q07889;", "P62993-1;") in signal_graph.edges

    assert ("sm_HS", "P22607-2;00048:579,00048:649,00048:650,00048:726,00048:762,00048:772") in signal_graph.edges
    assert ("P22607-2;00048:579,00048:649,00048:650,00048:726,00048:762,00048:772", "P05230;") in signal_graph.edges

    assert ("P55075-1;", "P09038;") in signal_graph.edges
    assert ("O60258-1;", "Q9GZV9;00164:178") in signal_graph.edges


def test_nodes_have_protein_attribute(glycolysis_graph):
    for node in glycolysis_graph.nodes:
        assert 'prevId' in glycolysis_graph.nodes[node]

def test_nodes_have_correct_protein_attribute(glycolysis_graph):
    assert glycolysis_graph.nodes['P16118;00046:33']['prevId'] == 'P16118'
    assert glycolysis_graph.nodes['sm_ADP']['prevId'] == 'sm_ADP'
    assert glycolysis_graph.nodes['P62714;']['prevId'] == 'P62714'