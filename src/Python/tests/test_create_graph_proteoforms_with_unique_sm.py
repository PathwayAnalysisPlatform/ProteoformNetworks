import pytest

from config import proteoforms, with_sm, with_unique_sm
from lib.networks import create_pathway_interaction_network


@pytest.fixture(scope="session")
def glycolysis_graph_with_unique_sm(tmpdir_factory):
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-9634600", proteoforms, with_unique_sm, graphs_path)


@pytest.fixture(scope="session")
def traffic_graph_with_unique_sm(tmpdir_factory):
    # Pathway "RHO sm_GTPases regulate P13569 trafficking" R-HSA-5627083
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-5627083", proteoforms, with_unique_sm, graphs_path)


@pytest.fixture(scope="session")
def signal_graph_with_unique_sm(tmpdir_factory):
    # Pathway "FRS-mediated FGFR3 signaling" R-HSA-5654706
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_pathway_interaction_network("R-HSA-5654706", proteoforms, with_unique_sm, graphs_path)

@pytest.fixture(scope="session")
def PI_3K_cascade_FGFR3_graph_with_unique_sm(tmpdir_factory):
    # Pathway "PI-3K cascade:FGFR3" R-HSA-5654710
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    g = create_pathway_interaction_network("R-HSA-5654710", proteoforms, with_unique_sm, graphs_path)
    return g


# Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
def test_connects_inputs_with_outputs(glycolysis_graph_with_unique_sm):
    # Input to output interactions for reaction R-HSA-163773:
    assert ("P16118;", "sm_R-HSA-163773_ADP") in glycolysis_graph_with_unique_sm.edges
    assert ("P16118;", "P16118;00046:33") in glycolysis_graph_with_unique_sm.edges
    assert ("sm_R-HSA-163773_ATP", "sm_R-HSA-163773_ADP") in glycolysis_graph_with_unique_sm.edges
    assert not ("sm_R-HSA-163773_ATP", "P16118;") in glycolysis_graph_with_unique_sm.edges
    assert ("sm_R-HSA-163773_ATP", "P16118;00046:33") in glycolysis_graph_with_unique_sm.edges

    # Input to output interactions for reaction R-HSA-163750
    assert not ("P16118;", "sm_R-HSA-163750_Pi") in glycolysis_graph_with_unique_sm.edges
    assert ("P16118;00046:33", "sm_R-HSA-163750_Pi") in glycolysis_graph_with_unique_sm.edges
    assert ("P16118;00046:33", "P16118;") in glycolysis_graph_with_unique_sm.edges
    assert ("sm_R-HSA-163750_H2O", "sm_R-HSA-163750_Pi") in glycolysis_graph_with_unique_sm.edges
    assert ("sm_R-HSA-163750_H2O", "P16118;") in glycolysis_graph_with_unique_sm.edges

    # Input to output interaction for reaction R-HSA-71802
    assert ("sm_R-HSA-71802_Fru(6)P", "sm_R-HSA-71802_D-Fructose_2,6-bisphosphate") in glycolysis_graph_with_unique_sm.edges
    assert ("sm_R-HSA-71802_Fru(6)P", "sm_R-HSA-71802_ADP") in glycolysis_graph_with_unique_sm.edges
    assert ("sm_R-HSA-71802_ATP", "sm_R-HSA-71802_D-Fructose_2,6-bisphosphate") in glycolysis_graph_with_unique_sm.edges
    assert ("sm_R-HSA-71802_ATP", "sm_R-HSA-71802_ADP") in glycolysis_graph_with_unique_sm.edges


def test_does_not_connect_input_to_output_when_same_molecule(traffic_graph_with_unique_sm):

    # In reaction R-HSA-5627071 there are same set of proteoforms as input and output, since they are in
    # different sub cellular locations. Nevertheless, there are no self interactions in our network
    assert not ("P13569;", "P13569;") in traffic_graph_with_unique_sm.edges
    assert not ("P17081;", "P17081;") in traffic_graph_with_unique_sm.edges
    assert not ("Q9HD26;", "Q9HD26;") in traffic_graph_with_unique_sm.edges
    assert not ("sm_R-HSA-5627071_GTP", "sm_R-HSA-5627071_GTP") in traffic_graph_with_unique_sm.edges


def test_connects_catalysts_with_outputs(glycolysis_graph_with_unique_sm):

    # Catalysts to output interactions for reaction R-HSA-163773:
    assert not ("P17612;", "P16118;") in glycolysis_graph_with_unique_sm.edges
    assert ("P17612;", "P16118;00046:33") in glycolysis_graph_with_unique_sm.edges
    assert ("P17612;", "sm_R-HSA-163773_ADP") in glycolysis_graph_with_unique_sm.edges
    assert not ("P22694;", "P16118;") in glycolysis_graph_with_unique_sm.edges
    assert ("P22694;", "P16118;00046:33") in glycolysis_graph_with_unique_sm.edges
    assert ("P22694;", "sm_R-HSA-163773_ADP") in glycolysis_graph_with_unique_sm.edges

    # Catalysts to output interactions for reaction R-HSA-163750
    assert ("P30153;", "P16118;") in glycolysis_graph_with_unique_sm.edges
    assert ("P30153;", "sm_R-HSA-163750_Pi") in glycolysis_graph_with_unique_sm.edges
    assert ("P30154;", "P16118;") in glycolysis_graph_with_unique_sm.edges
    assert ("P30154;", "sm_R-HSA-163750_Pi") in glycolysis_graph_with_unique_sm.edges
    assert ("P67775;", "P16118;") in glycolysis_graph_with_unique_sm.edges
    assert ("P67775;", "sm_R-HSA-163750_Pi") in glycolysis_graph_with_unique_sm.edges
    assert ("P62714;", "P16118;") in glycolysis_graph_with_unique_sm.edges
    assert ("P62714;", "sm_R-HSA-163750_Pi") in glycolysis_graph_with_unique_sm.edges
    assert ("Q14738;", "P16118;") in glycolysis_graph_with_unique_sm.edges
    assert ("Q14738;", "sm_R-HSA-163750_Pi") in glycolysis_graph_with_unique_sm.edges

    # Catalysts to output interactions for reaction R-HSA-70262
    assert not ("P16118;", "sm_R-HSA-70262_Fru(6)P") in glycolysis_graph_with_unique_sm.edges
    assert ("P16118;00046:33", "sm_R-HSA-70262_Fru(6)P") in glycolysis_graph_with_unique_sm.edges
    assert not ("P16118;", "sm_R-HSA-70262_Pi") in glycolysis_graph_with_unique_sm.edges
    assert ("P16118;00046:33", "sm_R-HSA-70262_Pi") in glycolysis_graph_with_unique_sm.edges

    # Catalysts to output interactions for reaction R-HSA-71802
    assert ("Q16875;", "sm_R-HSA-71802_D-Fructose_2,6-bisphosphate") in glycolysis_graph_with_unique_sm.edges
    assert ("Q16875;", "sm_R-HSA-71802_ADP") in glycolysis_graph_with_unique_sm.edges
    assert ("Q16877;", "sm_R-HSA-71802_D-Fructose_2,6-bisphosphate") in glycolysis_graph_with_unique_sm.edges
    assert ("Q16877;", "sm_R-HSA-71802_ADP") in glycolysis_graph_with_unique_sm.edges
    assert ("P16118;", "sm_R-HSA-71802_D-Fructose_2,6-bisphosphate") in glycolysis_graph_with_unique_sm.edges
    assert ("P16118;", "sm_R-HSA-71802_ADP") in glycolysis_graph_with_unique_sm.edges
    assert ("O60825;", "sm_R-HSA-71802_D-Fructose_2,6-bisphosphate") in glycolysis_graph_with_unique_sm.edges
    assert ("O60825;", "sm_R-HSA-71802_ADP") in glycolysis_graph_with_unique_sm.edges


def test_connects_regulators_with_outputs(traffic_graph_with_unique_sm):

    # Regulator to output interactions for reaction R-HSA-5627071:
    assert ("P17081;", "P13569;") in traffic_graph_with_unique_sm.edges
    assert ("Q9HD26;", "P13569;") in traffic_graph_with_unique_sm.edges
    assert not ("P13569;", "P13569;") in traffic_graph_with_unique_sm.edges
    assert ("sm_R-HSA-5627071_GTP", "P13569;") in traffic_graph_with_unique_sm.edges

    assert ("P17081;", "Q9HD26;") in traffic_graph_with_unique_sm.edges
    assert not ("Q9HD26;", "Q9HD26;") in traffic_graph_with_unique_sm.edges
    assert ("P13569;", "Q9HD26;") in traffic_graph_with_unique_sm.edges
    assert ("sm_R-HSA-5627071_GTP", "Q9HD26;") in traffic_graph_with_unique_sm.edges

    assert not ("P17081;", "P17081;") in traffic_graph_with_unique_sm.edges
    assert ("Q9HD26;", "P17081;") in traffic_graph_with_unique_sm.edges
    assert ("P13569;", "P17081;") in traffic_graph_with_unique_sm.edges
    assert ("sm_R-HSA-5627071_GTP", "P17081;") in traffic_graph_with_unique_sm.edges

    assert ("P17081;", "sm_R-HSA-5627071_GTP") in traffic_graph_with_unique_sm.edges
    assert ("Q9HD26;", "sm_R-HSA-5627071_GTP") in traffic_graph_with_unique_sm.edges
    assert ("P13569;", "sm_R-HSA-5627071_GTP") in traffic_graph_with_unique_sm.edges
    assert not ("sm_R-HSA-5627071_GTP", "sm_R-HSA-5627071_GTP") in traffic_graph_with_unique_sm.edges


def test_connects_regulators_with_outputs_2(signal_graph_with_unique_sm):

    # In reaction R-HSA-5654408
    assert ("Q9NP95;", "O43320;") or ("O43320;", "Q9NP95;") in signal_graph_with_unique_sm.edges
    assert ("O60258-1;", "Q9NP95;") in signal_graph_with_unique_sm.edges
    assert ("P22607-1;00048:577,00048:647,00048:648,00048:724,00048:760,00048:770", "Q9NP95;") in signal_graph_with_unique_sm.edges
    assert ("Q9GZV9;00164:178", "sm_R-HSA-5654408_HS") in signal_graph_with_unique_sm.edges
    assert not ("Q9GZV9;00164:178", "Q9GZV9;00164:178") in signal_graph_with_unique_sm.edges
    assert ("Q9GZV9;00164:178", "Q8WU20;00048:196,00048:306,00048:349,00048:392,00048:436,00048:471") in signal_graph_with_unique_sm.edges
    assert ("sm_R-HSA-5654408_HS", "O60258-1;") in signal_graph_with_unique_sm.edges


def test_connects_components_of_same_complex(glycolysis_graph_with_unique_sm):

    # As part of PP2A-ABdeltaC complex R-HSA-165961
    assert not ("P30153;", "P30153;") in glycolysis_graph_with_unique_sm.edges
    assert ("P30153;", "P30154;") in glycolysis_graph_with_unique_sm.edges
    assert ("P30154;", "P30153;") in glycolysis_graph_with_unique_sm.edges
    assert ("P30153;", "P67775;") in glycolysis_graph_with_unique_sm.edges
    assert ("P30153;", "P62714;") in glycolysis_graph_with_unique_sm.edges
    assert ("P30153;", "Q14738;") in glycolysis_graph_with_unique_sm.edges
    assert ("Q14738;", "P30154;") in glycolysis_graph_with_unique_sm.edges
    assert ("Q14738;", "P30153;") in glycolysis_graph_with_unique_sm.edges
    assert ("Q14738;", "P62714;") in glycolysis_graph_with_unique_sm.edges


def test_connect_components_and_participants2(PI_3K_cascade_FGFR3_graph_with_unique_sm):

    # For reaction R-HSA-177931
    assert ("P62993-1;", "P27986;") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("P62993-1;", "Q13480;") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert not ("P62993-1;", "P62993-1;") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("Q13480;", "P27986;") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert not ("Q13480;", "Q13480;") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("Q13480;", "P62993-1;") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges

    # Reaction R-HSA-5654705 "FGFR3-associated PI3K phosphorylates PIP2 to PIP3"
    # Catalysts to output
    assert ("Q8WU20;00048:196,00048:306,00048:349,00048:392,00048:436,00048:471", "sm_R-HSA-5654705_ADP") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("Q9GZV9;00164:178", "sm_R-HSA-5654705_ADP") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("P55075-1;", "sm_R-HSA-5654705_ADP") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges

    # Catalysts to output
    assert ("Q8WU20;00048:196,00048:306,00048:349,00048:392,00048:436,00048:471", "sm_R-HSA-5654705_PI(3,4,5)P3") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("Q9GZV9;00164:178", "sm_R-HSA-5654705_PI(3,4,5)P3") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("P55075-1;", "sm_R-HSA-5654705_PI(3,4,5)P3") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges

    # Input to Output
    assert ("sm_R-HSA-5654705_ATP", "sm_R-HSA-5654705_ADP") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("sm_R-HSA-5654705_PI(4,5)P2", "sm_R-HSA-5654705_ADP") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges

    # Components of the same complex
    assert ("P05230;", "O60258-1;") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("P12034-1;", "P55075-1;") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges

    # For Reaction R-HSA-5654709 "FGFR3- and PTPN11-associated PI3K phosphorylates PIP2 to PIP3"
    assert ("sm_R-HSA-5654709_ATP", "sm_R-HSA-5654709_ADP") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("sm_R-HSA-5654709_ATP", "sm_R-HSA-5654709_PI(3,4,5)P3") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("sm_R-HSA-5654709_PI(4,5)P2", "sm_R-HSA-5654709_ADP") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges
    assert ("sm_R-HSA-5654709_PI(4,5)P2", "sm_R-HSA-5654709_PI(3,4,5)P3") in PI_3K_cascade_FGFR3_graph_with_unique_sm.edges


def test_nodes_have_protein_attribute(glycolysis_graph_with_unique_sm):
    for node in glycolysis_graph_with_unique_sm.nodes:
        assert 'prevId' in glycolysis_graph_with_unique_sm.nodes[node]

def test_nodes_have_correct_protein_attribute(glycolysis_graph_with_unique_sm):
    assert glycolysis_graph_with_unique_sm.nodes['P16118;00046:33']['prevId'] == 'P16118'
    assert glycolysis_graph_with_unique_sm.nodes['sm_R-HSA-71802_ADP']['prevId'] == 'sm_R-HSA-71802_ADP'
    assert glycolysis_graph_with_unique_sm.nodes['P62714;']['prevId'] == 'P62714'