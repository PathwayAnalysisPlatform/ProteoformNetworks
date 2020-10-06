import pytest

from networks import create_graph


@pytest.fixture(scope="session")
def glycolysis_graph(tmpdir_factory):
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_graph("R-HSA-9634600", "proteoforms", False, graphs_path)


@pytest.fixture(scope="session")
def traffic_graph(tmpdir_factory):
    # Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_graph("R-HSA-5627083", "proteoforms", False, graphs_path)


@pytest.fixture(scope="session")
def signal_graph(tmpdir_factory):
    # Pathway "FRS-mediated FGFR3 signaling" R-HSA-5654706
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_graph("R-HSA-5654706", "proteoforms", False, graphs_path)


def test_create_graph_num_edges(glycolysis_graph):
    print(glycolysis_graph.edges)
    assert len(glycolysis_graph.edges) == 19


def test_create_graph_num_vertices(glycolysis_graph):
    print(glycolysis_graph.nodes)
    assert len(glycolysis_graph.nodes) == 13


# Test: Receive a Reaction with a direct participant EWAS input and a Simple entity output --> connects them
def test_connects_inputs_with_outputs(glycolysis_graph):
    # Input to output interactions for reaction R-HSA-163773:
    assert not ("P16118", "ADP") in glycolysis_graph.edges
    assert ("P16118", "P16118;00046:33") in glycolysis_graph.edges
    assert ("P16118;00046:33", "P16118") in glycolysis_graph.edges
    assert not ("ATP", "ADP") in glycolysis_graph.edges
    assert not ("ATP", "P16118") in glycolysis_graph.edges
    assert not ("ATP", "P16118;00046:33") in glycolysis_graph.edges

    # Input to output interactions for reaction R-HSA-163750
    assert not ("P16118", "Pi") in glycolysis_graph.edges
    assert not ("P16118;00046:33", "Pi") in glycolysis_graph.edges
    assert ("P16118;00046:33", "P16118") in glycolysis_graph.edges
    assert not ("H2O", "Pi") in glycolysis_graph.edges
    assert not ("H2O", "P16118") in glycolysis_graph.edges

    # Input to output interaction for reaction R-HSA-71802
    assert not ("Fru(6)P", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("Fru(6)P", "ADP") in glycolysis_graph.edges
    assert not ("ATP", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("ATP", "ADP") in glycolysis_graph.edges


def test_does_not_connect_input_to_output_when_same_molecule(traffic_graph):
    # In reaction R-HSA-5627071 there are same set of proteoforms as input and output, since they are in
    # different sub cellular locations. Nevertheless, there are no self interactions in our network
    assert not ("P13569", "P13569") in traffic_graph.edges
    assert not ("P17081", "P17081") in traffic_graph.edges
    assert not ("Q9HD26", "Q9HD26") in traffic_graph.edges
    assert not ("GTP", "GTP") in traffic_graph.edges


def test_connects_catalysts_with_outputs(glycolysis_graph):
    # Catalysts to output interactions for reaction R-HSA-163773:
    assert not ("P17612", "P16118") in glycolysis_graph.edges
    assert ("P17612", "P16118;00046:33") in glycolysis_graph.edges
    assert not ("P17612", "ADP") in glycolysis_graph.edges
    assert not ("P22694", "P16118") in glycolysis_graph.edges
    assert ("P22694", "P16118;00046:33") in glycolysis_graph.edges
    assert not ("P22694", "ADP") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-163750
    assert ("P30153", "P16118") in glycolysis_graph.edges
    assert not ("P30153", "Pi") in glycolysis_graph.edges
    assert ("P30154", "P16118") in glycolysis_graph.edges
    assert not ("P30154", "Pi") in glycolysis_graph.edges
    assert ("P67775", "P16118") in glycolysis_graph.edges
    assert not ("P67775", "Pi") in glycolysis_graph.edges
    assert ("P62714", "P16118") in glycolysis_graph.edges
    assert not ("P62714", "Pi") in glycolysis_graph.edges
    assert ("Q14738", "P16118") in glycolysis_graph.edges
    assert not ("Q14738", "Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-70262
    assert not ("P16118", "Fru(6)P") in glycolysis_graph.edges
    assert not ("P16118;00046:33", "Fru(6)P") in glycolysis_graph.edges
    assert not ("P16118", "Pi") in glycolysis_graph.edges
    assert not ("P16118;00046:33", "Pi") in glycolysis_graph.edges

    # Catalysts to output interactions for reaction R-HSA-71802
    assert not ("Q16875", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("Q16875", "ADP") in glycolysis_graph.edges
    assert not ("Q16877", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("Q16877", "ADP") in glycolysis_graph.edges
    assert not ("P16118", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("P16118", "ADP") in glycolysis_graph.edges
    assert not ("O60825", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert not ("O60825", "ADP") in glycolysis_graph.edges


def test_connects_regulators_with_outputs(traffic_graph):
    # Regulator to output interactions for reaction R-HSA-5627071:
    assert ("P17081", "P13569") in traffic_graph.edges
    assert ("Q9HD26", "P13569") in traffic_graph.edges
    assert not ("P13569", "P13569") in traffic_graph.edges
    assert not ("GTP", "P13569") in traffic_graph.edges

    assert ("P17081", "Q9HD26") in traffic_graph.edges
    assert not ("Q9HD26", "Q9HD26") in traffic_graph.edges
    assert ("P13569", "Q9HD26") in traffic_graph.edges
    assert not ("GTP", "Q9HD26") in traffic_graph.edges

    assert not ("P17081", "P17081") in traffic_graph.edges
    assert ("Q9HD26", "P17081") in traffic_graph.edges
    assert ("P13569", "P17081") in traffic_graph.edges
    assert not ("GTP", "P17081") in traffic_graph.edges

    assert not ("P17081", "GTP") in traffic_graph.edges
    assert not ("Q9HD26", "GTP") in traffic_graph.edges
    assert not ("P13569", "GTP") in traffic_graph.edges
    assert not ("GTP", "GTP") in traffic_graph.edges


def test_connects_regulators_with_outputs_2(signal_graph):
    # In reaction R-HSA-5654408
    assert ("Q9NP95", "O43320") or ("O43320", "Q9NP95") in signal_graph.edges
    assert ("O60258-1", "Q9NP95") in signal_graph.edges
    assert ("P22607-1;00048:577,00048:647,00048:648,00048:724,00048:760,00048:770", "Q9NP95") in signal_graph.edges
    assert not ("Q9GZV9;00164:178", "HS") in signal_graph.edges
    assert not ("Q9GZV9;00164:178", "Q9GZV9;00164:178") in signal_graph.edges
    assert ("Q9GZV9;00164:178",
            "Q8WU20;00048:196,00048:306,00048:349,00048:392,00048:436,00048:471") in signal_graph.edges
    assert not ("HS", "O60258-1") in signal_graph.edges


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


def test_connect_components_of_same_complex_2(signal_graph):
    assert not ("P01116-1;00115:179,01116:186", "GTP") in signal_graph.edges
    assert not ("P01111;00115:181,01116:186", "GTP") in signal_graph.edges
    assert not ("P01112;00115:181,00115:184,01116:186", "GTP") in signal_graph.edges
    assert not ("P01116-2;01116:186", "GTP") in signal_graph.edges

    assert ("P01112;00115:181,00115:184,01116:186", "P01111;00115:181,01116:186") in signal_graph.edges
    assert ("P01112;00115:181,00115:184,01116:186", "P01116-2;01116:186") in signal_graph.edges

    assert ("P01116-2;01116:186", "P01111;00115:181,01116:186") in signal_graph.edges
    assert ("P01116-2;01116:186", "P01112;00115:181,00115:184,01116:186") in signal_graph.edges

    assert ("P01111;00115:181,01116:186", "P01116-2;01116:186") in signal_graph.edges
    assert ("P01111;00115:181,01116:186", "P01112;00115:181,00115:184,01116:186") in signal_graph.edges

    assert ("Q07889", "P62993-1") in signal_graph.edges

    assert not ("HS", "P22607-2;00048:579,00048:649,00048:650,00048:726,00048:762,00048:772") in signal_graph.edges
    assert ("P22607-2;00048:579,00048:649,00048:650,00048:726,00048:762,00048:772", "P05230") in signal_graph.edges

    assert ("P55075-1", "P09038") in signal_graph.edges
    assert ("O60258-1", "Q9GZV9;00164:178") in signal_graph.edges
    assert ("Q9GZV9;00164:178", "O60258-1") in signal_graph.edges
