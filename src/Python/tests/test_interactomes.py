
import pytest

# R-HSA-70171  Glycolysis
# Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
# Pathway "RHO GTPases regulate P13569 trafficking" R-HSA-5627083

# Reaction "Phosphorylation of PF2K-Pase by PKA catalytic subunit" R-HSA-163773
# Reaction "KRAS4B is phosphorylated on serine 181" R-HSA-9653503

#region Genes No small molecules

# Test get_participants()
# Test get_components()
# Test add_nodes()
# Test create_interaction_network()

#endregion

#region Genes With small molecules

# Test get_participants()
# Test get_components()
# Test add_nodes()
# Test create_interaction_network()

#endregion

#region Genes With reaction-specific ids for small molecules

# Test get_participants()
# Test get_components()
# Test add_nodes()
# Test create_interaction_network()

#endregion

#region Proteins No small molecules

# Test get_participants()
# Test get_components()
# Test add_nodes()
# Test create_interaction_network()

#endregion

#region Proteins With small molecules

# Test get_participants()
# Test get_components()
# Test add_nodes()
# Test create_interaction_network()

#endregion

#region Proteins With reaction-specific ids for small molecules

# Test get_participants()
# Test get_components()
# Test add_nodes()
# Test create_interaction_network()

#endregion

#region Proteoforms No small molecules

# Test get_participants()
# Test get_components()
# Test add_nodes()
# Test create_interaction_network()

#endregion

#region Proteoforms With small molecules

# Test get_participants()
# Test get_components()
# Test add_nodes()
# Test create_interaction_network()

#endregion

#region Proteoforms With reaction-specific ids for small molecules

# Test get_participants()
# Test get_components()
# Test add_nodes()
# Test create_interaction_network()

#endregion

@pytest.fixture(scope="session")
def interactome_gene_no_sm(tmpdir_factory):
    # Pathway "Regulation of glycolysis by fructose" R-HSA-9634600
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return create_interaction_network("R-HSA-9634600", "genes", True, graphs_path)

def test_connects_inputs_with_outputs(glycolysis_graph):
    # Input to output interactions for reaction R-HSA-163773:
    assert ("PFKFB1", "ADP") in glycolysis_graph.edges
    assert ("ATP", "ADP") in glycolysis_graph.edges
    assert ("ATP", "PFKFB1") in glycolysis_graph.edges

    # Input to output interactions for reaction R-HSA-163750
    assert ("PFKFB1", "Pi") in glycolysis_graph.edges
    assert ("H2O", "Pi") in glycolysis_graph.edges
    assert ("H2O", "PFKFB1") in glycolysis_graph.edges

    # Input to output interaction for reaction R-HSA-71802
    assert ("Fru(6)P", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("Fru(6)P", "ADP") in glycolysis_graph.edges
    assert ("ATP", "D-Fructose 2,6-bisphosphate") in glycolysis_graph.edges
    assert ("ATP", "ADP") in glycolysis_graph.edges
