import pytest

# Pathway "Regulation of glycolysis by fructose" R-HSA-9634600

#region Genes No small molecules
from config import genes
from lib.graph_database_access import get_participants_by_pathway


@pytest.fixture(scope="session")
def glycolysis_get_participants(tmpdir_factory):
    return get_participants_by_pathway(genes, "R-HSA-9634600")

@pytest.fixture(scope="session")
def glycolysis_genes_no_sm(glycolysis_get_participants):
    participants = get_participants_by_pathway(genes, "R-HSA-9634600")
    components = get_components(genes)
    return

# Test get_participants()
def test_genes_no_sm_get_participants():
    participants = get_participants_by_pathway(genes, "R-HSA-9634600")


# Test get_components()
# Test add_nodes()
# Test create_interaction_network()
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