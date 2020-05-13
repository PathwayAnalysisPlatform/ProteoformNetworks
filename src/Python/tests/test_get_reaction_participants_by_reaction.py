import pytest
from pandas import DataFrame

from network_topology_queries import get_reaction_participants_by_reaction


@pytest.fixture(scope="session")
def get_o2_graph(tmpdir_factory):
    # Pathway "Hemoglobin A binds oxygen and releases protons and carbon dioxide"
    graphs_pah = tmpdir_factory.mktemp("tmpdir")
    return get_reaction_participants_by_reaction("R-HSA-1247668", "genes", True, graphs_pah)

@pytest.fixture(scope="session")
def get_o2_graph_no_small_molecules(tmpdir_factory):
    # Pathway "Hemoglobin A binds oxygen and releases protons and carbon dioxide"
    graphs_pah = tmpdir_factory.mktemp("tmpdir")
    return get_reaction_participants_by_reaction("R-HSA-1247668", "genes", False, graphs_pah)

@pytest.fixture(scope="session")
def map2ks_graph(tmpdir_factory):
    # Pathway "MAP2Ks phosphorylate MAPKs"
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return get_reaction_participants_by_reaction("R-HSA-5672973", "proteins", True, graphs_path)

@pytest.fixture(scope="session")
def map2ks_graph_no_small_molecules(tmpdir_factory):
    # Pathway "MAP2Ks phosphorylate MAPKs"
    graphs_path = tmpdir_factory.mktemp("tmpdir")
    return get_reaction_participants_by_reaction("R-HSA-5672973", "proteins", False, graphs_path)


def test_get_reaction_participants_genes_correct_entities(get_o2_graph):
    df = get_o2_graph
    assert type(df) == DataFrame
    assert len(df) == 10
    assert ((df['Entity'] == 'R-HSA-1237013') & (df['Id'] == 'HBB')).any()
    assert ((df['Entity'] == 'R-HSA-1237026') & (df['Id'] == 'HBA1')).any()
    assert ((df['Entity'] == 'R-ALL-71185') & (df['Id'] == 'heme')).any()
    assert ((df['Entity'] == 'R-ALL-29368') & (df['Id'] == 'O2')).any()


def test_get_reaction_participants_genes_correct_roles(get_o2_graph):
    df = get_o2_graph
    assert ((df['Id'] == 'HBB') & (df['Role'] == 'input')).any()
    assert ((df['Id'] == 'HBB') & (df['Role'] == 'output')).any()
    assert ((df['Id'] == 'HBA1') & (df['Role'] == 'input')).any()
    assert ((df['Id'] == 'HBA1') & (df['Role'] == 'output')).any()
    assert ((df['Id'] == 'CO2') & (df['Role'] == 'output')).any()


def test_get_reaction_participants_genes_fake_reaction_returns_empty_list():
    df = get_reaction_participants_by_reaction("fake_reaction", "genes", True)
    assert type(df) == DataFrame
    assert len(df) == 0


def test_get_reaction_participants_genes_no_small_molecules_correct_entities(get_o2_graph_no_small_molecules):
    df = get_o2_graph_no_small_molecules
    assert len(df) == 4
    assert ((df['Id'] == 'HBA1') & (df['Role'] == 'input')).any()
    assert ((df['Id'] == 'HBA1') & (df['Role'] == 'output')).any()
    assert ((df['Id'] == 'HBB') & (df['Role'] == 'input')).any()
    assert ((df['Id'] == 'HBB') & (df['Role'] == 'output')).any()
    assert not ((df['Entity'] == 'R-ALL-29368') & (df['Id'] == 'O2')).any()


def test_get_reaction_participants_proteins_correct_entities(map2ks_graph):
    df = map2ks_graph
    assert type(df) == DataFrame
    assert len(df) == 184
    assert ((df['Entity'] == 'R-HSA-109837') & (df['Id'] == 'Q02750')).any()
    assert ((df['Entity'] == 'R-HSA-5672632') & (df['Id'] == 'Q8WXI2')).any()
    assert ((df['Entity'] == 'R-HSA-140585') & (df['Id'] == 'P02679')).any()
    assert ((df['Entity'] == 'R-HSA-350713') & (df['Id'] == 'Q9Y490')).any()


def test_get_reaction_participants_proteins_correct_roles(map2ks_graph):
    df = map2ks_graph
    assert ((df['Id'] == 'Q02750') & (df['Role'] == 'catalystActivity')).any()
    assert ((df['Id'] == 'GTP') & (df['Role'] == 'catalystActivity')).any()
    assert ((df['Id'] == 'P31946') & (df['Role'] == 'input')).any()
    assert ((df['Id'] == 'P27448') & (df['Role'] == 'output')).any()
    assert ((df['Id'] == 'Ca2+') & (df['Role'] == 'output')).any()
    assert ((df['Id'] == 'Q02750') & (df['Role'] == 'regulatedBy')).any()


def test_get_reaction_participants_proteins_fake_reaction_returns_empty_list():
    df = get_reaction_participants_by_reaction("fake_reaction", "proteins", True)
    assert type(df) == DataFrame
    assert len(df) == 0


def test_get_reaction_participants_proteins_no_small_molecules_correct_entities(map2ks_graph_no_small_molecules):
    df = map2ks_graph_no_small_molecules
    assert len(df) == 162
    assert not (df['Id'] == 'Ca2+').any()
    assert not (df['Id'] == 'GTP').any()
    assert ((df['Id'] == 'Q02750') & (df['Role'] == 'catalystActivity')).any()
    assert ((df['Id'] == 'P31946') & (df['Role'] == 'input')).any()
    assert ((df['Id'] == 'P27448') & (df['Role'] == 'output')).any()
