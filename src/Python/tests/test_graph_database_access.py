import pandas as pd
import pytest

from config import proteoforms, proteins, genes, sm
from lib.graph_database_access import fix_neo4j_values, get_low_level_pathways, \
    get_reactions_by_pathway, get_pathways, get_reactions, get_complexes, get_complex_components_by_complex, \
    make_proteoform_string, get_participants_by_pathway, get_components_by_pathway, get_complexes_by_pathway


@pytest.fixture(scope="session")
def glycolysis_participants_genes():
    return get_participants_by_pathway("R-HSA-70171", genes)


@pytest.fixture(scope="session")
def glycolysis_participants_proteins():
    return get_participants_by_pathway("R-HSA-70171", proteins)


@pytest.fixture(scope="session")
def glycolysis_components_proteins():
    return get_components_by_pathway("R-HSA-70171", proteins)


@pytest.fixture(scope="session")
def glycolysis_participants_proteoforms():
    return get_participants_by_pathway("R-HSA-70171", proteoforms)


@pytest.fixture(scope="session")
def glycolysis_participants_sm():
    return get_participants_by_pathway("R-HSA-70171", sm)


@pytest.fixture(scope="session")
def ras_processing_participants_genes():
    return get_participants_by_pathway("R-HSA-9648002", genes)


@pytest.fixture(scope="session")
def ras_processing_participants_proteins():
    return get_participants_by_pathway("R-HSA-9648002", proteins)


@pytest.fixture(scope="session")
def ras_processing_participants_proteoforms():
    return get_participants_by_pathway("R-HSA-9648002", proteoforms)


@pytest.fixture(scope="session")
def ras_processing_participants_sm():
    return get_participants_by_pathway("R-HSA-9648002", sm)


def test_participant_records_columns_genes(ras_processing_participants_genes):
    df = ras_processing_participants_genes
    assert "Pathway" in df.columns
    assert "Reaction" in df.columns
    assert "Entity" in df.columns
    assert "Name" in df.columns
    assert "Type" in df.columns
    assert "Id" in df.columns
    assert "Database" in df.columns
    assert "Role" in df.columns


def test_participant_records_columns_proteins(ras_processing_participants_proteins):
    df = ras_processing_participants_proteins
    assert "Pathway" in df.columns
    assert "Reaction" in df.columns
    assert "Entity" in df.columns
    assert "Name" in df.columns
    assert "Type" in df.columns
    assert "Id" in df.columns
    assert "Database" in df.columns
    assert "Role" in df.columns


def test_participant_records_columns_proteoforms(ras_processing_participants_proteoforms):
    df = ras_processing_participants_proteoforms
    assert "Pathway" in df.columns
    assert "Reaction" in df.columns
    assert "Entity" in df.columns
    assert "Name" in df.columns
    assert "Type" in df.columns
    assert "Id" in df.columns
    assert "Database" in df.columns
    assert "Role" in df.columns


def test_participant_records_columns_sm(ras_processing_participants_sm):
    df = ras_processing_participants_sm
    assert "Pathway" in df.columns
    assert "Reaction" in df.columns
    assert "Entity" in df.columns
    assert "Name" in df.columns
    assert "Type" in df.columns
    assert "Id" in df.columns
    assert "Database" in df.columns
    assert "Role" in df.columns


def test_fix_neo4j_values_empty_dataframe():
    df = pd.DataFrame()
    result = fix_neo4j_values(df, level="proteins")
    assert len(result) == 0
    assert type(result) == pd.DataFrame
    assert len(result.columns) == 0


def test_pathway_not_exists_returns_empty_dataframe():
    result = get_participants_by_pathway("blabla", genes)
    assert len(result) == 0
    assert type(result) == pd.DataFrame
    assert len(result.columns) == 0


# Test: Cypher query to get all pathways gets all of them correctly
def test_get_pathways():
    pathways = get_low_level_pathways()
    assert len(pathways) == 1657
    assert type(pathways) == pd.DataFrame
    assert ((pathways['stId'] == "R-HSA-110056") & (pathways['displayName'] == "MAPK3 (ERK1) activation")).any()
    assert ((pathways['stId'] == "R-HSA-69200") & (pathways[
                                                       'displayName'] == "Phosphorylation of proteins involved in G1/S transition by active Cyclin E:Cdk2 complexes")).any()
    assert ((pathways['stId'] == "R-HSA-6782135") & (pathways['displayName'] == "Dual incision in TC-NER")).any()
    assert not ((pathways['stId'] == 'R-HSA-9612973') & (pathways['displayName'] == 'Autophagy')).any()
    assert not ((pathways['stId'] == 'R-HSA-1640170') & (pathways['displayName'] == 'Cell Cycle')).any()
    assert not ((pathways['stId'] == 'R-HSA-70171') & (pathways['displayName'] == 'Glycolysis')).any()


# Test: Query to get reactions of pathway returns the correct list of reactions
def test_get_reactions_of_pathway():
    reactions = get_reactions_by_pathway("R-HSA-170822")
    assert len(reactions) == 5
    assert (reactions['reaction'] == "R-HSA-170824").any()
    assert (reactions['reaction'] == "R-HSA-170796").any()
    assert (reactions['reaction'] == "R-HSA-170825").any()


# Test: Query to get participants of a pathway gets all participant reactions.
def test_query_for_pathway_participants_has_all_reactions(glycolysis_participants_proteins):
    df = glycolysis_participants_proteins
    assert len(df['Reaction'].unique()) == 15
    assert (df['Reaction'] == "R-HSA-8955794").any()
    assert (df['Reaction'] == "R-HSA-6799604").any()
    assert (df['Reaction'] == "R-HSA-70467").any()


# Test: Query to get participants of a pathway returns all small molecule participants
def test_query_for_pathway_participans_returns_all_simple_molecules(glycolysis_participants_sm):
    df = glycolysis_participants_sm
    print(df.loc[df['Type'] == 'SimpleEntity'])
    assert len(df.loc[df['Type'] == 'SimpleEntity']['Entity'].unique()) == 32


# Test: Query to get participants of a pathway returns all gene participants
def test_query_for_pathway_participans_returns_all_ewas(glycolysis_participants_proteins):
    df = glycolysis_participants_proteins
    print(df.loc[df['Type'] == 'EntityWithAccessionedSequence'])
    assert len(df.loc[df['Type'] == 'EntityWithAccessionedSequence']['Entity'].unique()) == 31


@pytest.fixture(scope="session")
def participant_proteins_Erythrocytes_take_up_carbon_dioxide_and_release_oxygen():
    return get_participants_by_pathway("R-HSA-1237044", proteins)


# Test: Query for participants of a pathway returns all the gene and small molecule participants decomposing complexes
# Pathway R-HSA-1237044 hast complex participants in its reactions, some complexes also have complex components.
# Checks if the component molecules of the complex components are in the result
def test_query_for_pathway_participants_decomposes_complexes(
        participant_proteins_Erythrocytes_take_up_carbon_dioxide_and_release_oxygen):
    df = participant_proteins_Erythrocytes_take_up_carbon_dioxide_and_release_oxygen
    # Reaction R-HSA-1237325 in the pathway has complex participants: "R-HSA-1237320"
    # Complex R-HSA-1237320 has 6 participant molecules:
    assert (df['Id'] == "P22748").any()
    assert (df['Id'] == "P29972").any()
    assert (df['Id'] == "P02730").any()
    assert (df['Id'] == "P68871").any()
    assert (df['Id'] == "P69905").any()


# Test: Query for participants of a pathway returns all the gene and small molecule participants decomposing sets
# Pathway R-HSA-70171 has reactions with EntitySets as participants, like reaction R-HSA-70420 with set R-HSA-450097
def test_query_for_pathway_participants_decomposes_sets(
        participant_proteins_Erythrocytes_take_up_carbon_dioxide_and_release_oxygen):
    df = participant_proteins_Erythrocytes_take_up_carbon_dioxide_and_release_oxygen
    # DefinedSet R-HSA-450097 has 4 members
    assert (df['Id'] == "P00915").any()
    assert (df['Id'] == "P00918").any()


# Test: If small molecules disabled, cypher query returns only the gene direct participants
def test_query_for_pathway_participants_disable_small_molecules(glycolysis_participants_proteins):
    df = glycolysis_participants_proteins
    assert len(df) == 31  # Only EWAS
    # Has EWAS participants
    assert (df['Entity'] == "R-HSA-5696062").any()
    assert (df['Entity'] == "R-HSA-6798334").any()
    assert (df['Entity'] == "R-HSA-6799597").any()
    # Does not have Small molecule participants
    assert (df['Entity'] != "R-ALL-449701").all()
    assert (df['Entity'] != "R-ALL-70113").all()
    assert (df['Entity'] != "R-ALL-217314").all()


# Test: If small molecules disabled, cypher query returns only the gene complex decomposed participants
def test_query_for_pathway_participants_complexes_show_only_ewas(glycolysis_participants_proteins):
    df = glycolysis_participants_proteins
    # The pathway has the reaction R-HSA-5696021 which has the complex R-HSA-5696043 as participant
    # THe components of the complex are 2, one EWAS and one small molecule:
    assert (df['Entity'] == 'R-HSA-5696062').any()  # EWAS ADPGK
    assert (df['Entity'] != 'R-ALL-217314').all()  # Small molecule Mg2+


# Test: Get pathway participants as genes
def test_query_for_pathway_participants_as_genes(glycolysis_participants_genes):
    df = glycolysis_participants_genes
    assert len(df) == 31

    # Participant protein Q9BRR6 should be in the result as a gene name: ADPGK not as UniProt accession
    assert (df['Id'] == 'ADPGK').any()
    assert (df['Id'] != 'Q9BRR6').all()

    # Simmilar for the next proteins
    assert (df['Id'] == 'BPGM').any()
    assert (df['Id'] != 'P07738').all()

    assert (df['Id'] == 'BPGM').any()
    assert (df['Id'] != 'P07738').all()

    assert (df['Id'] == 'PKLR').any()
    assert (df['Id'] != 'P07738').all()


def test_query_for_pathway_participants_as_genes_trims_gene_id(glycolysis_participants_genes):
    df = glycolysis_participants_genes
    assert ((df['Entity'] == 'R-HSA-70097') & (df['Id'] == "PKLR")).any()
    assert ((df['Entity'] == 'R-HSA-450658') & (df['Id'] == "PKM")).any()
    assert not ((df['Entity'] == 'R-HSA-450658') & (df['Id'] == "PKM-2 [cytosol]")).any()
    assert not ((df['Entity'] == 'R-HSA-450658') & (df['Id'] == "P62993")).any()
    print(df.loc[(df['Entity'] == 'R-HSA-70097') & (df['Id'] == "PKLR")])
    assert not ((df['Entity'] == 'R-HSA-70097') & (df['Id'] == "PKLR-1 [cytosol]")).any()
    assert not ((df['Entity'] == 'R-HSA-211388') & (df['Id'] == "PKLR-2")).any()


def test_query_for_pathway_participants_replaces_small_molecule_names(glycolysis_participants_sm):
    df = glycolysis_participants_sm
    assert ((df['Entity'] == 'R-ALL-29370') & (df['Id'] == 'sm_ADP')).any()
    assert ((df['Entity'] == 'R-ALL-29926') & (df['Id'] == 'sm_Mg2+')).any()
    assert ((df['Entity'] == 'R-ALL-70106') & (df['Id'] == 'sm_H+')).any()


# Test: Get pathway participants as proteins
def test_query_for_pathway_participants_as_proteins(glycolysis_participants_proteins):
    df = glycolysis_participants_proteins
    assert not ((df['Entity'] == 'R-HSA-5696062') & (df['Id'] == 'ADPGK')).any()
    assert ((df['Entity'] == 'R-HSA-5696062') & (df['Id'] == 'Q9BRR6')).any()

    assert not ((df['Entity'] == "R-HSA-6798334") & (df['Id'] == 'BPGM')).any()
    assert ((df['Entity'] == "R-HSA-6798334") & (df['Id'] == 'P07738')).any()

    assert not ((df['Entity'] == 'R-HSA-70412') & (df['Id'] == 'HK3')).any()
    assert ((df['Entity'] == 'R-HSA-70412') & (df['Id'] == 'P52790')).any()


def test_query_for_pathway_participants_as_proteins_returns_genes(glycolysis_participants_proteins):
    df = glycolysis_participants_proteins
    assert 'PrevId' in df.columns
    assert ((df['Id'] == 'Q9BRR6') & (df['PrevId'] == 'ADPGK')).any()
    assert ((df['Id'] == 'P07738') & (df['PrevId'] == 'BPGM')).any()
    assert not ((df['Id'] == 'Mg2+')).any()
    assert not ((df['Id'] == 'ADP')).any()


def test_query_for_pathway_participants_as_proteins_implicit_parameter(glycolysis_participants_proteins):
    df = glycolysis_participants_proteins
    assert not ((df['Entity'] == 'R-HSA-5696062') & (df['Id'] == 'ADPGK')).any()
    assert ((df['Entity'] == 'R-HSA-5696062') & (df['Id'] == 'Q9BRR6')).any()

    assert not ((df['Entity'] == "R-HSA-6798334") & (df['Id'] == 'BPGM')).any()
    assert ((df['Entity'] == "R-HSA-6798334") & (df['Id'] == 'P07738')).any()

    assert not ((df['Entity'] == 'R-HSA-70412') & (df['Id'] == 'HK3')).any()
    assert ((df['Entity'] == 'R-HSA-70412') & (df['Id'] == 'P52790')).any()


def test_query_for_pathway_participants_as_proteins_complex_should_not_be_in_records(glycolysis_participants_proteins):
    df = glycolysis_participants_proteins
    assert not (df["Entity"] == "R-HSA-5696043").any()
    assert not (df["Name"] == "BPGM dimer [cytosol]").any()
    assert not (df["Entity"] == "R-HSA-6799598").any()


def test_get_pathways():
    df = get_pathways()

    assert type(df) == pd.DataFrame
    assert len(df) == 2516
    assert "stId" in df.columns
    assert "displayName" in df.columns
    assert ((df['stId'] == "R-HSA-9612973") & (df["displayName"] == "Autophagy")).any()
    assert ((df['stId'] == "R-HSA-3000480") & (df["displayName"] == "Scavenging by Class A Receptors")).any()


def test_get_reactions():
    result = get_reactions()
    assert len(result) == 13661


def test_get_complexes():
    result = get_complexes()
    assert len(result) == 13362


def test_get_complex_components_genes_returns_components():
    df = get_complex_components_by_complex("R-HSA-983126", genes)
    assert type(df) == pd.DataFrame
    assert len(df) == 264
    assert ((df['Entity'] == 'R-HSA-141412') & (df['Id'] == 'CDC20')).any()
    assert ((df['Entity'] == 'R-HSA-174229') & (df['Id'] == 'ANAPC2')).any()
    assert ((df['Entity'] == 'R-HSA-939239') & (df['Id'] == 'UBC')).any()


def test_get_complex_components_genes_returns_components_2():
    df = get_complex_components_by_complex("R-HSA-2168879", genes)
    assert df['Id'].nunique() == 4
    assert (df['Id'] == 'HP').any()
    assert (df['Id'] == 'HBA1').any()
    assert (df['Id'] == 'HBB').any()
    assert (df['Id'] == 'CD163').any()


def test_get_complex_components_genes_with_non_existent_complex_returns_empty_list():
    df = get_complex_components_by_complex("fake_complex", genes)
    assert len(df) == 0


def test_get_complex_components_genes_without_small_molecules():
    df = get_complex_components_by_complex("R-HSA-2168879", genes)
    assert len(df) == 5
    assert not ((df['Entity'] == 'R-ALL-352327') & (df['Id'] == 'O2')).any()
    assert not ((df['Entity'] == 'R-ALL-917877') & (df['Id'] == 'heme')).any()
    assert ((df['Entity'] == 'R-HSA-2168862') & (df['Id'] == 'HBA1')).any()


def test_get_complex_components_proteins_returns_components():
    df = get_complex_components_by_complex("R-HSA-983126", proteins)
    assert type(df) == pd.DataFrame
    assert len(df) == 264
    assert not ((df['Entity'] == 'R-HSA-141412') & (df['Id'] == 'CDC20')).any()
    assert ((df['Entity'] == 'R-HSA-141412') & (df['Id'] == 'Q12834')).any()
    assert not ((df['Entity'] == 'R-HSA-174229') & (df['Id'] == 'ANAPC2')).any()
    assert ((df['Entity'] == 'R-HSA-174229') & (df['Id'] == 'Q9UJX6')).any()
    assert not ((df['Entity'] == 'R-HSA-939239') & (df['Id'] == 'UBC')).any()
    assert ((df['Entity'] == 'R-HSA-939239') & (df['Id'] == 'P0CG48')).any()


def test_get_complex_components_proteins_returns_components_2():
    df = get_complex_components_by_complex("R-HSA-2168879", proteins)
    assert df['Id'].nunique() == 4
    assert not ((df['Entity'] == 'R-ALL-352327') & (df['Id'] == 'O2')).any()
    assert not ((df['Entity'] == 'R-ALL-917877') & (df['Id'] == 'heme')).any()
    assert not ((df['Entity'] == 'R-HSA-2168862') & (df['Id'] == 'HBA1')).any()
    assert (df['Id'] == 'P00738').any()
    assert (df['Id'] == 'Q86VB7').any()
    assert (df['Id'] == 'P69905').any()
    assert (df['Id'] == 'P68871').any()


def test_get_complex_components_proteins_with_non_existent_complex_returns_empty_list():
    df = get_complex_components_by_complex("fake_complex", proteins)
    assert len(df) == 0


def test_get_complex_components_proteins_without_small_molecules():
    df = get_complex_components_by_complex("R-HSA-2168879", proteins)
    assert len(df) == 5
    assert not ((df['Entity'] == 'R-ALL-352327') & (df['Id'] == 'O2')).any()
    assert not ((df['Entity'] == 'R-ALL-917877') & (df['Id'] == 'heme')).any()
    assert ((df['Entity'] == 'R-HSA-2168862') & (df['Id'] == 'P69905')).any()
    assert ((df['Entity'] == 'R-HSA-2168872') & (df['Id'] == 'P68871')).any()


def test_get_complex_components_proteoforms_returns_components():
    df = get_complex_components_by_complex("R-HSA-983126", proteoforms)
    assert type(df) == pd.DataFrame
    assert len(df) == 264
    assert ((df['Entity'] == 'R-HSA-141412') & (df['Id'] == 'Q12834;')).any()
    assert ((df['Entity'] == 'R-HSA-174229') & (df['Id'] == 'Q9UJX6;')).any()
    assert ((df['Entity'] == 'R-HSA-939239') & (df['Id'] == 'P0CG48;')).any()


def test_get_complex_components_proteoforms_returns_components_2():
    df = get_complex_components_by_complex("R-HSA-2168879", proteoforms)
    assert df["Id"].nunique() == 5
    assert not (df['Id'] == 'O2').any()
    assert not (df['Id'] == 'heme').any()
    assert (df['Id'] == 'P00738;00034:52,00034:111,00034:149').any()
    assert (df['Id'] == 'P00738;00034:266,00034:309,00034:351').any()
    assert (df['Id'] == 'P69905;').any()
    assert (df['Id'] == 'P68871;').any()
    assert (df['Id'] == 'Q86VB7;').any()


def test_get_complex_components_proteoforms_with_non_existent_complex_returns_empty_list():
    df = get_complex_components_by_complex("fake_complex", proteoforms)
    assert len(df) == 0


def test_get_complex_components_proteoforms_without_small_molecules():
    df = get_complex_components_by_complex("R-HSA-2168879", proteoforms)
    assert len(df) == 5
    assert not ((df['Entity'] == 'R-ALL-352327') & (df['Id'] == 'O2')).any()
    assert not ((df['Entity'] == 'R-ALL-917877') & (df['Id'] == 'heme')).any()
    assert ((df['Entity'] == 'R-HSA-2168862') & (df['Id'] == 'P69905;')).any()


def test_make_proteoform_string_multiple_values():
    values = ["P40189-1", "00048:759", "00048:767", "00048:814", "00048:905", "00048:915"]
    ans = make_proteoform_string(values)
    assert ans == "P40189-1;00048:759,00048:767,00048:814,00048:905,00048:915"


def test_make_proteoform_string_single_value():
    values = ["P40189-1"]
    ans = make_proteoform_string(values)
    assert ans == "P40189-1;"


def test_make_proteoform_string_one_ptm():
    values = ["P40189-1", "00048:759"]
    ans = make_proteoform_string(values)
    assert ans == "P40189-1;00048:759"


def test_make_proteoform_string_string():
    value = "P40189-1"
    ans = make_proteoform_string(value)
    assert ans == "P40189-1;"


# The selected pathway is R-HSA-9634600 Regulation of glycolysis by fructose 2,6-bisphosphate metabolism
# It has 6 participant complexes

# Reaction R-HSA-163750 Dephosphorylation of phosphoPFKFB1 by PP2A complex
def test_get_participants_by_pathway_returns_components_of_reaction_R_HSA_163750():
    df = get_participants_by_pathway("R-HSA-9634600", proteins)
    assert ((df['Id'] == "P67775") & (df['Role'] == "catalystActivity")).any()
    assert ((df['Id'] == "P30153") & (df['Role'] == "catalystActivity")).any()
    assert ((df['Id'] == "P16118") & (df['Role'] == "input")).any()
    assert ((df['Id'] == "P16118") & (df['Role'] == "output")).any()

def test_get_complexes_by_pathway1():
    df = get_complexes_by_pathway("R-HSA-9634600")
    assert df['Complex'].nunique() == 3
    assert (df['Complex'] == "R-HSA-165961").any()
    assert (df['Complex'] == "R-HSA-163745").any()
    assert (df['Complex'] == "R-HSA-71786").any()

def test_get_participants_by_pathway_returns_components_of_reaction_R_HSA_163750():
    df = get_components_by_pathway("R-HSA-9634600", proteins)

    # Check it has all the complexes
    assert df['Complex'].nunique() == 3
    assert df[df['Complex'] == "R-HSA-163745"]['Id'].nunique() == 1
    assert df[df['Complex'] == "R-HSA-165961"]['Id'].nunique() == 5
    assert df[df['Complex'] == "R-HSA-71786"]['Id'].nunique() == 1

    assert ((df['Complex'] == "R-HSA-163745") & (df['Id'] == "P16118")).any()
    assert ((df['Complex'] == "R-HSA-71786") & (df['Id'] == "P16118")).any()
    assert ((df['Complex'] == "R-HSA-165961") & (df['Id'] == "P67775")).any()
    assert ((df['Complex'] == "R-HSA-165961") & (df['Id'] == "P62714")).any()
    assert ((df['Complex'] == "R-HSA-165961") & (df['Id'] == "P30153")).any()
    assert ((df['Complex'] == "R-HSA-165961") & (df['Id'] == "P30154")).any()
    assert ((df['Complex'] == "R-HSA-165961") & (df['Id'] == "Q14738")).any()
