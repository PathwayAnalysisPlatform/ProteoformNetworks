import config
from config import no_sm, with_sm, with_unique_sm, sm, LEVELS
from lib.graph_database_access import get_pathway_name, get_participants_by_pathway, get_components_by_pathway, \
    get_pathways, get_participants, get_components


def get_or_create_interaction_network(level, method, participants, components, out_path="", label=""):
    """
    Returns a networkx graph instance of the selected interaction network.
    Tries to read from a json file with the network. If the file does not exists it creates it.

    :param level: genes, proteins or proteoforms. This attribute is just to set it as graph property.
    :param method: "no_sm", "with_sm" or "with_unique_sm". This is just to set is as graph property.
    :param participants: pandas dataframe with the reaction participants
    :param components: pandas dataframe with the complex components
    :param out_path: path to directory to store the json file
    :param label: Any text to distinguish the graph file, ex. Pathway name
    :return: The networkx interaction network
    """

    filename = get_json_filename(level, method, out_path, label)

    if not Path(filename).exists():
        create_interaction_network(
            level, method, participants, components, out_path, label)
    g = read_graph(filename)

    if not 'Articulation Points' in g.graph:
        set_articulation_points(g)
        set_num_articulation_points(g)
        set_bridges(g)
        set_num_bridges(g)
        update_json_file(g, level, method, out_path)

    return g


def get_interactomes(input_data_path, output_networks_path):
    """ Get or create interactomes: With only accessioned sequence nodes, with small molecule nodes and with reaction-unique-small molecules"""
    participant_records = {l: get_participants(
        l, input_data_path) for l in [*LEVELS, sm]}
    components_records = {l: get_components(
        l, input_data_path) for l in [*LEVELS, sm]}

    interactomes_no_sm = {
        l: get_or_create_interaction_network(l, no_sm, participant_records, components_records, data_path) for l in
        LEVELS}
    interactomes_with_sm = {
        l: get_or_create_interaction_network(l, with_sm, participant_records, components_records, data_path) for l in
        LEVELS}
    interactomes_with_unique_sm = {
        l: get_or_create_interaction_network(
            l, with_unique_sm, participant_records, components_records, data_path)
        for l in LEVELS}
    return interactomes_no_sm, interactomes_with_sm, interactomes_with_unique_sm
