import pytest

from config import no_sm, genes
from lib.networks import get_or_create_pathway_interaction_network


def test_create_pathway_network_has_components_and_participants():
    g = get_or_create_pathway_interaction_network("R-HSA-9634600", genes, no_sm)
    assert g.nodes > 0
