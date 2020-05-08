import pytest
import pandas as pd

from lib.graph_database import get_query_result


def test_get_query_result():
    query = "MATCH (p:Pathway{stId:\"R-HSA-9634600\"}) RETURN p.displayName as Name"
    result = get_query_result(query)
    assert len(result) == 1
    assert type(result) == pd.DataFrame
    assert "Name" in result.columns
    assert result.iloc[0]['Name'] == "Regulation of glycolysis by fructose 2,6-bisphosphate metabolism"


def test_get_query_result_empty_result():
    query = "MATCH (p:Pathway{stId:\"something_wrong\"}) RETURN p.displayName as Name"
    result = get_query_result(query)
    assert len(result) == 0
    assert type(result) == pd.DataFrame
    assert "Name" not in result.columns
    assert len(result.columns) == 0