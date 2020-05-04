import re

from Python.config import LEVELS
from Python.lib.graph_database import get_query_result


def get_pathways():
    query = "MATCH (p:Pathway{speciesName:\"Homo sapiens\"})\nRETURN p.stId as stId, p.displayName as displayName"
    return get_query_result(query)


def get_reactions_by_pathway(pathway):
    query = f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:Reaction{{speciesName:'Homo sapiens'}}) RETURN rle.stId as reaction"
    return get_query_result(query)


def get_reactions_and_participants_by_pathway(pathway, showSmallMolecules=True, level="proteins"):
    """
    Notice that records are sorted by reaction for easy traversal later.
    :param pathway: string id of a pathway
    :return: pandas dataframe with the records
    """
    if level not in LEVELS:
        raise Exception

    print(f"\n\nQuerying {level} participants of pathway {pathway}...\n\n")

    query = f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:ReactionLikeEvent{{speciesName:'Homo sapiens'}})"
    query += "\nWITH rle"
    query += "\nMATCH p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity)"
    query += "\nWHERE last(labels(pe)) IN [\"EntityWithAccessionedSequence\""
    if showSmallMolecules:
        query += ", \"SimpleEntity\""
    query += "]"
    query += "\nOPTIONAL MATCH (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)"
    query += "\nRETURN DISTINCT rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name, last(labels(pe)) as Type,"
    query += "\nCASE WHEN last(labels(pe)) = \"SimpleEntity\" THEN pe.displayName "
    if level == "genes":
        query += " WHEN last(labels(pe)) = \"EntityWithAccessionedSequence\" THEN head(re.geneName) "
    else:
        query += " WHEN last(labels(pe)) = \"EntityWithAccessionedSequence\" THEN re.identifier "
    query += " ELSE re.identifier END as Id,"
    query += " rd.displayName AS Database, head([scores IN relationships(p) | type(scores)]) as Role"
    query += "\nORDER BY Reaction, Role, Type"
    print(query)

    df = get_query_result(query)

    df['Id'] = df.apply(lambda x: re.sub(r'\s*\[[\w\s]*\]\s*', '', x.Id) if x.Type == 'SimpleEntity' else x.Id, axis=1)

    return df
