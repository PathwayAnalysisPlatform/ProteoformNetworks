from Python.lib.graph_database import get_query_result

def get_pathways():
    query = "MATCH (p:Pathway{speciesName:\"Homo sapiens\"})\nRETURN p.stId as stId, p.displayName as displayName"
    return get_query_result(query)

def get_reactions_by_pathway(pathway):
    query = f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:Reaction{{speciesName:'Homo sapiens'}}) RETURN rle.stId as reaction"
    return get_query_result(query)

def get_reactions_and_participants_by_pathway(pathway):
    """
    Notice that records are sorted by reaction for easy traversal later.
    :param pathway: string id of a pathway
    :return: pandas dataframe with the records
    """
    query = f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:ReactionLikeEvent{{speciesName:'Homo sapiens'}})"
    query += "\nWITH rle"
    query += "\nMATCH p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity)"
    query += "\nWHERE last(labels(pe)) IN [\"EntityWithAccessionedSequence\", \"SimpleEntity\"]"
    query += "\nOPTIONAL MATCH (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)"
    query += "\nRETURN DISTINCT rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name, last(labels(pe)) as Type, re.identifier as Id, rd.displayName AS Database, head([scores IN relationships(p) | type(scores)]) as Role"
    query += "\nORDER BY Reaction, Role, Type"

    return get_query_result(query)