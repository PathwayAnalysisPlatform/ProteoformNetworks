from Python.lib.graph_database import get_query_result


def get_reactions_by_pathway(pathway):
    query = f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:Reaction{{speciesName:'Homo sapiens'}}) RETURN rle.stId"
    return get_query_result(query)

def get_reactions_and_participants_by_pathway(pathway):
    query = f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:Reaction{{speciesName:'Homo sapiens'}})"
    query += "\nWITH rle"
    query += "\nMATCH p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity)"
    query += "\nWHERE last(labels(pe)) IN [\"EntityWithAccessionedSequence\", \"SimpleEntity\"]"
    query += "\nOPTIONAL MATCH (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)"
    query += "\nRETURN DISTINCT rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name, last(labels(pe)) as Type, re.identifier as Id, rd.displayName AS Database, head([x IN relationships(p) | type(x)]) as Role"
    query += "\nORDER BY Role, Type"

    return get_query_result(query)