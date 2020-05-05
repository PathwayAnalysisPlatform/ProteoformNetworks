import re

from Python.config import LEVELS
from Python.lib.graph_database import get_query_result


def get_pathways():
    query = "MATCH (p:Pathway{speciesName:\"Homo sapiens\"})\nRETURN p.stId as stId, p.displayName as displayName"
    return get_query_result(query)


def get_reactions_by_pathway(pathway):
    query = f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:Reaction{{speciesName:'Homo sapiens'}}) RETURN rle.stId as reaction"
    return get_query_result(query)


def make_proteoform_string(value):
    if type(value) == str:
        return value
    proteoform = value[0]
    if len(value) == 2:
        proteoform = ";".join(value)
    elif len(value) > 2:
        isoform = value[0] + ";"
        ptms = ",".join(value[1:])
        proteoform = isoform + ptms
    return proteoform


def get_reactions_and_participants_by_pathway(pathway, showSmallMolecules=True, level="proteins"):
    """
    Notice that records are sorted by reaction for easy traversal later.
    :param pathway: string id of a pathway
    :return: pandas dataframe with the records
    """
    if level not in LEVELS:
        raise Exception

    print(f"\n\nQuerying {level} participants of pathway {pathway}...\n\n")
    query = ""
    if level in ["genes", "proteins"]:
        query += f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:ReactionLikeEvent{{speciesName:'Homo sapiens'}})"
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
    else:
        query += f"""
// Get pathway proteoform participants
MATCH (p:Pathway{{stId:"{pathway}"}})-[:hasEvent]->(rle:ReactionLikeEvent{{speciesName:'Homo sapiens'}})
WITH rle
MATCH p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity)
WHERE last(labels(pe)) IN ["EntityWithAccessionedSequence" """
        if showSmallMolecules:
            query += ", \"SimpleEntity\""
        query += """]
OPTIONAL MATCH (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)
WITH DISTINCT rle, pe, re, 
              last(labels(pe)) as Type,
              rd.displayName AS Database, 
              head([x IN relationships(p) | type(x)]) as Role
ORDER BY rle.stId, Role, Type
OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
WITH DISTINCT rle.stId as Reaction, 
			  pe.stId as Entity, 
              pe.displayName as Name,
              Type,
              CASE 
              	WHEN last(labels(pe)) = "SimpleEntity" THEN pe.displayName  
                WHEN last(labels(pe)) = "EntityWithAccessionedSequence" THEN 
                	CASE 
              			WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier 
                		ELSE re.identifier
              		END
              END  as Id,
              mod.identifier as ptm_type,
              tm.coordinate as ptm_coordinate,
              Database,
              Role
ORDER BY ptm_type, ptm_coordinate
WITH DISTINCT Reaction, Entity, Name, Type, Id,
                COLLECT(
              		ptm_type + ":" + CASE WHEN ptm_coordinate IS NOT NULL THEN ptm_coordinate ELSE "null" END
              	) AS ptms,
                Database, Role
RETURN DISTINCT Reaction, Entity, Name, Type, CASE WHEN Type = "SimpleEntity" THEN Id ELSE (Id+ptms) END as Id,
                Database, Role
ORDER BY Reaction, Role, Id"""

    print(query)

    df = get_query_result(query)

    df['Id'] = df.apply(lambda x: re.sub(r'\s*\[[\w\s]*\]\s*', '', x.Id) if x.Type == 'SimpleEntity' else x.Id, axis=1)
    if level == "proteoforms":
        df['Id'] = df['Id'].apply(make_proteoform_string)

    return df
