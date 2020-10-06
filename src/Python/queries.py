import re

import pandas as pd

from config import LEVELS
from lib.graph_database import get_query_result

QUERIES_PARTICIPANTS = {
    LEVELS[0]: """
    MATCH (pw:Pathway{speciesName:'Homo sapiens'})-[:hasEvent]->(rle:ReactionLikeEvent{speciesName:'Homo sapiens'}),
          p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'}),
          (pe)-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name, 
                    last(labels(pe)) as Type, head(re.geneName) as Id, re.databaseName AS Database, 
                    head([scores IN relationships(p) | type(scores)]) as Role
    ORDER BY Pathway, Reaction, Role, Type
    """,
    LEVELS[1]: """
    MATCH (pw:Pathway{speciesName:'Homo sapiens'})-[:hasEvent]->(rle:ReactionLikeEvent{speciesName:'Homo sapiens'}),
            p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'}),
            (pe)-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name,
                last(labels(pe)) as Type, re.identifier as Id, head(re.geneName) as PrevId,
                re.databaseName AS Database, head([scores IN relationships(p) | type(scores)]) as Role
    ORDER BY Pathway, Reaction, Role, Type
    """,
    LEVELS[2]: """
    MATCH (pw:Pathway{speciesName:'Homo sapiens'})-[:hasEvent]->(rle:ReactionLikeEvent{speciesName:'Homo sapiens'}),
              p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'}),
              (pe)-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    WITH DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe, re, head([x IN relationships(p) | type(x)]) as Role
    OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
    WITH DISTINCT Pathway, Reaction, pe.stId as Entity, pe.displayName as Name, last(labels(pe)) as Type,                
    CASE 
    WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier 
    ELSE re.identifier
    END as Id, re.identifier as PrevId,
    mod.identifier as ptm_type, tm.coordinate as ptm_coordinate, re.databaseName as Database, Role
    ORDER BY ptm_type, ptm_coordinate
    WITH DISTINCT Pathway, Reaction, Entity, Name, Type, Id, PrevId,
    COLLECT(ptm_type + ":" + CASE WHEN ptm_coordinate IS NOT NULL THEN ptm_coordinate ELSE "null" END) AS ptms,
    Database, Role
    RETURN DISTINCT Pathway, Reaction, Entity, Name, Type, (Id+ptms) as Id, PrevId, Database, Role
    ORDER BY Pathway, Reaction, Role
    """,
    "sm": """
    MATCH (pw:Pathway{speciesName:'Homo sapiens'})-[:hasEvent]->(rle:ReactionLikeEvent{speciesName:'Homo sapiens'}),
              p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:SimpleEntity),
              (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)
        RETURN DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name,
                        last(labels(pe)) as Type, pe.displayName as Id, rd.displayName AS Database,  
                        head([scores IN relationships(p) | type(scores)]) as Role
        ORDER BY Pathway, Reaction, Role, Type
    """
}

QUERIES_COMPONENTS = {
    LEVELS[0]: """
    MATCH (c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT c.stId as Complex, pe.stId AS Entity, pe.displayName AS Name, last(labels(pe)) as Type, head(re.geneName) as Id
    ORDER BY Complex
    """,
    LEVELS[1]: """
    MATCH (c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT c.stId as Complex, pe.stId AS Entity, pe.displayName AS Name, last(labels(pe)) as Type, re.identifier as Id, head(re.geneName) as PrevId
    ORDER BY Complex
    """,
    LEVELS[2]: """
    MATCH (c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    WITH DISTINCT c, pe, last(labels(pe)) as Type, re
    OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
    WITH DISTINCT c.stId as Complex, 
                  pe.stId AS Entity, 
                  pe.displayName AS Name,
                  Type,
                  CASE 
                    WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier 
                    ELSE re.identifier
                  END as Id, re.identifier as PrevId,
                  mod.identifier as ptm_type,
                  tm.coordinate as ptm_coordinate
    ORDER BY ptm_type, ptm_coordinate
    WITH DISTINCT Complex, Entity, Name, Type, Id, PrevId,
                    COLLECT(
                        ptm_type + ":" + CASE WHEN ptm_coordinate IS NOT NULL THEN ptm_coordinate ELSE "null" END
                    ) AS ptms   
    RETURN DISTINCT Complex, Entity, Name, Type, (Id+ptms) as Id, PrevId
    ORDER BY Complex
    """,
    "sm": """
    MATCH (c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:SimpleEntity)
    RETURN DISTINCT c.stId as Complex, pe.stId AS Entity, pe.displayName AS Name, last(labels(pe)) as Type, pe.displayName as Id
    ORDER BY Complex
    """
}

QUERY_REACTIONS_ONLY_WITH_EWAS_PARTICIPANTS = """
MATCH p = (rle:ReactionLikeEvent{speciesName:"Homo sapiens"})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity)
WITH DISTINCT rle.stId as Reaction, collect(pe.stId) as Entity, collect(last(labels(pe))) as Type, collect( pe.displayName) as names
WHERE size(Type) = 1 AND "EntityWithAccessionedSequence" in Type
RETURN Reaction, Entity, Type, names
"""

QUERY_REACTIONS_WITH_ONLY_SMALL_MOLECULE_PARTICIPANTS = """
MATCH p = (rle:ReactionLikeEvent{speciesName:"Homo sapiens"})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity)
WITH DISTINCT rle.stId as Reaction,  collect(DISTINCT pe.stId) as Entity, collect(DISTINCT last(labels(pe))) as Type, collect(DISTINCT pe.displayName) as names
WHERE size(Type) <= 1 AND "SimpleEntity" in Type
RETURN Reaction, Entity, Type, names
"""


def get_pathway_name(pathway):
    query = f"MATCH (p:Pathway{{stId:\"{pathway}\", speciesName:\"Homo sapiens\"}}) RETURN p.displayName as Name"
    return get_query_result(query)


def get_pathways():
    query = """
    MATCH (p:Pathway{speciesName:"Homo sapiens"})
    RETURN p.stId as stId, p.displayName as displayName
    """
    return get_query_result(query)


def get_low_level_pathways():
    query = """
    // Gets all low level pathways for human
    MATCH (p:Pathway{speciesName:"Homo sapiens"})
    WHERE NOT (p)-[:hasEvent]->(:Pathway)
    RETURN p.stId as stId, p.displayName as displayName LIMIT 5
    """
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


def fix_neo4j_values(df, level):
    """
    Corrects format of some fields of the resulting records of the query to match the text structure

    For the 'Name' and 'Id' column: remove the subcellular location and spaces
    For the proteoforms: convert the list of isoform + ptms into a single string

    :param df: Dataframe with at least columns 'Id', 'Type, 'Name'
    :param level: "genes", "proteins" or "proteoforms"
    :return: Pandas dataframe with the values fixed
    """
    if len(df) == 0:
        return df

    df['Id'] = df.apply(lambda x: re.sub(r'\s*\[[\w\s]*\]\s*', '', x.Id) if x.Type == 'SimpleEntity' else x.Id, axis=1)
    if level == "proteoforms":
        df['Id'] = df['Id'].apply(make_proteoform_string)
    df['Name'] = df['Name'].apply(lambda x: re.sub("\s*\[[\s\w]*\]\s*", '', x))
    return df


def get_reaction_participants(level, v=False):
    """
    Get list of participant molecules in the reactions for all pathways in Reactome

    :param sm: Bool to show simple entities or not
    :param level: String "genes", "proteins" or "proteoforms"
    :param v: Bool Show extra console messages
    :return: pandas dataframe with one participant per record
    Columns: Pathway, Reaction, Entity, Name, Type, Id, Database, Role

    * Notice that records are sorted by reaction for easy traversal later.
    """
    if level not in LEVELS:
        raise Exception

    query = QUERIES_PARTICIPANTS[level]
    if (v):
        print(query)

    df = get_query_result(query)
    df = fix_neo4j_values(df, level)

    return df


def get_reaction_participants_by_pathway(pathway, level, sm, v=False):
    """
    Get list of participant molecules in the reactions of a pathway from the graph database

    :param pathway: String of the stId of the pathway
    :param sm: Bool to show simple entities or not
    :param level: String "genes", "proteins" or "proteoforms"
    :param v: Bool Show extra console messages
    :return: pandas dataframe with one participant per record
    Columns: Reaction, Entity, Name, Type, Id, Database, Role

    * Notice that records are sorted by reaction for easy traversal later.
    """
    if level not in LEVELS:
        raise Exception

    if (v):
        print(f"\n\nQuerying {level} participants of pathway {pathway}...\n\n")
    query = ""
    if level in ["genes", "proteins"]:
        query = f"""
        MATCH (pw:Pathway{{stId:"{pathway}"}})-[:hasEvent]->(rle:ReactionLikeEvent),
              p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity),
              (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)
        WHERE last(labels(pe)) IN ["EntityWithAccessionedSequence" """
        if sm:
            query += ", \"SimpleEntity\""
        query += """
        ]
        RETURN DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name, 
                        last(labels(pe)) as Type,
                        CASE 
                            WHEN last(labels(pe)) = \"SimpleEntity\" THEN pe.displayName """
        if level == "genes":
            query += "      WHEN last(labels(pe)) = \"EntityWithAccessionedSequence\" THEN head(re.geneName) ELSE re.identifier END as Id,"
        else:
            query += "      WHEN last(labels(pe)) = \"EntityWithAccessionedSequence\" THEN re.identifier ELSE re.identifier END as Id, head(re.geneName) as PrevId, "
        query += """ 
                        rd.displayName AS Database, 
                        head([scores IN relationships(p) | type(scores)]) as Role
        ORDER BY Pathway, Reaction, Role, Type
        """
    else:
        query += f"""
        MATCH (pw:Pathway{{stId:"{pathway}"}})-[:hasEvent]->(rle:ReactionLikeEvent),
              p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity),
              (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)
        WHERE last(labels(pe)) IN ["EntityWithAccessionedSequence" """
        if sm:
            query += ", \"SimpleEntity\""
        query += """
        ]
        WITH DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe, re, head([x IN relationships(p) | type(x)]) as Role
        OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
        WITH DISTINCT Pathway, Reaction, pe.stId as Entity, pe.displayName as Name, last(labels(pe)) as Type,
                      CASE 
                        WHEN last(labels(pe)) = "SimpleEntity" THEN pe.displayName  
                        WHEN last(labels(pe)) = "EntityWithAccessionedSequence" THEN 
                            CASE 
                                WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier 
                                ELSE re.identifier
                            END
                      END  as Id,
                      mod.identifier as ptm_type, tm.coordinate as ptm_coordinate, re.databaseName as Database, Role
        ORDER BY ptm_type, ptm_coordinate
        WITH DISTINCT Pathway, Reaction, Entity, Name, Type, Id, Id as PrevId,
                      COLLECT(
                          ptm_type + ":" + CASE WHEN ptm_coordinate IS NOT NULL THEN ptm_coordinate ELSE "null" END
                      ) AS ptms,
                      Database, Role
        RETURN DISTINCT Pathway, Reaction, Entity, Name, Type, CASE WHEN Type = "SimpleEntity" THEN Id ELSE (Id+ptms) END as Id, PrevId, Database, Role
		ORDER BY Pathway, Reaction, Role
        """

    if (v):
        print(query)

    df = get_query_result(query)
    df = fix_neo4j_values(df, level)

    return df


def get_complex_components(level, v=False):
    """
        Get list of complex components participating in all pathways of Reactome

        :param sm: Bool to show simple entities or not
        :param level: String "genes", "proteins" or "proteoforms"
        :param v: Bool Show extra console messages
        :return: pandas dataframe with one component per record
        Columns: Pathway, Reaction, Entity, Name, Type, Id, Database, Role

        * Notice that records are sorted by complex for easy traversal later.
        """
    if (v):
        print(QUERIES_COMPONENTS[level])

    df = get_query_result(QUERIES_COMPONENTS[level])
    df = fix_neo4j_values(df, level)

    return df


def get_complex_components_by_pathway(pathway, level, sm, v):
    """
    Get list of complex components participating in the pathway from the graph database

    :param pathway: String of the stId of the pathway
    :param sm: Bool to show simple entities or not
    :param level: String "genes", "proteins" or "proteoforms"
    :param v: Bool Show extra console messages
    :return: pandas dataframe with one component per record
    Columns: Reaction, Entity, Name, Type, Id, Database, Role

    * Notice that records are sorted by complex for easy traversal later.
    """

    if level in ["genes", "proteins"]:
        query = f"""
        // Get Pathway complex participants
        MATCH (p:Pathway{{stId:"{pathway}"}})-[:hasEvent]->(rle:ReactionLikeEvent)
        WITH rle
        MATCH (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(c:Complex)
        WITH c
        MATCH (c)-[:hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity)-[:referenceEntity]->(re:ReferenceEntity)
        WHERE last(labels(pe)) in ["EntityWithAccessionedSequence" """
        if sm:
            query += ", \"SimpleEntity\""
        query += """
        ]
        RETURN DISTINCT c.stId as Complex, pe.stId AS Entity, pe.displayName AS Name, last(labels(pe)) as Type, 
        CASE WHEN last(labels(pe)) = \"SimpleEntity\" THEN pe.displayName """
        if level == "genes":
            query += " WHEN last(labels(pe)) = \"EntityWithAccessionedSequence\" THEN head(re.geneName) ELSE re.identifier END as Id "
        else:
            query += " WHEN last(labels(pe)) = \"EntityWithAccessionedSequence\" THEN re.identifier ELSE re.identifier END as Id, head(re.geneName) as PrevId "
        query += """ 
        ORDER BY Complex
        """
    else:
        query = f"""
        MATCH (:Pathway{{stId:"{pathway}"}})-[:hasEvent]->(rle:Reaction{{speciesName:'Homo sapiens'}}),
        p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(c:Complex),
        (c)-[:hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity)-[:referenceEntity]->(re:ReferenceEntity)
        WHERE last(labels(pe)) in ["EntityWithAccessionedSequence" """
        if sm:
            query += ", \"SimpleEntity\""
        query += """
        ]
        WITH DISTINCT c, pe, last(labels(pe)) as Type, re
        OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
        WITH DISTINCT c.stId as Complex, 
                      pe.stId AS Entity, 
                      pe.displayName AS Name,
                      Type,
                      CASE 
                        WHEN Type = "SimpleEntity" THEN pe.displayName  
                        WHEN Type = "EntityWithAccessionedSequence" THEN 
                            CASE 
                                WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier 
                                ELSE re.identifier
                            END
                      END  as Id,
                      mod.identifier as ptm_type,
                      tm.coordinate as ptm_coordinate
        ORDER BY ptm_type, ptm_coordinate
        WITH DISTINCT Complex, Entity, Name, Type, Id, Id as PrevId,
                        COLLECT(
                            ptm_type + ":" + CASE WHEN ptm_coordinate IS NOT NULL THEN ptm_coordinate ELSE "null" END
                        ) AS ptms   
        RETURN DISTINCT Complex, Entity, Name, Type, CASE WHEN Type = "SimpleEntity" THEN Id ELSE (Id+ptms) END as Id, PrevId
        ORDER BY Complex
        """

    if (v):
        print(query)

    df = get_query_result(query)
    df = fix_neo4j_values(df, level)

    return df


def get_reactions():
    query = "MATCH (rle:ReactionLikeEvent{speciesName:\"Homo sapiens\"}) RETURN rle.stId as stId"
    return get_query_result(query)


def get_complexes():
    query = "MATCH (c:Complex{speciesName:\"Homo sapiens\"}) RETURN c.stId as stId"
    return get_query_result(query)


def get_complex_components_by_complex(complex, level, showSmallMolecules, verbose=True):
    if level in ["genes", "proteins"]:
        query = f"""
        MATCH (c:Complex{{stId:"{complex}"}})-[:hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity)-[:referenceEntity]->(re:ReferenceEntity)
        WHERE last(labels(pe)) in ["EntityWithAccessionedSequence" """
        if showSmallMolecules:
            query += ", \"SimpleEntity\""
        query += """
        ]
        RETURN DISTINCT c.stId as Complex, pe.stId AS Entity, pe.displayName AS Name, last(labels(pe)) as Type, 
        CASE WHEN last(labels(pe)) = \"SimpleEntity\" THEN pe.displayName """
        if level == "genes":
            query += " WHEN last(labels(pe)) = \"EntityWithAccessionedSequence\" THEN head(re.geneName) ELSE re.identifier END as Id "
        else:
            query += " WHEN last(labels(pe)) = \"EntityWithAccessionedSequence\" THEN re.identifier ELSE re.identifier END as Id, head(re.geneName) as PrevId "
        query += """ 
        ORDER BY Complex
        """
    else:
        query = f"""
        MATCH (c:Complex{{stId:"{complex}"}})-[:hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity)-[:referenceEntity]->(re:ReferenceEntity)
        WHERE last(labels(pe)) in ["EntityWithAccessionedSequence" """
        if showSmallMolecules:
            query += ", \"SimpleEntity\""
        query += """
        ]
        WITH DISTINCT c, pe, last(labels(pe)) as Type, re
        OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
        WITH DISTINCT c.stId as Complex, 
                      pe.stId AS Entity, 
                      pe.displayName AS Name,
                      Type,
                      CASE 
                        WHEN Type = "SimpleEntity" THEN pe.displayName  
                        WHEN Type = "EntityWithAccessionedSequence" THEN 
                            CASE 
                                WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier 
                                ELSE re.identifier
                            END
                      END  as Id,
                      mod.identifier as ptm_type,
                      tm.coordinate as ptm_coordinate
        ORDER BY ptm_type, ptm_coordinate
        WITH DISTINCT Complex, Entity, Name, Type, Id, Id as PrevId,
                        COLLECT(
                            ptm_type + ":" + CASE WHEN ptm_coordinate IS NOT NULL THEN ptm_coordinate ELSE "null" END
                        ) AS ptms   
        RETURN DISTINCT Complex, Entity, Name, Type, CASE WHEN Type = "SimpleEntity" THEN Id ELSE (Id+ptms) END as Id, PrevId
        ORDER BY Complex
        """

    if (verbose):
        print(query)

    df = get_query_result(query)
    df = fix_neo4j_values(df, level)

    return df


def get_reaction_participants_by_reaction(reaction, level, showSmallMolecules, verbose=False):
    if level not in LEVELS:
        raise Exception

    if (verbose):
        print(f"\n\nQuerying {level} participants of reaction {reaction}...\n\n")
    query = ""
    if level in ["genes", "proteins"]:
        query = f"""
            MATCH p = (rle:ReactionLikeEvent{{stId:"{reaction}"}})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity),
                  (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)
            WHERE last(labels(pe)) IN ["EntityWithAccessionedSequence" """
        if showSmallMolecules:
            query += ", \"SimpleEntity\""
        query += """
            ]
            RETURN DISTINCT rle.stId as Reaction, 
                            pe.stId as Entity, 
                            pe.displayName as Name, 
                            last(labels(pe)) as Type,
                            CASE 
                                WHEN last(labels(pe)) = \"SimpleEntity\" THEN pe.displayName """
        if level == "genes":
            query += "      WHEN last(labels(pe)) = \"EntityWithAccessionedSequence\" THEN head(re.geneName) ELSE re.identifier END as Id, "
        else:
            query += "      WHEN last(labels(pe)) = \"EntityWithAccessionedSequence\" THEN re.identifier ELSE re.identifier END as Id, head(re.geneName) as PrevId, "
        query += """ 
                            rd.displayName AS Database, 
                            head([scores IN relationships(p) | type(scores)]) as Role
            ORDER BY Reaction, Role, Type
            """
    else:
        query += f"""
            MATCH p = (rle:ReactionLikeEvent{{stId:"{reaction}"}})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity),
                  (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)
            WHERE last(labels(pe)) IN ["EntityWithAccessionedSequence" """
        if showSmallMolecules:
            query += ", \"SimpleEntity\""
        query += """
            ]
            WITH DISTINCT rle.stId as Reaction,
    			  pe, re, 
                  head([x IN relationships(p) | type(x)]) as Role
            ORDER BY Reaction, Role
            OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
            WITH DISTINCT Reaction, 
                          pe.stId as Entity, 
                          pe.displayName as Name,
                          last(labels(pe)) as Type,
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
                          re.databaseName as Database,
                          Role
            ORDER BY ptm_type, ptm_coordinate
            WITH DISTINCT Reaction,
                          Entity,
                          Name,
                          Type,
                          Id, Id as PrevId, 
                          COLLECT(
                              ptm_type + ":" + CASE WHEN ptm_coordinate IS NOT NULL THEN ptm_coordinate ELSE "null" END
                          ) AS ptms,
                          Database,
                          Role
            ORDER BY Reaction, Role, Id
            RETURN DISTINCT Reaction, Entity, Name, Type, CASE WHEN Type = "SimpleEntity" THEN Id ELSE (Id+ptms) END as Id, PrevId, Database, Role
    		ORDER BY Reaction, Role
            """

    if (verbose):
        print(query)

    df = get_query_result(query)
    df = fix_neo4j_values(df, level)

    return df
