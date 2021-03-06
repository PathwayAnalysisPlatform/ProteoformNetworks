import re
import os.path
from os import path

import pandas as pd
from neo4j import GraphDatabase

from config import proteoforms, LEVELS
from queries import QUERIES_PARTICIPANTS, QUERIES_COMPONENTS, get_query_participants_by_pathway


def get_query_result(query):
    """
    Get pandas dataframe with the result of the query to Reactome Graph Database
    :param query: Cypher string with the query to get the data
    :return: pandas dataframe with the result records
    """
    db = GraphDatabaseAccess("bolt://localhost", "neo4j", "")
    df = db.get_result(query)
    db.close()
    return df;


class GraphDatabaseAccess:

    def __init__(self, uri, user, password):
        self.driver = GraphDatabase.driver(uri, auth=(user, password), encrypted=False)

    def close(self):
        self.driver.close()

    def get_result(self, query):
        with self.driver.session() as session:
            records = session.read_transaction(self.get_records, query)
            if records:
                return pd.DataFrame([r.values() for r in records], columns=records[0].keys())
        return pd.DataFrame()

    @staticmethod
    def get_records(tx, query):
        result = []
        for record in tx.run(query):
            result.append(record)
        return result


def make_proteoform_string(value):
    """
    Create proteoform string in the simple format: isoform;ptm1,ptm2...
    Adds a ';' at the end of the proteoform when there are no ptms, to make sure the string represents a proteoform.
    Examples: 	["P36507", "00046:null", "00047:null"] or 	["P28482"]

    :param value: array of strings
    :return:
    """
    if type(value) == str:
        return value + ";"
    if len(value) == 1:
        return value[0] + ";"
    if len(value) == 2:
        return ";".join(value)
    else:
        isoform = value[0] + ";"
        ptms = ",".join(value[1:])
        return isoform + ptms
    print(f"Got a weird value: {value}")
    return value[0]


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
    df['Id'] = df.apply(lambda x: str(x.Id).replace(" ", "_") if x.Type == 'SimpleEntity' else x.Id, axis=1)

    # df['UniqueId'] = df.apply(lambda x: re.sub(r'\s*\[[\w\s]*\]\s*', '', x.UniqueId) if x.Type == 'SimpleEntity' else x.Id, axis=1)
    # df['UniqueId'] = df.apply(lambda x: str(x.UniqueId).replace(" ", "_") if x.Type == 'SimpleEntity' else x.Id, axis=1)

    if level == proteoforms:
        df['Id'] = df['Id'].apply(make_proteoform_string)
    df['Name'] = df['Name'].apply(lambda x: re.sub("\s*\[[\s\w]*\]\s*", '', x))
    return df


def get_participants(level, location=""):
    """
    Gets reaction participants from Reactome in a table with the columns: [Pathway, Reaction, Entity, Name, Type, Id, Database, Role]

    :param location: directory where to search for the csv files
    :param level: genes, proteins or proteoforms
    :return: pandas dataframe
    """

    filename = location + "records_reaction_participants_" + level + ".csv"

    if not path.exists(filename):
        participants = get_query_result(QUERIES_PARTICIPANTS[level])
        participants = fix_neo4j_values(participants, level)
        participants.to_csv(filename)
        return participants
    else:
        return pd.read_csv(filename)

def get_participants_by_pathway(level, pathway):
    query = get_query_participants_by_pathway(level, pathway)
    participants = get_query_result(query)
    participants = fix_neo4j_values(participants, level)
    return participants


def get_components(level, location=""):

    filename = location + "records_complex_components_" + level + ".csv"

    if not path.exists(filename):
        components = get_query_result(QUERIES_COMPONENTS[level])
        components = fix_neo4j_values(components, level)
        components.to_csv(filename)
        return components
    else:
        return pd.read_csv(filename)

def get_components_by_pathway(level, pathway):
        query = get_query_participants_by_pathway(level, pathway)
        participants = get_query_result(query)
        participants = fix_neo4j_values(participants, level)
        return participants


if __name__ == "__main__":
    pathway = "R-HSA-70171"
    query = f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:Reaction{{speciesName:'Homo sapiens'}}) RETURN rle.stId"
    df = get_query_result(query)
    print(df)


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
    return f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:Reaction{{speciesName:'Homo sapiens'}}) RETURN rle.stId as reaction"

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
