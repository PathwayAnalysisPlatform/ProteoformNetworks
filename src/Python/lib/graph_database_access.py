import re
from os import path

import pandas as pd
from neo4j import GraphDatabase

from config import proteoforms, LEVELS, proteins, genes
from queries import QUERIES_PARTICIPANTS, QUERIES_COMPONENTS, get_query_participants_by_pathway, \
    QUERY_GET_COMPLEXES_BY_PATHWAY_OR_REACTION


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
    df['Id'] = df.apply(lambda x: str(x.Id).replace(" ", "_").strip() if x.Type == 'SimpleEntity' else x.Id, axis=1)

    if "UniqueId" in df.columns:
        df['UniqueId'] = df.apply(
            lambda x: re.sub(r'\s*\[[\w\s]*\]\s*', '', x.UniqueId) if x.Type == 'SimpleEntity' else x.Id, axis=1)
        df['UniqueId'] = df.apply(lambda x: str(x.UniqueId).replace(" ", "_") if x.Type == 'SimpleEntity' else x.Id,
                                  axis=1)

    if level == proteoforms:
        df['Id'] = df.apply(
            lambda x: make_proteoform_string(x.Id) if x.Type == 'EntityWithAccessionedSequence' else x.Id, axis=1)
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


def get_empty_participants_dataframe(level):
    if level == genes:
        return pd.DataFrame(
            columns=["Pathway", "Reaction", "Entity", "Name", "Type", "Id", "Database", "Role"])
    elif level == proteins:
        return pd.DataFrame(
            columns=["Pathway", "Reaction", "Entity", "Name", "Type", "Id", "PrevId", "Database", "Role"])
    elif level == proteoforms:
        return pd.DataFrame(
            columns=["Pathway", "Reaction", "Entity", "Name", "Type", "Id", "PrevId", "Database", "Role"])
    else:
        return pd.DataFrame(
            columns=["Pathway", "Reaction", "Entity", "Name", "Type", "Id", "UniqueId", "Database", "Role"])


def get_participants_by_pathway(pathway, level, out_path=""):
    # print(f"Getting participants for pathway {pathway} for level {level}")

    filename = out_path + "participants/pathway_" + pathway + "_" + level + ".csv"

    participants = pd.DataFrame()
    if not path.exists(filename):
        query = get_query_participants_by_pathway(level, pathway)
        participants = get_query_result(query)
        participants = fix_neo4j_values(participants, level)
        if len(participants) == 0:
            participants = get_empty_participants_dataframe(level)
        participants.to_csv(filename)
    else:
        participants = pd.read_csv(filename)
        if len(participants) == 0:
            participants = get_empty_participants_dataframe(level)
    return participants


def get_empty_components_dataframe(level):
    if level == genes:
        return pd.DataFrame(columns=['Complex', 'Entity', 'Name', 'Type', 'Id'])
    elif level == proteins:
        return pd.DataFrame(columns=['Complex', 'Entity', 'Name', 'Type', 'Id', 'PrevId'])
    elif level == proteoforms:
        return pd.DataFrame(columns=['Complex', 'Entity', 'Name', 'Type', 'Id', 'PrevId'])
    else:
        return pd.DataFrame(columns=['Complex', 'Entity', 'Name', 'Type', 'Id', 'UniqueId'])


def get_components(level, location=""):
    filename = location + "records_complex_components_" + level + ".csv"

    if not path.exists(filename):
        components = get_query_result(QUERIES_COMPONENTS[level])
        components = fix_neo4j_values(components, level)
        components.to_csv(filename)
        if len(components) == 0:
            return get_empty_components_dataframe(level)
        return components
    else:
       components = pd.read_csv(filename)
       if len(components) == 0:
           components = get_empty_components_dataframe(level)
       return components


def get_complexes_by_pathway(pathway):
    query = QUERY_GET_COMPLEXES_BY_PATHWAY_OR_REACTION.replace("Pathway{speciesName:'Homo sapiens'}",
                                                               f"Pathway{{speciesName:'Homo sapiens', stId:'{pathway}'}}")
    complexes = get_query_result(query)
    return complexes


def get_components_by_pathway(pathway, level, out_path=""):
    # print(f"Getting components for pathway {pathway} for level {level}")
    if len(pathway) > 0:
        # Get list of participating complexes
        df_complexes = get_complexes_by_pathway(pathway)

        if len(df_complexes) > 0:
            # Get the components of each complex
            dfs_components = [get_complex_components_by_complex(complex, level, out_path) for complex in
                              df_complexes["Complex"]]
            components = pd.concat(dfs_components)
            return components
        else:
            return get_empty_components_dataframe(level)
    else:
        query = QUERIES_COMPONENTS[level]
        components = get_query_result(query)
        components = fix_neo4j_values(components, level)
        if len(components) == 0:
            components = get_empty_components_dataframe(level)
        return components


if __name__ == "__main__":
    pathway = "R-HSA-70171"
    query = f"MATCH (p:Pathway{{stId:\"{pathway}\"}})-[:hasEvent]->(rle:Reaction{{speciesName:'Homo sapiens'}}) RETURN rle.stId"
    df = get_query_result(query)
    print(df)


def get_pathway_name(pathway):
    query = f"MATCH (p:Pathway{{stId:\"{pathway}\", speciesName:\"Homo sapiens\"}})-[:hasEvent]->(rle:ReactionLikeEvent{{speciesName:\"Homo sapiens\"}})" \
            f" RETURN DISTINCT p.displayName as Name"
    return get_query_result(query)


def get_pathways():
    query = """
    MATCH (p:Pathway{speciesName:"Homo sapiens"})-[:hasEvent]->(rle:ReactionLikeEvent{speciesName:"Homo sapiens"})
    RETURN DISTINCT p.stId as stId, p.displayName as displayName
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


def get_reactions():
    query = "MATCH (rle:ReactionLikeEvent{speciesName:\"Homo sapiens\"}) RETURN rle.stId as stId"
    return get_query_result(query)


def get_complexes():
    query = "MATCH (c:Complex{speciesName:\"Homo sapiens\"}) RETURN c.stId as stId"
    return get_query_result(query)


def get_complex_components_by_complex(complex, level, out_path=""):
    # print(f"\tGetting components of complex: {complex}")
    filename = out_path + "complexes/complex_" + complex + "_" + level + ".csv"

    if not path.exists(filename):
        query = QUERIES_COMPONENTS[level].replace(
            "Complex{speciesName:'Homo sapiens'}",
            f"Complex{{speciesName:'Homo sapiens', stId:'{complex}'}}")
        components = get_query_result(query)
        components = fix_neo4j_values(components, level)
        if len(components) == 0:
            components = get_empty_components_dataframe(level)
        components.to_csv(filename)
        return components
    else:
        components = pd.read_csv(filename)
        if len(components) == 0:
            components = get_empty_components_dataframe(level)
        components.to_csv(filename)
        return components


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
