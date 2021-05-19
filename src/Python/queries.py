from config import LEVELS, proteoforms, genes, proteins, sm

QUERIES_PARTICIPANTS = {
    genes: """
    MATCH (pw:Pathway{speciesName:'Homo sapiens'})-[:hasEvent]->(rle:ReactionLikeEvent{speciesName:'Homo sapiens'}),
          p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'}),
          (pe)-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name, 
                    last(labels(pe)) as Type, head(re.geneName) as Id, re.databaseName AS Database, 
                    head([scores IN relationships(p) | type(scores)]) as Role
    ORDER BY Pathway, Reaction, Role, Type LIMIT 100
    """,
    proteins: """
    MATCH (pw:Pathway{speciesName:'Homo sapiens'})-[:hasEvent]->(rle:ReactionLikeEvent{speciesName:'Homo sapiens'}),
            p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'}),
            (pe)-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name,
                last(labels(pe)) as Type, re.identifier as Id, head(re.geneName) as PrevId,
                re.databaseName AS Database, head([scores IN relationships(p) | type(scores)]) as Role
    ORDER BY Pathway, Reaction, Role, Type LIMIT 100
    """,
    proteoforms: """
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
    ORDER BY Pathway, Reaction, Role LIMIT 100
    """,
    sm: """
    MATCH (pw:Pathway{speciesName:'Homo sapiens'})-[:hasEvent]->(rle:ReactionLikeEvent{speciesName:'Homo sapiens'}),
              p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:SimpleEntity),
              (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)
        RETURN DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name,
                        last(labels(pe)) as Type, pe.displayName as Id, rd.displayName AS Database,  
                        head([scores IN relationships(p) | type(scores)]) as Role
        ORDER BY Pathway, Reaction, Role, Type LIMIT 100
    """
}

QUERIES_COMPONENTS = {
    genes: """
    MATCH (c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT c.stId as Complex, pe.stId AS Entity, pe.displayName AS Name, last(labels(pe)) as Type, head(re.geneName) as Id
    ORDER BY Complex LIMIT 100
    """,
    proteins: """
    MATCH (c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT c.stId as Complex, pe.stId AS Entity, pe.displayName AS Name, last(labels(pe)) as Type, re.identifier as Id, head(re.geneName) as PrevId
    ORDER BY Complex LIMIT 100
    """,
    proteoforms: """
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
    ORDER BY Complex LIMIT 100
    """,
    sm: """
    MATCH (c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:SimpleEntity)
    RETURN DISTINCT c.stId as Complex, pe.stId AS Entity, pe.displayName AS Name, last(labels(pe)) as Type, pe.displayName as Id
    ORDER BY Complex LIMIT 100
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
