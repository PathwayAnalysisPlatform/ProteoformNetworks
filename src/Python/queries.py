from config import LEVELS, proteoforms, genes, proteins, sm

QUERIES_PARTICIPANTS = {
    genes: """
    MATCH (pw:Pathway{speciesName:'Homo sapiens'})-[:hasEvent]->(rle:ReactionLikeEvent{speciesName:'Homo sapiens'}),
          p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'}),
          (pe)-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name, 
                    last(labels(pe)) as Type, head(re.geneName) as Id, re.databaseName AS Database, 
                    head([scores IN relationships(p) | type(scores)]) as Role
    ORDER BY Pathway, Reaction, Role, Type
    """,
    proteins: """
    MATCH (pw:Pathway{speciesName:'Homo sapiens'})-[:hasEvent]->(rle:ReactionLikeEvent{speciesName:'Homo sapiens'}),
            p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'}),
            (pe)-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name,
                last(labels(pe)) as Type, re.identifier as Id, head(re.geneName) as PrevId,
                re.databaseName AS Database, head([scores IN relationships(p) | type(scores)]) as Role
    ORDER BY Pathway, Reaction, Role, Type
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
    ORDER BY Pathway, Reaction, Role
    """,
    sm: """
    MATCH (pw:Pathway{speciesName:'Homo sapiens'})-[:hasEvent]->(rle:ReactionLikeEvent{speciesName:'Homo sapiens'}),
              p = (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:SimpleEntity),
              (pe)-[:referenceEntity]->(re:ReferenceEntity)-[:referenceDatabase]->(rd:ReferenceDatabase)
        RETURN DISTINCT pw.stId as Pathway, rle.stId as Reaction, pe.stId as Entity, pe.displayName as Name, last(labels(pe)) as Type, "sm_" + pe.displayName as Id, rd.displayName AS Database,  
                        head([scores IN relationships(p) | type(scores)]) as Role
        ORDER BY Pathway, Reaction, Role, Type
    """
}

def get_query_participants_by_pathway(level, pathway="", reaction=""):
    query = QUERIES_PARTICIPANTS[level]
    if len(pathway) > 0:
        query = query.replace("Pathway{speciesName:'Homo sapiens'}", f"Pathway{{speciesName:'Homo sapiens', stId:'{pathway}'}}")
    if len(reaction) > 0:
        query = query.replace("ReactionLikeEvent{speciesName:'Homo sapiens'}", f"ReactionLikeEvent{{speciesName:'Homo sapiens', stId:'{reaction}'}}")
    return query

def get_query_components(level, pathway="", reaction=""):
    query = QUERIES_COMPONENTS[level]
    if len(pathway) > 0:
        query = query.replace("Pathway{speciesName:'Homo sapiens'}", f"Pathway{{speciesName:'Homo sapiens', stId:'{pathway}'}}")
    if len(reaction) > 0:
        query = query.replace("ReactionLikeEvent{speciesName:'Homo sapiens'}", f"ReactionLikeEvent{{speciesName:'Homo sapiens, stId:'{reaction}'}}")
    return query

QUERIES_COMPONENTS = {
    genes: """
    MATCH (c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT c.stId as Complex, pe.stId AS Entity, pe.displayName AS Name, last(labels(pe)) as Type, head(re.geneName) as Id
    ORDER BY Complex
    
    MATCH (p:Pathway{speciesName:'Homo sapiens'})-[:hasEvent*]->(r:ReactionLikeEvent{speciesName:'Homo sapiens'})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator*]->(c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT r.displayName as Reaction, r.stId AS ReactionStId,  pe.stId as Participant, c.stId as Complex, labels(c) as ParticipantType, pe.stId AS Entity, pe.displayName AS Name, last(labels(pe)) as Type, head(re.geneName) as Id
    ORDER BY Reaction, Complex, Entity
    """,
    proteins: """
    MATCH (c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
    RETURN DISTINCT c.stId as Complex, pe.stId AS Entity, pe.displayName AS Name, last(labels(pe)) as Type, re.identifier as Id, head(re.geneName) as PrevId
    ORDER BY Complex
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
    ORDER BY Complex
    """,
    sm: """
    MATCH (c:Complex{speciesName:'Homo sapiens'})-[:hasComponent|hasMember|hasCandidate*]->(pe:SimpleEntity)
    RETURN DISTINCT c.stId as Complex, pe.stId AS Entity, pe.displayName as Name, last(labels(pe)) as Type, "sm_" + pe.displayName as Id, "sm_" + c.stId + "_" + pe.displayName as UniqueId
    ORDER BY Complex
    """
}

QUERY_GET_COMPLEXES_BY_PATHWAY_OR_REACTION = """
MATCH (p:Pathway{stId:"R-HSA-983169"})-[:hasEvent*]->(r:ReactionLikeEvent{stId:"R-HSA-8863895"})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator*]->(pe:Complex)
RETURN DISTINCT r.stId AS Reaction, pe.stId as Participant, pe.displayName AS DisplayName, labels(pe)
"""

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

QUERY_GET_ALL_PROTEOFORMS = """
MATCH (pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
WITH DISTINCT pe, re
OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
WITH DISTINCT pe.stId as Entity, 
			  pe.displayName as Name,              
              CASE WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier ELSE re.identifier END as Id, 
    		  mod.identifier as ptm_type, 
              tm.coordinate as ptm_coordinate
ORDER BY ptm_type, ptm_coordinate
WITH DISTINCT Entity, Name, Id, COLLECT(ptm_type + ":" + CASE WHEN ptm_coordinate IS NOT NULL THEN ptm_coordinate ELSE "null" END) AS ptms
WITH DISTINCT Entity, Name, (Id+ptms) as Id
RETURN DISTINCT Id ORDER BY Id
"""

QUERY_GET_PROTEOFORMS_OF_EACH_PROTEIN = """
MATCH (pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
WITH DISTINCT pe, re
OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
WITH DISTINCT pe.stId as Entity, 
			  pe.displayName as Name, 
              re.identifier as Protein,
              CASE WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier ELSE re.identifier END as Id, 
    		  mod.identifier as ptm_type, 
              tm.coordinate as ptm_coordinate
ORDER BY Protein, ptm_type, ptm_coordinate
WITH DISTINCT Entity, Name, Protein, Id, COLLECT(ptm_type + ":" + CASE WHEN ptm_coordinate IS NOT NULL THEN ptm_coordinate ELSE "null" END) AS ptms
WITH DISTINCT Entity, Name, Protein, (Id+ptms) as Proteoform ORDER BY Proteoform
WITH DISTINCT Protein, COLLECT(DISTINCT Proteoform) as Proteoforms ORDER By Protein
RETURN DISTINCT Protein, Proteoforms
"""

QUERY_GET_NUM_PROTEOFORMS_PER_PROTEIN = """
MATCH (pe:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
WITH DISTINCT pe, re
OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
WITH DISTINCT pe.stId as Entity, 
			  pe.displayName as Name, 
              re.identifier as Protein,
              CASE WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier ELSE re.identifier END as Id, 
    		  mod.identifier as ptm_type, 
              tm.coordinate as ptm_coordinate
ORDER BY Protein, ptm_type, ptm_coordinate
WITH DISTINCT Entity, Name, Protein, Id, COLLECT(ptm_type + ":" + CASE WHEN ptm_coordinate IS NOT NULL THEN ptm_coordinate ELSE "null" END) AS ptms
WITH DISTINCT Entity, Name, Protein, (Id+ptms) as Proteoform ORDER BY Proteoform
WITH DISTINCT Protein, COLLECT(DISTINCT Proteoform) as Proteoforms ORDER By Protein
WITH Protein, Proteoforms, size(Proteoforms) as NumProteoforms
WHERE NumProteoforms > 1
RETURN DISTINCT Protein, Proteoforms, NumProteoforms ORDER BY NumProteoforms DESC
"""

QUERY_GET_PATHWAYS_BY_PROTEIN = """
MATCH (p:Pathway{speciesName:"Homo sapiens"})-[:hasEvent*]->(rle:ReactionLikeEvent{speciesName:"Homo sapiens"}),
      (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity),
      (pe)-[:referenceEntity]->(re:ReferenceEntity{identifier:"P04049", databaseName:"UniProt"})
RETURN DISTINCT p.stId AS PathwayId, p.displayName AS Pathway, re.identifier AS Identifier
ORDER BY PathwayId, Identifier
"""


