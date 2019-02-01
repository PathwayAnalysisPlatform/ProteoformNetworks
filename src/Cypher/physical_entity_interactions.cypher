MATCH (re1:ReferenceEntity)<-[:referenceEntity]-(pe1:PhysicalEntity)<-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]-(r:ReactionLikeEvent)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe2:PhysicalEntity)-[:referenceEntity]->(re2:ReferenceEntity)
WHERE r.stId = "R-HSA-8863895" AND re1.databaseName = re2.databaseName = "UniProt" AND r.speciesName = pe1.speciesName = pe2.speciesName = "Homo sapiens" AND re1.identifier <> re2.identifier
RETURN DISTINCT pe1.stId, pe1.displayName, r.stId, pe2.stId, pe2.displayName
ORDER BY r.stId, pe1.stId, pe2.stId
