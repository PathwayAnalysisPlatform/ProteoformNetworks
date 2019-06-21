MATCH (re1:ReferenceEntity)<-[:referenceEntity]-(pe1:PhysicalEntity)<-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]-(r:ReactionLikeEvent)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe2:PhysicalEntity)-[:referenceEntity]->(re2:ReferenceEntity),
(peie1:InstanceEdit)-[:created]-(pe1), 
(peie2:InstanceEdit)-[:created]-(pe2)
WHERE re1.databaseName = re2.databaseName = "UniProt" AND r.speciesName = pe1.speciesName = pe2.speciesName = "Homo sapiens" AND re1.identifier <> re2.identifier
RETURN DISTINCT pe1.stId as pe1, peie1.dateTime as pe1_date, pe2.stId as pe2, peie2.dateTime as pe2_date
