MATCH (re:ReferenceEntity{databaseName:"UniProt"})<-[:referenceEntity]-(pe:PhysicalEntity{speciesName:"Homo sapiens"})
WITH re, pe
OPTIONAL MATCH (p:Pathway{speciesName:'Homo sapiens'})-[:hasEvent*]->(r:Reaction{speciesName:'Homo sapiens'})-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe)
RETURN re.identifier as protein, count(DISTINCT r) as reactionCount, count(DISTINCT p) as pathwayCount
