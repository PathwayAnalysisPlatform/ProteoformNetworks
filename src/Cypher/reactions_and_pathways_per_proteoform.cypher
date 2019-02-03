MATCH (pe:PhysicalEntity{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})
WITH DISTINCT pe, re
OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
WITH DISTINCT pe, re, tm.coordinate as coordinate, mod.identifier as type 
ORDER BY type, coordinate
WITH DISTINCT pe, re.identifier as protein, CASE WHEN re.variantIdentifier IS NOT NULL THEN re.variantIdentifier ELSE re.identifier END as isoform, COLLECT(CASE WHEN coordinate IS NOT NULL THEN coordinate ELSE "null" END + ":" + type) AS ptms
OPTIONAL MATCH (p:Pathway)-[:hasEvent*]->(r:Reaction)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe)
WHERE p.speciesName = 'Homo sapiens' AND r.speciesName = 'Homo sapiens' AND pe.speciesName = 'Homo sapiens'
RETURN DISTINCT protein, isoform, reduce(ptmsString = ";", ptm IN ptms| ptmsString + ptm + "," ) as ptms, size(ptms) as num_ptms, count(DISTINCT r) as reactionCount, count(DISTINCT p) as pathwayCount
