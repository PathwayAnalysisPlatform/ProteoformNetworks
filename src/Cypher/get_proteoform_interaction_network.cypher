// Not sure if this query is correct
MATCH (re1:ReferenceEntity)<-[:referenceEntity]-(pe1:PhysicalEntity)<-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]-(r:ReactionLikeEvent)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe2:PhysicalEntity)-[:referenceEntity]->(re2:ReferenceEntity)
WHERE r.stId = "R-HSA-8863895" 
AND re1.databaseName = re2.databaseName = "UniProt" 
AND r.speciesName = pe1.speciesName = pe2.speciesName = "Homo sapiens" 
AND re1.identifier <> re2.identifier
WITH re1, pe1, r, pe2, re2
OPTIONAL MATCH (pe1)-[:hasModifiedResidue]->(tm1:TranslationalModification)-[:psiMod]->(mod1:PsiMod),
               (pe2)-[:hasModifiedResidue]->(tm2:TranslationalModification)-[:psiMod]->(mod2:PsiMod)
WITH DISTINCT r.stId as reaction, 
            pe1.stId as pe1, re1.identifier as protein1, re1.variantIdentifier as isoform1, tm1.coordinate as coordinate1, mod1.identifier as type1,
            pe2.stId as pe2, re2.identifier as protein2, re2.variantIdentifier as isoform2, tm2.coordinate as coordinate2, mod2.identifier as type2
            ORDER BY type1, type2, coordinate1, coordinate2
WITH DISTINCT reaction, 
               pe1, 
               protein1,
               CASE WHEN isoform1 IS NOT NULL THEN isoform1 ELSE protein1 END as isoform1,
               COLLECT(type1 + ":" + CASE WHEN coordinate1 IS NOT NULL THEN coordinate1 ELSE "null" END) AS ptms1,
               pe2,
               protein2,
               CASE WHEN isoform2 IS NOT NULL THEN isoform2 ELSE protein2 END as isoform2,
               COLLECT(type2 + ":" + CASE WHEN coordinate2 IS NOT NULL THEN coordinate2 ELSE "null" END) AS ptms2
RETURN DISTINCT pe1, (isoform1 + ptms1) as proteoform1, reaction, pe2, (isoform2 + ptms2) as proteoform2
                ORDER BY reaction, proteoform1, proteoform2
                