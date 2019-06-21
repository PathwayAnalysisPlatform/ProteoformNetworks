MATCH (peie:InstanceEdit)-[:created]-(pe:PhysicalEntity{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})
WITH pe, re, peie
OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod), 
               (tm)-[:created]-(tmie:InstanceEdit)
WITH DISTINCT pe.stId AS pe,
			  peie.dateTime as pe_date,
              re.identifier AS protein,
              re.variantIdentifier AS isoform,
              tm.coordinate as coordinate, 
              mod.identifier as type,
              CASE WHEN tmie.dateTime < peie.dateTime THEN peie.dateTime ELSE tmie.dateTime END as tm_date
ORDER BY type, coordinate
WITH DISTINCT pe,
              pe_date,
		        protein,
              CASE WHEN isoform IS NOT NULL THEN isoform ELSE protein END as isoform,
              collect(tm_date) as tm_dates,
              collect(type + ":" + CASE WHEN coordinate IS NOT NULL THEN coordinate ELSE "null" END) AS ptms
WITH DISTINCT pe, (isoform + ptms) as proteoform, 
                CASE WHEN size(tm_dates) = 0 THEN pe_date ELSE head(min(tm_dates)) END as date  
WITH DISTINCT pe, proteoform, collect(date) as date
RETURN DISTINCT pe, proteoform, head(min(date)) as date
ORDER BY date DESC
