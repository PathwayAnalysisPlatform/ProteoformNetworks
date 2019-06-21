# Steps to execute

## Examples of pathways with gene level only overlap

* Get the protein and proteoform lists from the graph database in NEO4J console.

Gene list:
~~~~
MATCH (ewas:EntityWithAccessionedSequence{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})
WITH re.identifier as protein, re.geneName as genes
WHERE size(genes) > 0  
UNWIND genes as gene
RETURN DISTINCT gene
~~~~

Protein list:
~~~~
MATCH (pe:PhysicalEntity{speciesName:"Homo sapiens"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:"UniProt"})
RETURN DISTINCT re.identifier as protein
~~~~

Proteoform list:
~~~~
MATCH (pe:PhysicalEntity{speciesName:'Homo sapiens'})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})
WITH DISTINCT pe, re
OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
WITH DISTINCT pe.stId AS physicalEntity,
                re.identifier AS protein,
                re.variantIdentifier AS isoform,
                tm.coordinate as coordinate, 
                mod.identifier as type ORDER BY type, coordinate
WITH DISTINCT physicalEntity,
				protein,
                CASE WHEN isoform IS NOT NULL THEN isoform ELSE protein END as isoform,
                COLLECT(type + ":" + CASE WHEN coordinate IS NOT NULL THEN coordinate ELSE "null" END) AS ptms
                RETURN DISTINCT isoform, ptms
                ORDER BY isoform, ptms
~~~~

Then convert the proteoform format from NEO4J to SIMPLE. Use PathwayMatcher class called ProteoformFormatConverter.
~~~~
java -cp PathwayMatcher.jar no.uib.pap.pathwaymatcher.tools.ProteoformFormatConverter reactome/ all_proteoforms_neo4j.csv all_proteoforms.csv
~~~~

* Find out the gene, protein and proteoform members of each pathway. For this we execute PathwayMatcher and get the whole search result.

Genes:
~~~~
java -jar PathwayMatcher.jar -t gene -i reactome/all_genes.csv -o reactome/all_genes/
~~~~

Proteins:
~~~~
java -jar PathwayMatcher.jar -t uniprot -i reactome/all_proteins.csv -o reactome/all_proteins/
~~~~

Proteoforms:
~~~~
java -jar PathwayMatcher.jar -t proteoform -i reactome/all_proteoforms.csv -o reactome/all_proteoforms/ -m strict
~~~~

* Execute main C++ program to create the pathway sets and calculate overlaps:
~~~~
g++ src/3_rule_out_gene_centric_overlap/rule_out_gene_centric_overlap.cpp src/main.cpp -O3 -o Debug/analysis.exe
~~~~

From the report file ("reports/3_rule_out_gene_centric_overlap_analysis.txt") choose a pair of pathways to see the variation in overlap.

* Get members of both pathways at the different levels of granularity. The pathway names, isoforms and post translational modifications responsible for decomposing the gene level only overlap will show up.

Genes:
~~~~
MATCH (pathway:Pathway{speciesName:"Homo sapiens"})-[:hasEvent*]->(rle:Reaction{speciesName:"Homo sapiens"}),
      (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:"Homo sapiens"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})
WITH pathway, rle, re.identifier as protein, re.geneName as genes
WHERE size(genes) > 0 AND pathway.stId IN ["R-HSA-110056", "R-HSA-6783783"]
UNWIND genes as gene
WITH DISTINCT pathway, gene, protein
RETURN DISTINCT collect(DISTINCT pathway.stId), gene
~~~~

Proteins:
~~~~
MATCH (pathway:Pathway{speciesName:"Homo sapiens"})-[:hasEvent*]->(rle:Reaction{speciesName:"Homo sapiens"}),
      (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:"Homo sapiens"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})
      WHERE pathway.stId IN ["R-HSA-110056", "R-HSA-6783783"]
RETURN DISTINCT collect(DISTINCT pathway.stId), re.identifier
~~~~

Proteoforms:
~~~~
MATCH (pathway:Pathway{speciesName:"Homo sapiens"})-[:hasEvent*]->(rle:Reaction{speciesName:"Homo sapiens"}),
      (rle)-[:input|output|catalystActivity|physicalEntity|regulatedBy|regulator|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity{speciesName:"Homo sapiens"})-[:referenceEntity]->(re:ReferenceEntity{databaseName:'UniProt'})
      WHERE pathway.stId IN ["R-HSA-109703", "R-HSA-111447"]
WITH DISTINCT pathway, rle, pe, re
OPTIONAL MATCH (pe)-[:hasModifiedResidue]->(tm:TranslationalModification)-[:psiMod]->(mod:PsiMod)
WITH DISTINCT pathway, rle.stId as reaction, pe.stId AS physicalEntity,
                re.identifier AS protein, re.variantIdentifier AS isoform,  tm.coordinate as coordinate, 
                mod.identifier as type 
ORDER BY type, coordinate
WITH DISTINCT pathway, reaction, physicalEntity, protein,
                CASE WHEN isoform IS NOT NULL THEN isoform ELSE protein END as isoform,
                COLLECT(type + ":" + CASE WHEN coordinate IS NOT NULL THEN coordinate ELSE "null" END) AS ptms
RETURN DISTINCT collect(DISTINCT pathway.stId), isoform, ptms
ORDER BY isoform, ptms
~~~~

## Examples of modified overlap

Follow the same steps of the gene level only overlap.

# Disease modules

* Download the data of all GWAS from https://www.ncbi.nlm.nih.gov/gap/phegeni and store it at resources/PheGenI/

