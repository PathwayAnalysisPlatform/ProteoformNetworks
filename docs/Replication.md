# Steps to execute

## Examples of pathways

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

* Execute C++ program to create the pathway sets and calculate overlaps:
~~~~

~~~~