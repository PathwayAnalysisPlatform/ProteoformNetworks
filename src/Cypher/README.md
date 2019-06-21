# Note:

* The cypher files always end with ".cypher".
* They can span over multiple lines.
* Prefer the use of single quotes (`'`). Sometimes there were errors with double quotes (`"`).
* Avoid the `|` character. It breaks the command string.
The error looks like this, where "gly" is the word right after the `|`:
~~~~
'gly' is not recognized as an internal or external command,
operable program or batch file.
Warning message:
In shell(command) :
  'cypher-shell --non-interactive "MATCH (pe:PhysicalEntity{speciesName:\"Homo sapiens\"})-[:hasModifiedResidue]->(tm:TranslationalModification)-[:created]-(tmie:InstanceEdit), (psimod:PsiMod)<-[:psiMod]-(tm)  WHERE psimod.displayName =~ \".*(gluc|gly|fucosyl|xylosyl|galactos).*\" OR  RETURN pe.displayName as physicalEntity, psimod.identifier as psimod, tm.displayName as modification, tmie.dateTime as date  ORDER BY date, psimod"  >  "resources/csv/glyco_annotations.csv"' execution failed with error code 255
~~~~
* The last line should be an empty line.
Otherwise, the error looks like this:
~~~~
In readLines(file.cypher) :
  incomplete final line found on 'src/Cypher/phospho_annotations.cypher'
~~~~
* If a cypher command has a mistake, the error could look like:
~~~~
                                                         ^
Error in read.table(file = file, header = header, sep = sep, quote = quote,  : 
  no lines available in input
In addition: Warning message:
In shell(command) :
~~~~