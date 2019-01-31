# Note:

* The cypher files always end with ".cypher".
* They can span over multiple lines.
* Prefer the use of single quotes (`'`). Sometimes there were errors with double quotes (`"`).

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