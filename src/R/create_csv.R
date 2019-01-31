# Function to create csv from a cypher query to Neo4j

source("src/R/paths.R")

create.csv <- function(query, file){
  command <- paste("cypher-shell",
                   "--non-interactive",
                   shQuote(query),
                   " > ",
                   shQuote(file)
  )
  cat("Executing command: ", command)
  shell(command)
}

# Return a dataframe with the data extracted from Neo4j
get.data <- function(name) {
  file.csv <- get.path.csv(name)
  file.cypher <- get.path.cypher(name)
  query <- paste(readLines(file.cypher), collapse=" ")
  create.csv(query, file.csv)
  read.csv(file.csv)
}