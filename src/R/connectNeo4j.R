# Set up connection to Neo4j

library(neo4r)

con <- neo4j_api$new(
  url = "http://localhost:7474",
  user = "", 
  password = ""
)
status <- con$ping()

if(status != 200){
  stop("Not connected to Neo4j")
}