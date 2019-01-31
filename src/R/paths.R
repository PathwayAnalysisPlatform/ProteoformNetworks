# Get the complete paths for each type of file

path.figure <- "figures/"
path.r <- "src/R/"
path.cypher <- "src/Cypher/"
path.csv <- "resources/csv/"

get.path.figure <- function(name) {
  paste0(path.figure, name, ".png")
}

get.path.r <- function(name) {
  paste0(path.R, name, ".R")
}

get.path.cypher <- function(name) {
  paste0(path.cypher, name, ".cypher")
}

get.path.csv <- function(name) {
  paste0(path.csv, name, ".csv")
}