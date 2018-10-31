# 
# This script demonstrates some functions to manipulate graphs using igraph.
#

# Libraries

library(igraph)


# Parameters

## Resources Files

intactEdgesFile <- "resources/intact/intact_18.08.17_edges"
reactomeEdgesFile <- "resources/reactome/proteinInternalEdges.tsv.gz"


# Functions

#' Merges the reactome and intact graphs.
#' 
#' @param reactomeGraph a vertice of the intact edge
#' @param verticeIntact2 the other vertice of the intact edge
#' @param edgesReactomeCommon the Reactome edge list
#' 
#' @return a boolean indicating whether the given intact edge is in Reactome
mergeGraphs <- function(reactomeGraph, intactGraph) {
    
    result <- union(reactomeGraph, intactGraph)
    
    if ("source_1" %in% list.edge.attributes(result)) {
        
        source1 <- get.edge.attribute(result, "source_1")
        source2 <- get.edge.attribute(result, "source_2")
        sourceMerge <- ifelse(!is.na(source1), source1, source2)
        
        result <- set.edge.attribute(graph = reactomeGraph, name = "source", value = sourceMerge)
        
        result <- remove.edge.attribute(result, "source_1")
        result <- remove.edge.attribute(result, "source_2")
        
    }
    
    return(result)
}


#' Returns the main component of a graph.
#' 
#' @param graph the entire graph
#' 
#' @return the main component of the graph
getMainComponent <- function(graph) {
    
    components <- decompose(graph)
    mainComponent <- NULL
    maxSize <- 0
    
    for (subGraph in components) {
        
        size <- length(V(subGraph))
        
        if (size > maxSize) {
            
            mainComponent <- subGraph
            maxSize <- size
            
        }
    }
    
    return(mainComponent)
    
}


# Main script

## Load tables as data frames

print(paste(Sys.time(), " Loading data", sep = ""))

edgesIntact <- read.table(intactEdgesFile, header = T, sep = " ", stringsAsFactors = F, quote = "", comment.char = "")

edgesReactome <- read.table(reactomeEdgesFile, header = T, sep = "\t", quote = "", comment.char = "", stringsAsFactors = F)
edgesReactome <- edgesReactome[, c("id1", "id2")]


## Make graphs

print(paste(Sys.time(), " Making graphs", sep = ""))

graphReactome <- graph_from_data_frame(edgesReactome)
graphIntact <- graph_from_data_frame(edgesIntact)


## Simplify graphs

print(paste(Sys.time(), " Simplifying graphs", sep = ""))

graphReactome <- simplify(graphReactome, remove.multiple = T, remove.loops = T, edge.attr.comb = "first")
graphIntact <- simplify(graphIntact, remove.multiple = T, remove.loops = T, edge.attr.comb = "first")


## Merge Reactome and Intact

print(paste(Sys.time(), " Merging Reactome and Intact", sep = ""))

graphMerged <- mergeGraphs(graphReactome, graphIntact)


## Extract main component

print(paste(Sys.time(), " Extracting main component", sep = ""))

graphReactome <- getMainComponent(graphReactome)
graphIntact <- getMainComponent(graphIntact)
graphMerged <- getMainComponent(graphMerged)


## Get degrees

print(paste(Sys.time(), " Getting degrees", sep = ""))

degreeReactome <- degree(graphReactome)
degreeIntact <- degree(graphIntact)
degreeMerged <- degree(graphMerged)


## Set weights according to the degree

print(paste(Sys.time(), " Setting weights", sep = ""))

reactomeDF <- as_data_frame(graphReactome)
reactomeEdgesDegree <- degreeReactome[reactomeDF$from] + degreeReactome[reactomeDF$to] - 1
reactomeP <- 1 / reactomeEdgesDegree
reactomePLog <- -log10(reactomeP)
E(graphReactome)$weight <- reactomePLog

intactDF <- as_data_frame(graphIntact)
intactEdgesDegree <- degreeIntact[intactDF$from] + degreeIntact[intactDF$to] - 1
intactP <- 1 / intactEdgesDegree
intactPLog <- -log10(intactP)
E(graphIntact)$weight <- intactPLog

mergedDF <- as_data_frame(graphMerged)
mergedEdgesDegree <- degreeMerged[mergedDF$from] + degreeMerged[mergedDF$to] - 1
mergedP <- 1 / mergedEdgesDegree
mergedPLog <- -log10(mergedP)
E(graphMerged)$weight <- mergedPLog


## Get random sets of proteins and make a matrix of the distance

sampleSize <- 10

graphList <- list(graphReactome, graphIntact, graphMerged)

for (graph in graphList) {
    
    simpleDistanceMatrix <- matrix(data = NA, nrow = sampleSize, ncol = sampleSize)
    weightedDistanceMatrix <- matrix(data = NA, nrow = sampleSize, ncol = sampleSize)
    
    proteinsI <- sample(V(graph)$name, size = sampleSize)
    proteinsJ <- sample(V(graph)$name, size = sampleSize)
    
    for (i in 1:sampleSize) {
        
        for (j in 1:sampleSize) {
            
            simpleDistanceMatrix[i, j] <- distances(graph = graph, v = proteinsI[i], to = proteinsJ[j], weights = NA)
            weightedDistanceMatrix[i, j] <- distances(graph = graph, v = proteinsI[i], to = proteinsJ[j], weights = NULL)
            
        }
    }
    
    plot(simpleDistanceMatrix, weightedDistanceMatrix)
    
}