# Plot protein interaction network evolution by PTM annotation years

library(igraph)

source("src/R/paths.R")
source("src/R/create_csv.R")

# Get all possible connections: pairs of physical entities
connections <- get.data("physical_entity_interactions")
connections$pe1 <- trimws(connections$pe1)
connections$pe2 <- trimws(connections$pe2)
connections$pe1_date <- as.Date(connections$pe1_date, format="%Y")
connections$pe2_date <- as.Date(connections$pe2_date, format="%Y")
winDialogString("Finished getting the connections", default = "Finished getting the connections")

# Get all possible nodes with Date of creation of a proteoform as date of creation of its last modification or itself
nodes <- get.data("proteoforms_creation_date")
nodes$pe <- trimws(nodes$pe)
nodes$date <- as.Date(nodes$date, format="%Y")
winDialogString("Finished getting the nodes", default = "Finished getting the nodes")

# Plot interaction network before date
# Requires: nodes$date format in %Y
# Requires: connections$pe1.stId, connections$pe2.stId
date <- "2003-01-01"
plot.network.by.date <- function(date = "2003-01-01", nodes, connections){
  # Filter nodes
  pe.set <- nodes[which(nodes$date < date),"pe"]
  pe.set <- unique(pe.set)
  
  # Filter connections 
  # pe.set <- c ("R-HSA-166025", "R-HSA-168104", "R-HSA-168108", "R-HSA-8863966", "R-HSA-166033")
  good1 <- which(connections$pe1.stId %in% pe.set)
  good2 <- which(connections$pe2.stId %in% pe.set)
  good <- good1[good1 %in% good2]
  
  # Pairs of interacting proteoforms as physical entities
  selected <- connections[good, ]
  colnames(selected) <- c("from", "to")
  
  # colnames(selected) <- c("from", "pe")
  # selected <- merge(selected, nodes, by="pe")
  
  net <- graph_from_data_frame( selected, directed=F) 
  p <- plot(net, edge.arrow.size = 0.5, 
       vertex.color = "gold", 
       vertex.size = 15, 
       vertex.frame.color = "gray",
       vertex.label.color = "black",
       vertex.label.cex = 0.8,
       edge.curved=.3)
  name <- paste0("proteoform_interaction_network_", date)
  # png(get.path.figure(name), height = 12, width = 15, units = "cm", res = 600)
  # plot(p)
  # dummy <- dev.off()
}

# -----------------------------------------------------------------------------

plot.2003.06 <- plot.network.by.date("2003-06-01", nodes, connections)
plot.2003.07 <- plot.network.by.date("2003-07-01", nodes, connections)
plot.2003.08 <- plot.network.by.date("2003-08-01", nodes, connections)
plot.2003.09 <- plot.network.by.date("2003-09-01", nodes, connections)
# plot.2009 <- plot.network.by.date("2009-01-01", nodes, connections)
# plot.2012 <- plot.network.by.date("2012-01-01", nodes, connections)
# plot.2015 <- plot.network.by.date("2015-01-01", nodes, connections)
# plot.2018 <- plot.network.by.date("2018-01-01", nodes, connections)

head(nodes)




