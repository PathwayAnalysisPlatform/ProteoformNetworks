
# Read GWAS table with all the phenotypes ----

path <- "../resources/PheGenI/"
phenotypes.csv <- read.csv(file = "../resources/PheGenI/PheGenI_Association_full.tab", sep = "\t")
phenotypes.csv$Trait <- as.character(phenotypes.csv$Trait)
phenotypes.csv$Trait <- gsub(",", "",  gsub(" ", "", gsub("-", "_", phenotypes.csv$Trait)))
phenotypes.csv$Gene <- as.character(phenotypes.csv$Gene)
phenotypes.csv$Gene.2 <- as.character(phenotypes.csv$Gene.2)

# Create gene list files for each trait ----
# store first record

trait <- phenotypes.csv$Trait[1]
gene1 <- phenotypes.csv$Gene[1]
gene2 <- phenotypes.csv$Gene.2[1]

gene.set <- c(gene1)
if(gene1 != gene2) {
  gene.set <- c(gene.set, gene2)
}
trait.previous <- trait

# store the rest of records
trait <- "nothing"
for(r in 1:length(phenotypes.csv$Trait)) {
  trait <- phenotypes.csv$Trait[r]
  gene1 <- phenotypes.csv$Gene[r]
  gene2 <- phenotypes.csv$Gene.2[r]
  if(trait.previous != trait) {
    cat("Storing: ", trait.previous, "\n")
    write.csv(unique(gene.set), paste(path, trait.previous, ".txt", sep = ""), quote = FALSE, row.names = FALSE)
    gene.set <- c()
  }
  
  gene.set <- c(gene.set, gene1)
  if(gene1 != gene2) {
    gene.set <- c(gene.set, gene2)
  }
  
  trait.previous <- trait
}
# store the last trait
cat("Storing: ", trait, "\n")
write.csv(unique(gene.set), paste(path, trait, ".txt", sep = ""), quote = FALSE, row.names = FALSE)

# Create an igraph network for each phenotype ----
## Execute PathwayMatcher to get the graph formed by internal edges in proteins and proteoforms

traits <- sort(unique(phenotypes.csv$Trait))

for(t in traits){
  cat("Creating networks for: ", t, "\n")
  command <- paste("java", "-jar", "../resources/PathwayMatcher.jar", 
                   "-i", paste(path, t, ".txt", sep = ""),
                   "-o", paste(path, t, sep = ""),
                   "-t", "gene",
                   "-gg", "-gp", "-gu", "-tlp"
              )
  cat(command, "\n")
  system(command, intern = TRUE)  
}

# Compare each pair of phenotypes for overlap in the protein and proteoform ----

# Create three matrices with the number of overlaping genes, proteins, proteoforms and mod
InitMatrix <- function() {
  m <- matrix(0L, nrow = length(traits), ncol = length(traits))
  row.names(m) <- traits
  colnames(m) <- traits
  m
}

m.g <- InitMatrix() # Matrix for gene overlap
m.pt <- InitMatrix() # Matrix for protein overlap
m.pf <- InitMatrix() # Matrix for proteoforms overlap
m.mp <- InitMatrix() # Matrix for modified proteoforms overlap

## Calculate the overlapping

### Create functions to read the graphs
library(igraph)

LoadGraph <- function(entity) {
  function(trait) {
    print(paste(" Loading ", entity," data for: ", trait, sep = ""))
    file <- paste(path, trait, "/", entity, "InternalEdges.tsv", sep = "")
    table <- read.table(file, header = T, sep = "\t", quote = "", comment.char = "", stringsAsFactors = F)
    graph_from_data_frame(table)  
  }
}

LoadGraphGene <- LoadGraph("gene")
LoadGraphProtein <- LoadGraph("protein")
LoadGraphProteoform <- LoadGraph("proteoform")

### Create functions to calculate graph intersections

Overlap <- function(g1, g2) {
  return(length(intersect(V(g1)$name, V(g2)$name)))
}

### Traverse each ordered pair

set.candidate.gene <- c()
set.candidate.protein <- c()
set.candidate.proteoform <- c()

hnf1a <- c(45, 207, 279, 370)

for(t1 in traits) {
  
  cat("----------")
    
  t1.graph.genes <- LoadGraphGene(t1)
  t1.graph.protein <- LoadGraphProtein(t1)
  t1.graph.proteoform <- LoadGraphProteoform(t1)
  
  for(t2 in traits) {
    
    if(t1 < t2) {
      cat("Comparing ", t1, " with ", t2, "\n")
      
      t2.graph.genes <- LoadGraphGene(t2)
      t2.graph.protein <- LoadGraphProtein(t2)
      t2.graph.proteoform <- LoadGraphProteoform(t2)

      m.g[t1, t2] <- Overlap(t1.graph.genes, t2.graph.genes)
      #cat(t1, " -- ", m.g[t1, t2], ": " , m.g[t1, t2], "\n")
      m.pt[t1, t2] <- Overlap(t1.graph.protein, t2.graph.protein)
      m.pf[t1, t2] <- Overlap(t1.graph.proteoform, t2.graph.proteoform)
      
      ## If the overlapping entities are proteoforms add to the candidates
      ## If they overlap just in the protein or the proteoform network
    }
  }
}
