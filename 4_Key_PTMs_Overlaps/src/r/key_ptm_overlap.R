# Read files with the data on the pathways overlap with only modified proteoforms

library(ggplot2)

path <- "resources/reactome/key_ptm_overlap/"

file.proteins <- paste(path, "proteins.txt", sep = "")
file.proteoforms <- paste(path, "proteoforms.txt", sep = "")
file.modifications <- paste(path, "modifications.txt", sep = "")
file.pathways <- paste(path, "pathway_pairs.txt", sep = "")

proteins <- read.csv(file = file.proteins, sep = "\t")
proteoforms <- read.csv(file = file.proteoforms, sep = "\t")
modifications <- read.csv(file = file.modifications, sep = "\t")
pathways <- read.csv(file = file.pathways, sep = "\t")

plotHistogram <- function(data, title) {
  ggplot(data = data, aes(FREQUENCY)) + geom_histogram(binwidth = 20) + ggtitle(title)
}

plot.proteins <- plotHistogram(proteins, "Protein frequencies")
plot.proteins
plot.proteoforms <- plotHistogram(proteoforms, "Proteoform frequencies")
plot.proteoforms
plot.modifications <- plotHistogram(modifiactions, "Modification frequencies")
plot.modifications

pathways.80 <- pathways[which(pathways$MODIFIED_RATIO == 1),]

head(proteoforms[order(proteoforms$FREQ, decreasing = TRUE),])

names(proteoforms) <- c("PROTEOFORM", "FREQ")

proteoforms[which(proteoforms$FREQUENCY == "45104"),]
