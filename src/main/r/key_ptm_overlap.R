# Read files with the data on the pathways overlap with only modified proteoforms

library(ggplot2)

path <- "resources/reactome/key_ptm_overlap/"

file.proteins <- paste(path, "proteins.txt", sep = "")
file.proteoforms <- paste(path, "proteoforms.txt", sep = "")
file.modifications <- paste(path, "modifications.txt", sep = "")
#file.pathways <- paste(path, "pathway_pairs.txt", sep = "")

proteins <- read.csv(file = file.proteins, sep = "\t")
proteoforms <- read.csv(file = file.proteoforms, sep = "\t")
modifications <- read.csv(file = file.modifications, sep = "\t")
#pathways <- read.csv(file = file.pathways, sep = "\t")

plot.proteins <- ggplot(proteins, aes)
