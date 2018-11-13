#!/usr/bin/env Rscript

# Read file with the modification frequency of the set overlap

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

stopifnot(length(args) >= 4)

file.pathways <- args[1]
file.modifications <- args[2]
file.proteins <- args[3]
file.proteoforms <- args[4]

proteins <- read.csv(file = file.proteins, sep = "\t")
proteoforms <- read.csv(file = file.proteoforms, sep = "\t")
modifications <- read.csv(file = file.modifications, sep = "\t")
pathways <- read.csv(file = file.pathways, sep = "\t")

plotHistogram <- function(data, title) {
  ggplot(data = data, aes(FREQUENCY)) + geom_histogram(binwidth = 20) + ggtitle(title)
}

png('figures/4_modified_overlap/protein_frequencies.png')
plot.proteins <- plotHistogram(proteins, "Protein frequencies")
plot.proteins

png('figures/4_modified_overlap/proteoforms_frequencies.png')
plot.proteoforms <- plotHistogram(proteoforms, "Proteoform frequencies")
plot.proteoforms

png('figures/4_modified_overlap/modifications_frequencies.png')
plot.modifications <- plotHistogram(modifications, "Modification frequencies")
plot.modifications
