# Calculate the improvement of search results using protein accessions or proteoforms

library(ggplot2)

source("src/R/paths.R")
source("src/R/create_csv.R")

# ----------------------------------------------------------------------------------------------
# Jitter plot showing the percentage of improvement for each proteoform.

## Get the hits of each protein
hits.proteins <- get.data("reactions_and_pathways_per_protein")

## Get hits of each proteoform with the protein
hits.proteoforms <- get.data("reactions_and_pathways_per_proteoform")

data <- merge(hits.proteins, hits.proteoforms, by = "protein")
colnames(data) <- c("proteins", "reactionCount.protein", "pathwayCount.protein", "isoform", "ptms", "num_ptms", "reactionCount.proteoform", "pathwayCount.proteoform")
data$reactionCount.diff <- 100 - (data$reactionCount.proteoform * 100 / data$reactionCount.protein)
data$pathwayCount.diff <- 100 - (data$pathwayCount.proteoform * 100 / data$pathwayCount.protein)

data$reactionCount.diff[is.nan(data$reactionCount.diff)] <- 0
data$pathwayCount.diff[is.nan(data$pathwayCount.diff)] <- 0

head(data)
nrow(data)

data[which(data$reactionCount.protein== 0), ]

## Plot improvement for proteoforms with n or more modifications
plot.improvement <- function(data, name, type = "reaction", n = 2) {
  ss <- data[which(data$num_ptms >= n),]
  if(type == "reaction"){
    colnames(ss)[9] <- "count.diff"
  } else if(type == "pathway"){
    colnames(ss)[10] <- "count.diff"
  } else{
    stop("Wrong entity type, must be: reaction | pathway")
  }
  p <- ggplot(data = ss, aes(x=count.diff)) + geom_histogram(binwidth=.5, colour="black", fill="white") + 
    geom_vline(data = ss, aes(xintercept=mean(count.diff)),
               linetype="dashed", size=1, colour="red")
  png(get.path.figure(name), height = 12, width = 15, units = "cm", res = 600)
  plot(p)
  dummy <- dev.off()
}

## Plot to showing the average percentage of improvement for proteoforms with 0 or more modifications
plot.improvement(data, "improvement_reactions_0", type = "reaction", n = 0)
plot.improvement(data, "improvement_reactions_1", type = "reaction", n = 1)
plot.improvement(data, "improvement_reactions_2", type = "reaction", n = 2)
plot.improvement(data, "improvement_reactions_3", type = "reaction", n = 3)
plot.improvement(data, "improvement_reactions_4", type = "reaction", n = 4)
plot.improvement(data, "improvement_reactions_5", type = "reaction", n = 5)

plot.improvement(data, "improvement_pathways_0", type = "pathway", n = 0)
plot.improvement(data, "improvement_pathways_1", type = "pathway", n = 1)
plot.improvement(data, "improvement_pathways_2", type = "pathway", n = 2)
plot.improvement(data, "improvement_pathways_3", type = "pathway", n = 3)
plot.improvement(data, "improvement_pathways_4", type = "pathway", n = 4)
plot.improvement(data, "improvement_pathways_5", type = "pathway", n = 5)


## Plot to showing the average percentage of improvement for proteoforms with 1 or more modifications
head(data[which(data$reactionCount.diff < 100 & data$reactionCount.diff > 90),])

