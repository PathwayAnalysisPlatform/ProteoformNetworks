degree.proteins <- read.csv("reports/node_degree_proteins.txt", stringsAsFactors = F, header = T, sep = "\t")
degree.mean.proteins <- mean(degree.proteins[,2])

degree.proteoforms <- read.csv("reports/node_degree_proteoforms.txt", stringsAsFactors = F, header = T, sep = "\t")
degree.mean.proteoforms <- mean(degree.proteoforms[,2])
