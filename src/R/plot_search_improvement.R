# Calculate the improvement of search results using protein accessions or proteoforms

library(ggplot2)

source("src/R/paths.R")
source("src/R/create_csv.R")

# ----------------------------------------------------------------------------------------------
# Jitter plot showing the percentage of improvement for each proteoform.



  # Plot to showing the average percentage of improvement for proteoforms with 0, 1, 2,... n modifications compared
# to searching using the pure accession number.
