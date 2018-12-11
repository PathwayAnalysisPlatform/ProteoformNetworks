# Plot the number of reactions and pathways where each gene, protein and proteoform participates

set.seed(1234)
hits = data.frame(
  entity.type = factor(rep(c("GENE", "PROTEIN", "PROTEOFORM"), each=200)),
  reactions = sample(1:300, 200),
  pathways = sample(1:200, 200)
)
head(hits, 4)