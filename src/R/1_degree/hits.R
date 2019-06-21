# Plot one continuous varible

library(dplyr)
library(ggplot2)
library(ggpubr)

set.seed(1234)
hits = data.frame(
  entity.type = factor(rep(c("GENE", "PROTEIN", "PROTEOFORM"), each=200)),
  reactions = c(rnorm(200, 75), rnorm(200, 64), rnorm(200, 55)),
  pathways = c(rnorm(200, 70), rnorm(200, 60), rnorm(200, 50))
)

hits[sample(1:length(hits$entity.type), 50, replace = F), ]

m <- hits %>%
  group_by(entity.type) %>%
  summarise(grp.reactions.mean = mean(reactions), grp.pathways.mean = mean(pathways))
m

colors.soft = c("#fff7bc", "#f0f0f0", "#deebf7");
colors.medium = c("#ffeda0", "#bdbdbd", "#9ecae1");
colors.strong = c("#fec44f", "#636363", "#3182bd");

# Density plot
p <- ggdensity(hits, x = "reactions",
               add = "mean", rug = TRUE,
               color = "entity.type", fill = "entity.type",
               palette = colors.strong)
p
png(paste0("figures/num_reactions_density.png"), height = 12, width = 12, units = "cm", res = 600)
plot(p)
dummy <- dev.off()

p <- gghistogram(hits, x = "reactions",
               add = "mean", rug = TRUE,
               color = "entity.type", fill = "entity.type",
               palette = colors.strong)
p
png(paste0("figures/num_reactions_count.png"), height = 12, width = 12, units = "cm", res = 600)
plot(p)
dummy <- dev.off()

p <- ggboxplot(hits, x = "entity.type", y = "reactions",
               add = "mean", rug = TRUE,
               color = "entity.type", fill = "entity.type",
               palette = colors.strong)
p
png(paste0("figures/num_reactions_box.png"), height = 12, width = 12, units = "cm", res = 600)
plot(p)
dummy <- dev.off()

p <- ggdensity(hits, x = "pathways",
               add = "mean", rug = TRUE,
               color = "entity.type", fill = "entity.type",
               palette = colors.strong)
p
png(paste0("figures/num_pathways_density.png"), height = 12, width = 12, units = "cm", res = 600)
plot(p)
dummy <- dev.off()

p <- gghistogram(hits, x = "pathways",
                 add = "mean", rug = TRUE,
                 color = "entity.type", fill = "entity.type",
                 palette = colors.strong)
p
png(paste0("figures/num_pathways_count.png"), height = 12, width = 12, units = "cm", res = 600)
plot(p)
dummy <- dev.off()

p <- ggboxplot(hits, x = "entity.type", y = "pathways",
               add = "mean", rug = TRUE,
               color = "entity.type", fill = "entity.type",
               palette = colors.strong)
p
png(paste0("figures/num_pathways_box.png"), height = 12, width = 12, units = "cm", res = 600)
plot(p)
dummy <- dev.off()

