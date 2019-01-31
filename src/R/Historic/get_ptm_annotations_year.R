library(ggplot2)

source("src/R/paths.R")
source("src/R/create_csv.R")

# Read data and correct format
name <- "ptm_annotations_type_and_date_aggregated"
data <- get.data(name)
data$date <- as.Date(data$date, format="%Y")

# Prepare and select the data
data$total <- cumsum(data$times)
data <- aggregate(x = data$times, by = list(data$date), FUN = sum)
colnames(data) <- c("date", "annotations")
data$total <- cumsum(data$annotations)

# Prepare the plot
png(get.path.figure(name), width = 800, height = 600)
y.ticks <- seq(0, data$total[nrow(data)], length = 10)
p <- ggplot(data = data, aes(x = date, y = total)) +
  geom_line(color = "#9ecae1", size = 2) +
  geom_point(color = "#3182bd", size = 2) +
  labs(title = "Evolution of PTM annotations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(name = "Years", limits = data$date, labels = as.numeric(format(data$date,'%Y'))) +
  scale_y_discrete(name = "# PTM annotations", limits = y.ticks, labels = y.ticks)
dev.off()
p

# ----------------------------------------------------------------------------

# Jitter plot

plot.ptm.annotations <- function(name, title = "PTM Annotation Evolution"){
  # Get and process the data
  data <- get.data(name)
  data$date <- as.Date(data$date, format="%Y")
  data$order <- 1:nrow(data)
  data$psimod <- as.factor(data$psimod)
  
  # Make the plot
  y.ticks <- as.integer(seq(0, data$order[nrow(data)], length = 10))
  p <- ggplot(data, aes(date, order)) +
    labs(title = title) +
    geom_jitter(aes(colour = psimod), width = 500) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_x_discrete(name = "Years", limits = data$date, labels = as.numeric(format(data$date,'%Y'))) +
    scale_y_discrete(name = "# PTM annotations", limits = y.ticks, labels = y.ticks)
  png(get.path.figure(name), height = 12, width = 12, units = "cm", res = 600)
  plot(p)
  dummy <- dev.off()
  return(p)
}

# ----------------------------------------------------------------------------
# All
p <- plot.ptm.annotations("ptm_annotations_type_and_date")
p

# ----------------------------------------------------------------------------
# Phosphorylated
p <- plot.ptm.annotations("phospho_annotations")
p

# ----------------------------------------------------------------------------
# Hydroxylated
p <- plot.ptm.annotations("hydro_annotations")
p

# ----------------------------------------------------------------------------
# Glycosylated
p <- plot.ptm.annotations("glyco_annotations")
p
