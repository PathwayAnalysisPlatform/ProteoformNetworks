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
y.ticks <- seq(0, data$total[nrow(data)], length = 10)
p <- ggplot(data = data, aes(x = date, y = total)) +
  geom_line(color = "#9ecae1", size = 2) +
  geom_point(color = "#3182bd", size = 2) +
  labs(title = "Evolution of PTM annotations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(name = "Years", limits = data$date, labels = as.numeric(format(data$date,'%Y'))) +
  ylab("# PTM annotations")
png(get.path.figure(name), height = 12, width = 15, units = "cm", res = 600)
plot(p)
dummy <- dev.off()
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
  png(get.path.figure(name), height = 12, width = 20, units = "cm", res = 600)
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

# ----------------------------------------------------------------------------

data <- get.data("ptm_type_frequency")
data$psimod <- sprintf("%05d", data$psimod)

data$psimod <- factor(data$psimod, levels = data$psimod[order(data$times)])
data$psimod  # notice the changed order of factor levels

data$group <- "other"

data[grepl("phospho", data$name),]$group <- "phosphorylated"

data[grepl("gly", data$name),]$group <- "glycosylated"
data[grepl("gluc", data$name),]$group <- "glycosylated"
data[grepl("fucosyl", data$name),]$group <- "glycosylated"
data[grepl("xylosyl", data$name),]$group <- "glycosylated"
data[grepl("galactos", data$name),]$group <- "glycosylated"

data[grepl("hydro", data$name),]$group <- "hydroxylated"

data[data$psimod %in% c("01148", "01149", "00064", "01381", "00115"),]$group <- "acylated"

data[data$psimod == "00041",]$group <- "carboxylated"

data[data$psimod %in% c("00084", "00083", "00085", "00342", "00599", "00076", "00078", "00077", "00239", "00080"),]$group <- "methylated"

y.ticks <- as.integer(seq(0, max(data$times), length = 10))

p <- ggplot(data[1:50,], aes(psimod, times)) + 
  geom_point(aes(colour = group)) +
  labs(title = "PTM Frequencies") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(name = "# PTM annotations", limits = y.ticks, labels = y.ticks)
png(get.path.figure("ptm_type_frequency"), height = 12, width = 20, units = "cm", res = 600)
plot(p)
dummy <- dev.off()
p

# ----------------------------------------------------------------------------
# PTM Frequencies by group type aggregation

data <- aggregate(x = data$times, by = list(data$group), FUN = sum)
colnames(data) <- c("group", "annotations")

data$group <- factor(data$group, levels = data$group[order(data$annotations)])

y.ticks <- as.integer(data$annotations)

p <- ggplot(data, aes(group, annotations, color = group)) + 
  geom_bar(stat="identity") +
  labs(title = "PTM Frequencies") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_discrete(name = "# PTM annotations", limits = y.ticks, labels = y.ticks)
png(get.path.figure("ptm_type_frequency_by_group"), height = 12, width = 12, units = "cm", res = 600)
plot(p)
dummy <- dev.off()
p
