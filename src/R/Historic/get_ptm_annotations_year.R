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

name <- "ptm_annotations_type_and_date"
data <- get.data(name)
data$date <- as.Date(data$date, format="%Y")
data$order <- 1:nrow(data)

# png(get.path.figure(name), width = 800, height = 600)
p <- ggplot(data, aes(date, order))
p + geom_jitter(aes(colour = psimod), width = 2000)
# dev.off()
p

# ----------------------------------------------------------------------------
# Phosphorylated
name <- "phospho_annotations"
data <- get.data(name)
data$date <- as.Date(data$date, format="%Y")
data$order <- 1:nrow(data)
data$psimod <- as.factor(data$psimod)

png(get.path.figure(name), width = 800, height = 600)
p <- ggplot(data, aes(date, order))
p + geom_jitter(aes(colour = psimod), width = 500)
dev.off()
p

# ----------------------------------------------------------------------------
# Hydroxylated
name <- "hydro_annotations"
data <- get.data(name)
data$date <- as.Date(data$date, format="%Y")
data$order <- 1:nrow(data)
data$psimod <- as.factor(data$psimod)

png(get.path.figure(name), width = 800, height = 600)
p <- ggplot(data, aes(date, order))
p + geom_jitter(aes(colour = psimod), width = 500)
dev.off()
p

# ----------------------------------------------------------------------------
# Glycosylated
name <- "glyco_annotations"
data <- get.data(name)
data$date <- as.Date(data$date, format="%Y")
data$order <- 1:nrow(data)
data$psimod <- as.factor(data$psimod)

png(get.path.figure(name), width = 800, height = 600)
p <- ggplot(data, aes(date, order))
p + geom_jitter(aes(colour = psimod), width = 500)
dev.off()
p

