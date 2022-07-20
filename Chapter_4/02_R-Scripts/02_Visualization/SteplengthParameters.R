################################################################################
#### Step Length Parameters
################################################################################
# Description: Plot of step length distribution parameters

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)     # For data wrangling

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load fitted distributions
dists_dynamic <- read_rds("03_Data/03_Results/StepLengthDynamic.rds")
dists_means <- read_rds("03_Data/03_Results/StepLengthMeans.rds")

# Function to make nice facet labels
labelfun <- function(l) {
  rl <- round(l, 3)
  l <- ifelse(rl == 0 & l != 0, format(l, scientific = T), l)
}

# Visualize everything
p <- dists_dynamic %>%
  pivot_longer(shape:mu, names_to = "Parameter", values_to = "Value") %>%
  ggplot(aes(x = duration, y = Value)) +
    geom_jitter(width = 0.1, alpha = 0.2, size = 0.5) +
    geom_point(data = dists_means, aes(x = duration, y = mean), col = "orange", size = 5) +
    geom_line(data = dists_means, aes(x = duration, y = mean), col = "orange") +
    scale_y_continuous(labels = labelfun) +
    facet_wrap(~ Parameter, nrow = 2, scales = "free") +
    theme_minimal() +
    xlab("Step Duration") +
    ylab("Parameter Estimate") +
    theme(strip.background = element_rect(fill = "gray95", color = NA))

# Store the plot
ggsave("04_Manuscript/99_DistributionParameters.png"
  , plot   = p
  , width  = 5
  , height = 5
  , scale  = 1.4
  , bg     = "white"
)
