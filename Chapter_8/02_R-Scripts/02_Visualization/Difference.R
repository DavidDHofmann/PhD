################################################################################
#### Difference Plots
################################################################################
# Description: Put together the individual difference plots

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
setwd(wd)

# Load required packages
library(tidyverse)      # To wrangle data
library(ggpubr)         # To arrange multiple plots

# Load individual plots
p1 <- read_rds("04_Manuscript/99_HeatmapsDifference.rds")
p2 <- read_rds("04_Manuscript/99_BetweennessDifference.rds")
p3 <- read_rds("04_Manuscript/99_DistanceDifference.rds")

# We want to get a similar look to the "facets"
p1 <- p1 + facet_wrap(~ "Heatmap")
p2 <- p2 + facet_wrap(~ "Betweenness")
p3 <- p3 + facet_wrap(~ "Distance")

# # Extract legends
# p1_legend <- get_legend(p1)
# p2_legend <- get_legend(p2)
# p3_legend <- get_legend(p3)
#
# # Remove legends from plots
# p1 <- p1 + theme(legend.position = "none")
# p2 <- p2 + theme(legend.position = "none")
# p3 <- p3 + theme(legend.position = "none")

# Also remove the axis from the second and third plot
p2 <- p2 + rremove("y.text")
p3 <- p3 + rremove("y.text")

# Adjust margins so plots can be put closer together
p1 <- p1 + theme(plot.margin = unit(c(0.1, 0, 0.1, 0), "cm"))
p2 <- p2 + theme(plot.margin = unit(c(0.1, 0, 0.1, -0.75), "cm"))
p3 <- p3 + theme(plot.margin = unit(c(0.1, 0, 0.1, -0.75), "cm"))

# Arrange them and add legends pack
p <- ggarrange(p1, p2, p3, nrow = 1, align = "hv")
ggsave("04_Manuscript/99_Difference.png", p, width = 10, height = 5)
