################################################################################
#### Random Network Plot
################################################################################
# Description: Simple plot of a random network that I can put in the background
# of a presentation

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(igraph)

# Create random network
net <- erdos.renyi.game(n = 12, 20, type = "gnm")

# Visualize it
plot(net, vertex.label = NA)
