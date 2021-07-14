################################################################################
#### How To Rasterize Stuff
################################################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(davidoff)
library(viridis)
library(tidyverse)
library(raster)
library(sf)

# Set seed
set.seed(123)

# Function to create random lines
randomLine <- function(nodes = 2){
  x <- runif(nodes)
  y <- runif(nodes)
  line <- spLines(cbind(x, y))
}

# Create a couple of lines
lines <- lapply(rep(10, 5), randomLine) %>% do.call(rbind, .)

# Add some info
lines$ID <- 1:length(lines)

# Plot the lines
plot(lines, col = viridis(length(lines)))

# Create an empty raster onto which we can rasterize
r <- raster(ncol = 20, nrow = 20, xmn = 0, xmx = 1, ymn = 0, ymx = 1, vals = 0)

# Rasterize lines
lines_r <- rasterizeVelox(lines, r)

# Visualize
plot(lines_r, col = "black", box = F, axes = F, legend = F)
plot(lines, col = "white", add = T, lty = 1:5)

plot(lines_r, col = magma(20), box = F, axes = F, legend = F)
plot(lines, add = T, col = "white", lty = 1:5)
text(lines_r)

# Visualize
png("test.png", width = 1080, height = 1080, pointsize = 30)
plot(lines_r, col = "white", box = F, axes = F, legend = F)
plot(lines, col = "black", add = T, lty = 1:5, lwd = 3)
dev.off()

png("test2.png", width = 1080, height = 1080, pointsize = 30)
plot(lines_r, col = magma(20), box = F, axes = F, legend = F, lwd = 2)
plot(lines, add = T, col = "white", lty = 1:5, lwd = 3)
text(lines_r)
dev.off()
