################################################################################
#### Visualization of the Social Landscape
################################################################################
# Clear R's brain
rm(list = ls())

# Set a seed
set.seed(123)

# Load required packages
library(raster)
library(tidyverse)
library(adehabitatHR)
library(viridis)

# Create an empty raster
r <- raster(ncol = 30, nrow = 30, xmn = 0, xmx = 900, ymn = 0, ymx = 900)
values(r) <- rnorm(mean = 0, sd = 1, n = ncell(r))

# Simulate 3 points of attraction
points <- sampleRandom(r, size = 3, sp = T)

# Calculate distance
distance <- distanceFromPoints(r, points)

# Create likelihood surface (with a bit of randomness)
like <- - distance + 75 * r
like <- (like - minValue(like)) / (maxValue(like) - minValue(like))

# Visualize it
plot(like)

# Function to simulate movement
move <- function(cell = NULL, like = NULL, steps = NULL, directions = 8){
  path <- c()
  for (i in 1:steps){
    path[i] <- cell
    neighbors <- adjacent(like, cell, directions = directions, pairs = F)
    values <- values(like)[neighbors]
    probs <- values / sum(values)
    cell <- sample(neighbors, prob = probs, size = 1)
  }
  return(path)
}

# Use it to simulate a couple of tracks
tracks <- tibble(TrackID = 1:5)
tracks$Coords <- lapply(1:5, function(x){
  track <- move(cell = 500, like = like, steps = 500, directions = 16)
  track <- SpatialPoints(xyFromCell(r, cell = track))
  return(track)
})

# Create utilization distribution
r <- extend(r, c(-200, 1100, -200, 1100))
r <- disaggregate(r, fact = 5)
tracks$UD <- lapply(tracks$Coords, function(x){
  kernelUD(x, h = "href", grid = as(r, "SpatialPixels"))
})

# Get 95% Home range
tracks$HR <- lapply(tracks$UD, function(x){
  getverticeshr(x, percent = 95)
})

# Convert UDs to rasters
tracks$UD <- lapply(tracks$UD, raster)

# Plot 1
png("plot1.png", width = 1080, height = 1080, bg = "transparent")
plot(tracks$UD[[1]], col = "transparent", bg = "transparent", box = F, axes = F, legend = F)
plot(spLines(tracks$Coords[[1]]), add = T, lwd = 5, col = "white")
dev.off()

# Plot 2
png("plot2.png", width = 1080, height = 1080, bg = "transparent")
plot(tracks$UD[[1]], col = magma(20), bg = "transparent", box = F, axes = F, legend = F)
dev.off()

# Plot 3
png("plot3.png", width = 1080, height = 1080, bg = "transparent")
plot(tracks$UD[[1]], col = "transparent", bg = "transparent", box = F, axes = F, legend = F)
plot(tracks$HR[[1]], add = T, border = "white", lty = 2, lwd = 5)
dev.off()

# Plot 4
png("plot4.png", width = 1080, height = 1080, bg = "transparent")
plot(tracks$UD[[2]], col = "transparent", bg = "transparent", box = F, axes = F, legend = F)
plot(spLines(tracks$Coords[[2]]), add = T, lwd = 5, col = "white")
dev.off()

# Plot 5
png("plot5.png", width = 1080, height = 1080, bg = "transparent")
plot(tracks$UD[[2]], col = magma(20), bg = "transparent", box = F, axes = F, legend = F)
dev.off()

# Plot 6
png("plot6.png", width = 1080, height = 1080, bg = "transparent")
plot(tracks$UD[[2]], col = "transparent", bg = "transparent", box = F, axes = F, legend = F)
plot(tracks$HR[[2]], add = T, border = "white", lty = 2, lwd = 5)
dev.off()
