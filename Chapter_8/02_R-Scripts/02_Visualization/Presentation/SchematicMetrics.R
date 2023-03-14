################################################################################
#### Schematic Visualization of the Different Metrics
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
setwd(wd)

# Load required packages
library(igraph)
library(raster)
library(tidyverse)
library(sf)
library(rgeos)
library(Rcpp)
library(pbmcapply)
library(spatstat)
library(maptools)

# Load custom functions
source("02_R-Scripts/00_Functions.R")
sourceCpp("02_R-Scripts/00_Functions.cpp")

################################################################################
#### "Simulated" Movement
################################################################################
# Generate raster
n <- 100
r <- raster(ncol = n, nrow = n, vals = 1:(n ** 2), xmn = 0, xmx = n, ymn = 0, ymx = n)

# Create two circles
circle1 <- gBuffer(SpatialPoints(matrix(c(20, 50), ncol = 2)), width = 10, quadsegs = 50)
circle2 <- gBuffer(SpatialPoints(matrix(c(80, 50), ncol = 2)), width = 10, quadsegs = 50)
circle1 <- as(circle1, "SpatialPolygonsDataFrame")
circle2 <- as(circle2, "SpatialPolygonsDataFrame")
circles <- rbind(circle1, circle2)

# Let's also suppose there is a settlement
pts <- rbind(c(45, 65), c(60, 35))
settl <- SpatialPoints(pts)
settl <- distanceFromPoints(r, settl)
pts <- as.data.frame(pts)
names(pts) <- c("x", "y")

# Create corridor
corr_x <- seq(20, 80, length.out = 100)
corr_y <- sin(seq(0, 2 * pi, length.out = 100)) * 10 + 50
corr <- spLines(cbind(corr_x, corr_y))
corr <- gBuffer(corr, width = 5)

# Plot all
plot(settl)
plot(circle1, add = T)
plot(circle2, add = T)
plot(corr, add = T)
points(pts, col = "red")

# Create paths
createPath <- function() {
  points1 <- spsample(circle1, 20, type = "random")
  points2 <- SpatialPoints(arrange(as.data.frame(coordinates(spsample(corr, 20, type = "random"))), x))
  points3 <- spsample(circle2, 20, type = "random")
  xy <- rbind(coordinates(points1), coordinates(points2), coordinates(points3))
  xy <- data.frame(xy)
  names(xy) <- c("x", "y")
  return(xy)
}
path <- lapply(1:200, function(x){
  xy <- createPath()
  xy$TrackID <- x
  return(xy)
})
path <- do.call(rbind, path)

# Compute the distance of each coordinate to the village
path$Distance <- raster::extract(settl, cbind(path$x, path$y))

# Add some pseudo columns
path$FloodLevel <- "Min"
path$StepNumber <- 0
path$Area       <- 0

# Raster for visitation history
rx <- aggregate(r, fact = 2)
values(rx) <- 1:ncell(rx)

# Compute heatmaps and betweenness
betw <- betweenSims(path, rx, eps = 0.05)
heat <- rasterizeSims(path, r)

# Calculate point density within vicinity of humans
hwc <- subset(path, Distance < 15)
hwc <- rasterize(cbind(hwc$x, hwc$y), r, fun = "count")
hwc <- focal(hwc, w = matrix(rep(1, 9), nrow = 3), fun = "sum")

# Plot track
p1 <- ggplot() +
  geom_sf(data = st_as_sf(circle1), fill = "white", col = "white", alpha = 0.2) +
  geom_sf(data = st_as_sf(circle2), fill = "white", col = "white", alpha = 0.2) +
  geom_path(data = path, aes(x = x, y = y, group = as.factor(TrackID)), col = "orange") +
  theme_void()
p2 <- ggplot() +
  geom_raster(data = as.data.frame(betw, xy = T), aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(option = "magma") +
  geom_sf(data = st_as_sf(circle1), fill = "white", col = "white", alpha = 0.1) +
  geom_sf(data = st_as_sf(circle2), fill = "white", col = "white", alpha = 0.1) +
  theme_void() +
  theme(legend.position = "none")
p3 <- ggplot() +
  geom_raster(data = as.data.frame(heat, xy = T), aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(option = "magma", na.value = "black") +
  geom_sf(data = st_as_sf(circle1), fill = "white", col = "white", alpha = 0.2) +
  geom_sf(data = st_as_sf(circle2), fill = "white", col = "white", alpha = 0.2) +
  theme_void() +
  theme(legend.position = "none")
p4 <- ggplot() +
  geom_raster(data = as.data.frame(hwc, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_point(data = pts, aes(x = x, y = y), col = "red") +
  scale_fill_viridis_c(option = "plasma", na.value = "black") +
  geom_sf(data = st_as_sf(circle1), fill = "white", col = "white", alpha = 0.2) +
  geom_sf(data = st_as_sf(circle2), fill = "white", col = "white", alpha = 0.2) +
  theme_void() +
  theme(legend.position = "none")


# Store plots
ggsave(plot = p1, "05_Presentation/99_SchematicSimulations.png")
ggsave(plot = p2, "05_Presentation/99_SchematicHeatmap.png")
ggsave(plot = p3, "05_Presentation/99_SchematicBetweenness.png")
ggsave(plot = p4, "05_Presentation/99_SchematicHWC.png")
