################################################################################
#### Graphical Abstract Plots
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(terra)       # For handling spatial data
library(raster)      # For handling spatial data
library(igraph)      # To handle networks
library(tidyverse)   # To wrangle data
library(sf)          # To plot spatial data
library(rgeos)       # For handling spatial data
library(Rcpp)        # Import C++ functions
library(pbmcapply)   # parallel mcapply with progress bar
library(spatstat)    # For handling spatial data
library(maptools)    # For handling spatial data

# Load custom functions
source("02_R-Scripts/00_Functions.R")
sourceCpp("02_R-Scripts/00_Functions.cpp")

# We'll put the plots into a subfolder, as we're going to compile them later in
# a libre impress file
dir.create("04_Manuscript/Figures/GraphicalAbstract", showWarnings = F)

################################################################################
#### Floodmap Single Extremes
################################################################################
# Find all floodmaps and their associated date
flood <- "03_Data/01_RawData/FLOODMAPS" %>%
  dir(pattern = ".tif$", full.names = T) %>%
  rast()

# Let's determine the flood and cloud extent in each image
flood_summary <- flood %>%
  expanse(unit = "km", byValue = T) %>%
  as.data.frame() %>%
  pivot_wider(
    , id_cols     = layer
    , names_from  = value
    , values_from = area
    , values_fill = 0
  ) %>%
  rename(Flood = "0", Dryland = "255", Cloud = "127") %>%
  mutate(Total = Flood + Dryland + Cloud) %>%
  mutate(Cloud = Cloud / Total)

# Let's remove layers where there is loads of clouds and find the maps with
# minimal and maximum flood
toplot <- flood_summary %>%
  subset(Cloud < 0.25) %>%
  subset(Flood %in% c(min(Flood), max(Flood))) %>%
  pull(layer)

# Plot them both and store the plots to file
extremes <- flood[[toplot]]
extremes <- classify(extremes, cbind(c(0, 127, 255), c(1, 0, 0)))
names(extremes) <- c("max", "min")
extremes <- as.data.frame(extremes, xy = T)

# Plot
p1 <- ggplot(extremes, aes(x = x, y = y, fill = as.factor(min))) +
  geom_raster() +
  scale_fill_manual(values = c("white", "cornflowerblue")) +
  theme_void() +
  theme(legend.position = "none")
p2 <- ggplot(extremes, aes(x = x, y = y, fill = as.factor(max))) +
  geom_raster() +
  scale_fill_manual(values = c("white", "cornflowerblue")) +
  theme_void() +
  theme(legend.position = "none")

# Store plots
ggsave(plot = p1, "04_Manuscript/Figures/GraphicalAbstract/FloodExtremes_Maximum.png")
ggsave(plot = p2, "04_Manuscript/Figures/GraphicalAbstract/FloodExtremes_Minimum.png")

################################################################################
#### Floodmap Composite Extremes
################################################################################
# Load composite extremes
water <- rast("03_Data/02_CleanData/WaterCover.tif")[[c("min", "max")]]
water <- crop(water, flood)
water <- as.data.frame(water, xy = T)

# Plot
p1 <- ggplot(water, aes(x = x, y = y, fill = as.factor(min))) +
  geom_raster() +
  scale_fill_manual(values = c("white", "cornflowerblue")) +
  theme_void() +
  theme(legend.position = "none")
p2 <- ggplot(water, aes(x = x, y = y, fill = as.factor(max))) +
  geom_raster() +
  scale_fill_manual(values = c("white", "cornflowerblue")) +
  theme_void() +
  theme(legend.position = "none")

# Store plots
ggsave(plot = p1, "04_Manuscript/Figures/GraphicalAbstract/FloodExtremesComposite_Maximum.png")
ggsave(plot = p2, "04_Manuscript/Figures/GraphicalAbstract/FloodExtremesComposite_Minimum.png")

################################################################################
#### Static Covariates
################################################################################
# Load static covariates
trees <- rast("03_Data/02_CleanData/TreeCover.tif")[[2]]
shrub <- rast("03_Data/02_CleanData/ShrubCover.tif")[[2]]
human <- rast("03_Data/02_CleanData/HumanInfluence.tif")

# Crop them to the extent of the delta
trees <- crop(trees, flood) %>% as.data.frame(xy = T)
shrub <- crop(shrub, flood) %>% as.data.frame(xy = T)
human <- crop(human, flood) %>% as.data.frame(xy = T)

# Plot
p1 <- ggplot(trees, aes(x = x, y = y, fill = mean)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "sqrt") +
  theme_void() +
  theme(legend.position = "none")
p2 <- ggplot(shrub, aes(x = x, y = y, fill = mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  theme_void() +
  theme(legend.position = "none")
p3 <- ggplot(human, aes(x = x, y = y, fill = HumanInfluence)) +
  geom_raster() +
  scale_fill_viridis_c() +
  theme_void() +
  theme(legend.position = "none")

# Store plots
ggsave(plot = p1, "04_Manuscript/Figures/GraphicalAbstract/CovariatesTrees.png")
ggsave(plot = p2, "04_Manuscript/Figures/GraphicalAbstract/CovariatesShrubs.png")
ggsave(plot = p3, "04_Manuscript/Figures/GraphicalAbstract/CovariatesHumans.png")

################################################################################
#### Schematic Visualization of the Different Metrics
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
pts <- rbind(c(45, 65), c(75, 35))
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
rx <- disaggregate(r, fact = 2)
# rx <- aggregate(r, fact = 2)
values(rx) <- 1:ncell(rx)

# Compute heatmaps and betweenness
betw <- betweenSims(path, rx, eps = 0.05)
heat <- rasterizeSims(path, r)

# Calculate point density within vicinity of humans
hwc <- subset(path, Distance < 15)
hwc <- rasterize(cbind(hwc$x, hwc$y), r, fun = "count", background = 0)

# Do some styling
hwc_styled  <- trim(focal(sqrt(hwc), w = focalWeight(hwc, d = 1, type = "Gauss"), fun = "mean"))
heat_styled <- trim(focal(sqrt(heat), w = focalWeight(heat, d = 2, type = "Gauss"), fun = "mean"))
betw_styled <- trim(focal(sqrt(betw), w = focalWeight(betw, d = 2, type = "Gauss"), fun = "mean"))

# Visualize all
p1 <- ggplot(as.data.frame(heat_styled, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colors = hcl.colors(n = 100, palette = "Spectral", rev = T)) +
  theme_void() +
  theme(legend.position = "none")
p2 <- ggplot(as.data.frame(betw_styled, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colors = viridis::magma(100)) +
  theme_void() +
  theme(legend.position = "none")
p3 <- ggplot(as.data.frame(hwc_styled, xy = T), aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_gradientn(colors = terrain.colors(100, rev = T)) +
  theme_void() +
  theme(legend.position = "none")

# Store to file
ggsave(plot = p1, "04_Manuscript/Figures/GraphicalAbstract/MetricsHeatmap.png")
ggsave(plot = p2, "04_Manuscript/Figures/GraphicalAbstract/MetricsBetweenness.png")
ggsave(plot = p3, "04_Manuscript/Figures/GraphicalAbstract/MetricsHWC.png")
