################################################################################
#### Visualization of Source Areas
################################################################################
# Description: Simulation of dispersal

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_10")

# Load required packages
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(tidyverse)      # To wrangle data
library(mapview)        # For interactive plots
library(sf)             # For plotting

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Load Required Data
################################################################################
# Load the required spatial layers layers
layers <- rast("03_Data/02_CleanData/SpatialLayers.tif")

# Define an extent and create a reference raster
ext <- ext(20.5, 26.5, -21.5, -17.5)
r <- rast(ext, resolution = 100 / 111000, crs = "epsg:4326")

# Load areas that we want to remove from the domain of "suitable" habitats
humans <- layers[["HumansBuff5000"]]
humans <- crop(humans, ext)
humans <- humans > 0
water  <- layers[["Water"]]
water  <- crop(water, ext)
major  <- vect("03_Data/02_CleanData/MajorWaters.shp")
water  <- mask(water, major, updatevalue = 1, inverse = T)

# Let's create a mask that we'll use to remove everything that is in one of the
# two layers
m <- water + humans
m <- m > 0
m <- as.polygons(m)
m <- m[m$Water == 1]
m <- buffer(m, width = 0)

# Specify simulation parameters
min_area <- 50
max_area <- 900
n        <- 250
buff     <- 10000

# Get an extent from the input raster, and buffer if desired
area <- ext(r)
area <- as.polygons(area, crs = crs(r))
area <- buffer(area, width = buff)
area <- area - m

# Randomly place points in areas that shall not be removed
pts <- spatSample(area, size = n)
crs(pts) <- CRS("+init=epsg:4326")

# Run triangulation to get polygons
pols    <- voronoi(pts)
pols    <- pols - m
pols    <- crop(pols, ext(r))
pols    <- disagg(pols)
crs(pols) <- CRS("+init=epsg:4326")

# Identify polygons that are too big
pols$Area <- expanse(pols, unit = "km")
pols$Big  <- pols$Area > max_area

# Split polygons that are too big and remove polygons that are too small
pols      <- splitPolyUntil(pols, max_area = max_area)
pols      <- pols[pols$Area >= min_area, ]
pols      <- disagg(pols)
pols$ID   <- sample(1:length(pols))

################################################################################
#### Prepare Plots
################################################################################
# Convert data to plot with sf
m <- st_as_sf(m)
r <- as.data.frame(r, xy = T)
p <- st_as_sf(pts)
pols <- st_as_sf(pols)

# Prepare plots
p1 <- ggplot() +
  geom_sf(data = m, fill = "gray10", color = NA) +
  coord_sf(
      crs  = 4326
    , xlim = c(xmin(ext), xmax(ext))
    , ylim = c(ymin(ext), ymax(ext))
  ) +
  theme_minimal()

p2 <- ggplot() +
  geom_sf(data = m, fill = "gray10", color = NA) +
  geom_sf(data = p, color = "orange", size = 0.2) +
  coord_sf(
      crs  = 4326
    , xlim = c(xmin(ext), xmax(ext))
    , ylim = c(ymin(ext), ymax(ext))
  ) +
  theme_minimal()

p3 <- ggplot() +
  geom_sf(data = m, fill = "gray10", color = NA) +
  geom_sf(data = pols, aes(fill = as.factor(ID)), color = NA) +
  geom_sf(data = p, color = "orange", size = 0.2) +
  coord_sf(
      crs  = 4326
    , xlim = c(xmin(ext), xmax(ext))
    , ylim = c(ymin(ext), ymax(ext))
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  ) +
scale_fill_viridis_d()

# Store the plots
ggsave("05_Presentation/HRSimulation_01.png"
  , plot   = p1
  , bg     = "transparent"
  , width  = 8
  , height = 4
  , scale  = 1.1
)
ggsave("05_Presentation/HRSimulation_02.png"
  , plot   = p2
  , bg     = "transparent"
  , width  = 8
  , height = 4
  , scale  = 1.1
)
ggsave("05_Presentation/HRSimulation_03.png"
  , plot   = p3
  , bg     = "transparent"
  , width  = 8
  , height = 4
  , scale  = 1.1
)
