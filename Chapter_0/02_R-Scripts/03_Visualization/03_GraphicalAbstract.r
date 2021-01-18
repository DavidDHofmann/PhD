################################################################################
#### Plot of Historic Range
################################################################################
# Description: A cartoon like plot of the historic range of african wild dogs

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(tidyverse)  # For data wrangling
library(raster)     # To handle spatial data
library(rgeos)      # To manipulate spatial data
library(rgdal)      # To read and write spatial data
library(tmap)       # To plot nice maps
library(smoothr)    # To smooth spatial objects
library(davidoff)   # Custom functions
library(viridis)    # For nice colors
library(rasterVis)  # For nice spatial plots

################################################################################
#### Data Preparation
################################################################################
# Load map of africa
africa    <- readOGR("03_Data/02_CleanData/00_General_Africa.shp")
africa2   <- readOGR("03_Data/02_CleanData/00_General_Africa.shp")
historic  <- readOGR("03_Data/01_RawData/DAVID/HistoricRange.shp")
dogs      <- readOGR("03_Data/02_CleanData/00_General_WildDogs_IUCN.shp")
kaza      <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")

# Buffer africa slightly
africa <- gBuffer(africa, width = 0.01)

# Disaggregate
africa <- disaggregate(africa)

# Identify area of each polygon
africa$Size <- gArea(africa, byid = T)

# Keep only the largest two
africa <- subset(africa, Size %in% sort(africa$Size, decreasing = T)[1:2])

# Remove madagascar
africa <- africa[1, ]

# Simplify and smoothen africa shape
africa <- gSimplify(africa, tol = 0.5)
africa <- smooth(africa, method = "ksmooth")

# Smoothen and crop historic range
historic <- smooth(historic, method = "ksmooth")
historic <- gIntersection(historic, africa, byid = F)

# Simplify and smoothen dog distribution
dogs <- gSimplify(dogs, tol = 0.2)
dogs <- smooth(dogs, method = "ksmooth")

# Create a buffered polygon of africa
africa2 <- gBuffer(africa, width = 100/111000)

# Simplify and smoothen kaza
kaza <- gSimplify(kaza, tol = 0.1)
kaza <- smooth(kaza, method = "ksmooth")

# Identify wild dog strongholds
strong <- rbind(disaggregate(dogs)[c(7, 19, 29), ])

# Prepare a map of Africa (with KAZA)
p1 <- tm_shape(africa2) +
    tm_polygons(col = "gray70", border.col = "gray70", lwd = 2) +
  tm_shape(africa) +
    tm_polygons(
        col = "black"
      , lwd = 10
      , border.col = "gray70"
    ) +
  tm_shape(historic) +
    tm_polygons(
        col = "orange"
      , lwd = 0.7
      , border.col = "black"
      , alpha = 0.2
    ) +
  tm_shape(dogs) +
    tm_polygons(
        col           = "orange"
      , alpha         = 0.8
      , border.alpha  = 0
    ) +
  tm_shape(kaza) +
    tm_borders(
        col = "white"
      , lwd = 15
    ) +
  tm_layout(
      asp         = 0.8
    , frame       = "white"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.color    = "transparent"
)

png("test.png", width = 1980, height = 2200, bg = "transparent")
p1
dev.off()

############################################################
#### Permeability Surface
############################################################
# Load required data
permeability  <- "03_Data/03_Results/99_PermeabilityMap.tif" %>%
  raster()
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>%
  readOGR() %>%
  as("SpatialLines")

# Plot
p <- tm_shape(permeability) +
    tm_raster(palette = "viridis", style = "cont", legend.show = F) +
  tm_shape(kaza) +
    tm_lines(col = "white", lwd = 15)

# Store
png("test.png", width = 1980, height = 2200, bg = "transparent")
p
dev.off()

############################################################
#### Path Map
############################################################
# Load required data
paths  <- "03_Data/03_Results/99_LeastCostPaths.tif" %>%
  raster()
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>%
  readOGR() %>%
  as("SpatialLines")

# Plot
p <- tm_shape(sqrt(paths)) +
    tm_raster(palette = "-Spectral", style = "cont", legend.show = F) +
  tm_shape(kaza) +
    tm_lines(col = "white", lwd = 15)

# Store
png("test2.png", width = 1980, height = 2200, bg = "transparent")
p
dev.off()

############################################################
#### Corridor Map
############################################################
# Load required data
corrs  <- "03_Data/03_Results/99_LeastCostCorridors2.tif" %>%
  raster()
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>%
  readOGR() %>%
  as("SpatialLines")

# Plot
p <- tm_shape(corrs) +
    tm_raster(palette = "-Spectral", style = "cont", legend.show = F) +
  tm_shape(kaza) +
    tm_lines(col = "white", lwd = 15)

# Store
png("test3.png", width = 1980, height = 2200, bg = "transparent")
p
dev.off()
