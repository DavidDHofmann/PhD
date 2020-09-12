################################################################################
#### Plot of Historic Range
################################################################################
# Description: A cartoon like plot of the historic range of african wild dogs

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)  # For data wrangling
library(raster)     # To handle spatial data
library(rgeos)      # To manipulate spatial data
library(rgdal)      # To read and write spatial data
library(tmap)       # To plot nice maps
library(smoothr)    # To smooth spatial objects
library(davidoff)   # Custom functions

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

################################################################################
#### Plot
################################################################################
# Prepare a map of Africa
p1 <- tm_shape(africa2) +
    tm_polygons(col = "gray70", border.col = "gray70", lwd = 2) +
  tm_shape(africa) +
    tm_polygons(
        col = "black"
      , lwd = 0.7
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
  tm_shape(strong) +
    tm_polygons(
        col           = lighten("orange", 1.3)
      , border.alpha  = 0
    ) +
  tm_layout(
      asp         = 0.8
    , frame       = "white"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.color    = "white"
)
