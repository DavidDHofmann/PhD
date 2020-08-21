############################################################
#### Combining all Sources for Water Layers
############################################################
# Description: In this script I combine the water layers from all different
# sources (Globeland, ORI, OSM, Dominik, MERIT)

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_0"
setwd(wd)

# load packages
library(tidyverse)
library(raster)
library(rgeos)
library(gdalUtils)
library(rgdal)
library(terra)
library(parallel)

# Allow for multicore use
beginCluster()

# Load the layers we want to merge (Globeland, Dynamic Floodmaps, Merit Rivers)
flood <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled" %>%
  dir(path = ., pattern = ".tif$", full.names = T) %>%
  stack()
globe <- raster("03_Data/02_CleanData/01_LandCoverClasses30_Globeland.tif")
merit <- raster("03_Data/02_CleanData/03_LandscapeFeatures_Rivers_Merit.tif")

# We need to remove the globeland waterbodies (Code 40) for the extent of the
# dynamic floodmaps. Let's get a polygon for the extent for which we have
# dynamic floodmaps first
p <- as(extent(flood[[1]]), "SpatialPolygons")

# Replace the values below the polygon to 0
globe[cellsFromExtent(globe, p)] <- 0

# Before we add the dynamic floodmaps, let's merge the globeland and merit data
globe <- max(globe, merit)

# We need to extend the floodmaps to match the extent of the globeland layer
flood <- suppressMessages(
  mclapply(1:nlayers(flood), mc.cores = detectCores() - 1, function(x){
    flood_extended <- extend(flood[[x]], globe, value = NA)
    flood_extended <- writeRaster(flood_extended, tempfile())
    return(flood_extended)
  }) %>% stack()
)

# Fill the globeland layer with floodmaps
globe <- suppressMessages(
  mclapply(1:nlayers(flood), mc.cores = detectCores() - 1, function(x){
    filled <- mask(globe, flood[[x]], maskvalue = 1, updatevalue = 1)
    filled <- writeRaster(filled, tempfile())
    return(filled)
  }) %>% stack()
)

# Convert to terra raster
globe <- rast(globe)

# Let's also transfer the layernames
names(globe) <- names(flood)

# Save the result to file
terra::writeRaster(
    globe
  , filename  = "03_Data/02_CleanData/01_LandCover_Water(Merged).tif"
  , overwrite = TRUE
)

# Terminate the cluster
endCluster()
