################################################################################
#### Combining all Sources for Water Layers
################################################################################
# Description: In this script I combine the water layers from all different
# sources (Globeland, ORI, OSM, Dominik, MERIT)

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# load packages
library(tidyverse)  # For data wrangling
library(raster)     # To handle spatial data
library(terra)      # To handle spatial data
library(rgdal)      # To handle spatial data
library(pbmcapply)  # To use multiple cores

# Load the layers we want to merge (Globeland, Dynamic Floodmaps, Merit Rivers)
flood <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled" %>%
  dir(path = ., pattern = ".tif$", full.names = T) %>%
  stack()
globe <- raster("03_Data/02_CleanData/01_LandCover_WaterCover_GLOBELAND.tif")
coper <- raster("03_Data/02_CleanData/01_LandCover_WaterCover_COPERNICUS.tif")
merit <- raster("03_Data/02_CleanData/03_LandscapeFeatures_Rivers_MERIT.tif")

# We need to remove the globeland and copernicus waterbodies for the extent of
# the dynamic floodmaps. Let's get a polygon for the extent for which we have
# dynamic floodmaps first
p <- as(extent(trim(flood[[1]])), "SpatialPolygons")

# Replace the values below the polygon to 0
globe[cellsFromExtent(globe, p)] <- 0
coper[cellsFromExtent(coper, p)] <- 0

# Before we add the dynamic floodmaps, let's merge the globeland, copernicus,
# and merit data
globe <- max(globe, coper, merit)

# Fill the globeland layer with data from the dynamic floodmaps
globe <- suppressMessages(
  pbmclapply(
      X                   = 1: nlayers(flood)
    , mc.cores            = detectCores() - 1
    , ignore.interactive  = T
    , FUN                 = function(x){
      filled <- mask(globe, flood[[x]], maskvalue = 1, updatevalue = 1)
      filled <- writeRaster(filled, tempfile())
      return(filled)
  }) %>% stack()
)

# Let's also transfer the layernames
names(globe) <- names(flood)

# Let's also write the dates as seperate file
df <- data.frame(
    Layer = names(flood)
  , Date  = as.Date(substr(names(flood), start = 2, stop = 11), format = "%Y.%m.%d")
)
write_csv(df, "03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.csv")

# Save the result to file. Note that compressing the file will result in very
# slow extraction time. I'll thus leave the layers uncompressed, although this
# substantially inflates file sizes
terra::writeRaster(
    x         = rast(globe)
  , filename  = "03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd"
  , overwrite = TRUE
  , options   = c("COMPRESSION=NONE")
)
