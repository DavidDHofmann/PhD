############################################################
#### Preparation of the Wild Dog Distribution Data
############################################################
# Description: Preparation of a shapefile of the distribution of wild dogs.

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(rgdal)
library(raster)

# Load the shapefile
dogs <- readOGR("03_Data/01_RawData/IUCN/data_0.shp")

# Visualize
plot(dogs, col = "purple")

# Store the file
writeOGR(
    dogs
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_WildDogs_IUCN"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
