################################################################################
#### Rasterization of Farms by Gabriele
################################################################################
# Description: Here I rasterize all of the farms from the shapefile that
# Gabriele provided

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)     # To handle raster data
library(terra)      # To handle raster data
library(rgdal)      # To handle vector data
library(gdalUtils)  # Some helpfull tools

# Import Gabriele's farms
crops <- readOGR("03_Data/01_RawData/GABRIELE/Farmland.shp")

# Write the layer to the cleaned data
writeOGR(
    crops
  , dsn       = "03_Data/02_CleanData"
  , layer     = "04_AnthropogenicFeatures_Farms_GABRIELE"
  , driver    = "ESRI Shapefile"
  , overwrite = T
)

# Load reference raster for 250 meters
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")

# Rasterize the farms
r <- rasterize(vect(crops), r, field = 1, background = 0)
writeRaster(
    x         = raster(r)
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_GABRIELE.tif"
  , overwrite = TRUE
)
