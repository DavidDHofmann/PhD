################################################################################
#### Rasterization of Farms by Gabriele
################################################################################
# Description: Here I rasterize all of the farms from the shapefile that
# Gabriele provided

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)     # To handle raster data
library(rgdal)      # To handle vector data
library(gdalUtils)  # Some helpfull tools

# Make use of multicore abilities
beginCluster()

# Import Gabriele's farms
crops <- readOGR("03_Data/01_RawData/GABRIELE/Farmland.shp")

# Write the layer to the cleaned data
writeOGR(
    crops
  , dsn       = "03_Data/02_CleanData"
  , layer     = "04_AnthropogenicFeatures_Farms_Gabriele"
  , driver    = "ESRI Shapefile"
  , overwrite = T
)

# Load reference raster for 250 meters. Fill its cells with only 0s for the
# rasterization
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")
values(r) <- 0

# I first tried to run the rasterization with the rasterize command. However, I
# found that gdal's gdal_rasterize works much faster. To use it, we need to
# save a new rasterfile.
writeRaster(
    r
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_Farms_Gabriele.tif"
  , overwrite = TRUE
)

# Rasterize the farms
gdal_rasterize(
    src_datasource  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Farms_Gabriele.shp"
  , dst_filename    = "03_Data/02_CleanData/04_AnthropogenicFeatures_Farms_Gabriele.tif"
  , burn            = 1
  , at              = TRUE
)

# End cluster
endCluster()
