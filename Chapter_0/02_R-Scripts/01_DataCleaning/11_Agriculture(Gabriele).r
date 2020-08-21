############################################################
#### Rasterization of Farms by Gabriele
############################################################
# Description: Here I rasterize all of the farms from the shapefile that
# Gabriele provided

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(raster)
library(gdalUtils)
library(rgdal)

# Make use of multicore abilities
beginCluster()

# Import Gabriele's farms
crops <- shapefile("03_Data/01_RawData/GABRIELE/Farmland")

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
r250 <- raster("03_Data/02_CleanData/00_General_Raster250.tif")
values(r250) <- 0

# I first tried to run the rasterization with the rasterize command. However, I
# found that gdal's gdal_rasterize works much faster. To use it, we need to
# save a new rasterfile.
writeRaster(
    r250
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_Farms_Gabriele.tif"
  , overwrite = TRUE
)

# Rasterize the farms
gdal_rasterize(
    src_datasource = "03_Data/02_CleanData/04_AnthropogenicFeatures_Farms_Gabriele.shp"
  , dst_filename = "03_Data/02_CleanData/04_AnthropogenicFeatures_Farms_Gabriele.tif"
  , burn = 1
  , at = TRUE
)

# End cluster
endCluster()
