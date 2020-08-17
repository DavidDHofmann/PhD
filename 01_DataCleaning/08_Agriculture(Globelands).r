############################################################
#### Resampling and Cropping of the Globeland Agriculture Data
############################################################
# Description: Although we will remove agricultural fields from the final land
# cover class layer we can still use globeland's agricultural fields to augment
# our croplands layer. In this script I therefore extract agricultural fields so
# we can merge them.

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(raster)

# Make use of multicore abilities
beginCluster()

# Import the Globelands dataset
crops <- raster("03_Data/01_RawData/GLOBELAND/Globeland.tif")

# Reclassify the layer so that only cultivated land (Code 10) is kept as 1, the
# rest as 0
rcl <- data.frame(
    old = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 255)
  , new = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
crops <- reclassify(crops, rcl)

# Load the reference raster
r250 <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Aggregate the globelands dataset to match the resolution of the reference
# raster
crops <- aggregate(crops, fact = round(250 / 30), fun = max)

# Resample the layer to match the reference raster
crops_res <- resample(crops, r250, "ngb")

# Save the result to file
writeRaster(crops_res
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_Globelands.tif"
  , overwrite = TRUE
)

# End cluster
endCluster()
