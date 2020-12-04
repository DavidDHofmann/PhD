################################################################################
#### Resampling and Cropping of the Globeland Agriculture Data
################################################################################
# Description: We can use globeland's agricultural fields to augment our
# croplands layer. In this script I therefore extract agricultural fields so we
# can merge them.

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)   # To handle raster data
library(terra)    # To handle raster data

# Import the Globelands dataset
crops <- rast("03_Data/01_RawData/GLOBELAND/Globeland.tif")

# Keep only crops
crops <- crops == 10

# Load the reference raster
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")

# Aggregate the globelands dataset to match the resolution of the reference
# raster
crops <- aggregate(crops, fact = round(250 / 30), fun = max)

# Resample the layer to match the reference raster
crops_res <- resample(crops, r, "near")

# Save the result to file
writeRaster(
    x         = raster(crops_res)
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_GLOBELAND.tif"
  , overwrite = TRUE
)
