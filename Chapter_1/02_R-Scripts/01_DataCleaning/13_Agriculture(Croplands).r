################################################################################
#### Resampling and Cropping of the Croplands Agriculture Data
################################################################################
# Description: In this script I prepare the layer downloaded from croplands. It
# depicts areas in which there are agricultural fields. Unfortunately the layer
# does not show all fields, which is why we have to combine it with other layers
# later.

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)  # To handle reaster data
library(terra)   # To handle reaster data

# Import the Croplands dataset
crops <- rast("03_Data/01_RawData/CROPLANDS/Croplands.tif")

# Load the reference raster
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")

# Crop to reference raster
crops <- crop(crops, r)

# Aggregate the croplands dataset to match the resolution of the reference
# raster
crops <- aggregate(crops, fact = round(250 / 30), fun = max)

# Now water is still included in the raster. Let's reclassify the values so we
# only keep the crops
rcl <- data.frame(old = c(1, 2), new = c(0, 1))
crops <- classify(crops, rcl)

# Resample the layer to match the reference raster
crops_res <- resample(crops, r, "near")

# Save the result to file
writeRaster(
    x         = raster(crops_res)
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_CROPLANDS.tif"
  , overwrite = TRUE
)
