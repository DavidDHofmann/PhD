############################################################
#### Resampling and Cropping of the Croplands Agriculture Data
############################################################
# Description: In this script I prepare the layer downloaded from croplands. It
# depicts areas in which there are agricultural fields. Unfortunately the layer
# does not show all fields, which is why we have to combine it with other layers
# later.

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)

# Make use of multicore abilities
beginCluster()

# Import the Croplands dataset
crops <- raster("03_Data/01_RawData/CROPLANDS/Croplands.tif")

# Load the reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# # Crop the croplands layer to the desired extent
# crops <- crop(crops, r)
#
# # Because the original file is so massive I will replace it with the cropped
# # version
# writeRaster(
#     x         = crops
#   , filename  = "03_Data/01_RawData/CROPLANDS/Croplands.tif"
#   , overwrite = T
# )

# Aggregate the croplands dataset to match the resolution of the reference
# raster
crops <- aggregate(crops, fact = round(250 / 30), fun = max)

# Now water is still included in the raster. Let's reclassify the values so we
# only keep the crops
rcl <- data.frame(old = c(1,2), new = c(0,1))
crops <- reclassify(crops, rcl)

# Resample the layer to match the reference raster
crops_res <- resample(crops, r, "ngb")

# Save the result to file
writeRaster(crops_res
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_Croplands.tif"
  , overwrite = TRUE
)

# End cluster
endCluster()
