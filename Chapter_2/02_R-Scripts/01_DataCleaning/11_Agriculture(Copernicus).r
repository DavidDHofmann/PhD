################################################################################
#### Resampling and Cropping of the Copernicus Agriculture Data
################################################################################
# Description: We can use copernicus's agricultural fields to augment our
# croplands layer. In this script I therefore extract agricultural fields so we
# can merge them.

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)   # To handle raster data

# Make use of multiple cores
beginCluster()

# Import the Copernicus dataset
crops <- raster("03_Data/01_RawData/COPERNICUS/Copernicus.tif")

# Load the reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Keep only crops
crops <- crops == 40

# Aggregate the copernicus dataset to match the resolution of the reference
# raster
fact <- res(r)[1] / res(crops)[1]
crops <- aggregate(crops, fact = round(fact), fun = max)

# Resample the layer to match the reference raster
crops_res <- resample(crops, r, "ngb")

# Check NAs
sum(is.na(values(crops_res)))

# Save the result to file
writeRaster(
    x         = crops_res
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_COPERNICUS.tif"
  , overwrite = TRUE
)

# End cluster
endCluster()
