################################################################################
#### Cleaning Facebook Human Density Data
################################################################################
# Description: Preparation of the tiles that were downloaded from facebook
# https://data.humdata.org/dataset/highresolutionpopulationdensitymaps
# Includes stitching, aggregating and reprojecting the tiles

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)     # To handle raster data
library(rgdal)      # To handle spatial data

# Make use of multicore abilities
beginCluster()

################################################################################
#### Stitch Tiles
################################################################################
# Load the data
dat <- dir(
    path       = "03_Data/01_RawData/FACEBOOK"
  , pattern    = "population.*tif$"
  , full.names = T
)
dat <- lapply(dat, raster)

# Merge all tiles together
dat <- do.call(merge, dat)

# Crop merged files to our extent
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")
dat <- crop(dat, r, snap = "out")

# Store the raster to file
writeRaster(
    x         = dat
  , filename  = "03_Data/01_RawData/FACEBOOK/HumanDensity.tif"
  , overwrite = TRUE
)

################################################################################
#### Aggregate & Crop
################################################################################
# Aggregate the layer to 250m
coarse <- aggregate(dat, fact = round(250 / 30), fun = sum)

# Finally we need to resample the layer. Lets load the reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Resample the population density layer to the reference raster
coarse <- resample(coarse, r, "bilinear")

# Replace NAs and values below 0 with 0s
coarse <- reclassify(coarse, rcl = c(NA, NA, 0))
coarse <- reclassify(coarse, rcl = c(-Inf, 0, 0))

# Look at final raster
coarse

# Store the result
writeRaster(
    x         = coarse
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanDensity_FACEBOOK.tif"
  , overwrite = TRUE
)

# End cluster
endCluster()
