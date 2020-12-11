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
library(terra)      # To handle raster data

################################################################################
#### Stitch Tiles
################################################################################
# Load the data
dat <- dir(
    path    = "03_Data/01_RawData/FACEBOOK"
  , pattern = "population.*tif$"
  , full.names = T
)
dat <- lapply(dat, raster)

# Merge the tiles
dat <- mosaic(dat[[1]], dat[[2]], dat[[3]], dat[[4]], dat[[5]], fun = max)

# Convert to terra
dat <- rast(dat)

# Crop the merged tiles to our extent (with a slight buffer of 1km)
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")
extent <- vect(as(extent(raster(r)) + c(-1, 1, -1, 1) / 111, "SpatialPolygons"))
dat <- crop(dat, extent, snap = "out")

# Replace NAs with 0s
dat <- classify(dat, rcl = matrix(c(NA, 0), nrow = 1))

# Store the raster to file
writeRaster(
    x         = raster(dat)
  , filename  = "03_Data/01_RawData/FACEBOOK/HumanDensity.tif"
  , overwrite = TRUE
)

################################################################################
#### Aggregate and Resample
################################################################################
# Aggregate the layer to 250m
coarse <- aggregate(dat, fact = round(250 / 30), fun = sum)

# Resample the population density layer to the reference raster
coarse <- resample(coarse, r, "bilinear")

# Store the result
writeRaster(
    x         = raster(coarse)
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanDensity_FACEBOOK.tif"
  , overwrite = TRUE
)
