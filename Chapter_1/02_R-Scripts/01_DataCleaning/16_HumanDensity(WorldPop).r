################################################################################
#### Cleaning Worldpop Human Density Data
################################################################################
# Description: In this script I'll clean the population density estimates
# downloaded from here: https://www.worldpop.org/

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
    path       = "03_Data/01_RawData/WORLDPOP"
  , full.names = T
  , pattern    = ".tif$"
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
  , filename  = "03_Data/01_RawData/WORLDPOP/HumanDensity.tif"
  , overwrite = TRUE
)

################################################################################
#### Aggregate and Resample
################################################################################
# Aggregate the layer to 250m
coarse <- aggregate(dat, fact = round(250 / 90), fun = sum)

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
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanDensity_WORLDPOP.tif"
  , overwrite = TRUE
)

# End cluster
endCluster()
