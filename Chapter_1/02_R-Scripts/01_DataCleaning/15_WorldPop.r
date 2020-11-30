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

# Load packages
library(terra)
library(raster)
library(sp)
library(tidyverse)
library(viridis)

# Load the data
dat <- dir(
    path       = "03_Data/01_RawData/WORLDPOP"
  , full.names = T
  , pattern    = ".tif$"
)
dat <- lapply(dat, raster)

# Merge the tiles
dat <- mosaic(dat[[1]], dat[[2]], dat[[3]], dat[[4]], dat[[5]], fun = max)

# Convert to terra
dat <- rast(dat)

# Crop the data to our study area
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")
dat <- crop(dat, r, snap = "out")

# Replace NAs with 0s
values(dat)[is.na(values(dat))] <- 0

# Store the raster to file
writeRaster(dat
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanDensity_WORLDPOP.tif"
  , overwrite = TRUE
)
