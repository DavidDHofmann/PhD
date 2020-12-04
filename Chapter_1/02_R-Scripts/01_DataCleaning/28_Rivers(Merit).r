################################################################################
#### Preparation of the River Widths By Merit
################################################################################
# Description: Preparation and cleaning of the Merit river layers

# Clear R's brain
rm(list = ls())

# Specify the working directories
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load packages
library(raster)
library(rgdal)

# Identify the files to load
files <- dir(
    path        = "03_Data/01_RawData/MERIT"
  , pattern     = ".tif$"
  , full.names  = T
)

# Load them into a list
rivers <- lapply(files, raster)

# Merge them into a single file
rivers <- do.call(merge, rivers)

# Reclassify pixel values and keep only those rivers with a desired width
rcl <- data.frame(from = c(-Inf, 10), to = c(10, Inf), new = c(0, 1))
rivers <- reclassify(rivers, rcl)

# Load the reference raster in order to reproject the river layer
r250 <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# First, we aggregate the layer to match the resolution of the reference raster
# The reference raster is resolve around 250 meters, the river layer at around
# 90 meters
rivers <- aggregate(rivers, fact = round(250 / 90), fun = max)

# Secondly, we ca resample the river layer to match the origin and extent of the
# reference raster
rivers <- resample(rivers, r250, "ngb")

# Store the final result to file
writeRaster(
    rivers
  , "03_Data/02_CleanData/03_LandscapeFeatures_Rivers_Merit.tif"
  , overwrite = TRUE
)
