################################################################################
#### Averaged Vegetation Layer
################################################################################
# Description: Prepare layer that depicts the average vegetation layer. Note
# that this is basically the counterpart of the average water layer

# clear r's brain
rm(list = ls())

# Define the working directories
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)  # For data wrangling
library(raster)     # For manipulating spatial data
library(terra)      # For manipulating spatial data

# Load the vegetation layers again
files <- dir(
    path        = "03_Data/01_RawData/MODIS/MOD44B/Stitched"
  , pattern     = ".*MODIS.tif$"
  , full.names  = T
)
names <- substr(
    x     = basename(files)
  , start = 14
  , stop  = nchar(basename(files)) - 10
)
modis <- rast(files)
names(modis) <- names

# Extract the separate layers
modis_shrub <- modis[[1]]
modis_noveg <- modis[[2]]
modis_trees <- modis[[3]]

# Values above 100 are water. Let's reclassify vegetation below water to 0
values(modis_shrub)[values(modis_shrub) > 100] <- 0
values(modis_noveg)[values(modis_noveg) > 100] <- 100
values(modis_trees)[values(modis_trees) > 100] <- 0

# Visualize again
plot(c(modis_shrub, modis_noveg, modis_trees))

# Load averaged water layer
water <- rast("03_Data/02_CleanData/01_LandCover_WaterCoverAveraged_MERGED.grd")

# Replace values below water to 0
modis_shrub <- mask(modis_shrub
  , mask        = water
  , maskvalue   = 1
  , updatevalue = 0
)
modis_noveg <- mask(modis_noveg
  , mask        = water
  , maskvalue   = 1
  , updatevalue = 100
)
modis_trees <- mask(modis_trees
  , mask        = water
  , maskvalue   = 1
  , updatevalue = 0
)

# Visualize layers again
plot(c(modis_shrub, modis_noveg, modis_trees))

# Make sure that the layers range from 0 to 1 rather than from 0 to 100
modis_shrub <- modis_shrub / 100
modis_noveg <- modis_noveg / 100
modis_trees <- modis_trees / 100

# Visualize layers again
plot(c(modis_shrub, modis_noveg, modis_trees))

# Put the stacks into a list
modis <- list(modis_shrub, modis_noveg, modis_trees)

# Prepare filenames
names <- c(
    "03_Data/02_CleanData/01_LandCover_NonTreeVegetationAveraged_MODIS.tif"
  , "03_Data/02_CleanData/01_LandCover_NonVegetatedAveraged_MODIS.tif"
  , "03_Data/02_CleanData/01_LandCover_TreeCoverAveraged_MODIS.tif"
)

# Store the rasterstacks
for (i in 1:length(names)){
  writeRaster(raster(modis[[i]]), names[i], overwrite = TRUE)
}
