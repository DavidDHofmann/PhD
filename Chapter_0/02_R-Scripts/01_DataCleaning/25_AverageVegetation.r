############################################################
#### Averaged Vegetation Layer
############################################################
# Description: Prepare layer that depicts the average vegetation layer. Note
# that this is basically the counterpart of the average water layer

# clear r's brain
rm(list = ls())

# Define the working directories
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(tidyverse)
library(raster)
library(terra)

# Make use of multiple cores
beginCluster()

# Load the vegetation layers again
files <- dir(
    path        = "03_Data/01_RawData/MODIS/MOD44B/Stitched"
  , pattern     = ".*MODIS.tif$"
  , full.names  = T
)
names <- substr(files, start = 55, stop = nchar(files) - 10)
modis <- stack(files)
names(modis) <- names

# Extract the separate layers
modis_shrub <- modis[[1]]
modis_noveg <- modis[[2]]
modis_trees <- modis[[3]]

# Values above 100 are water. Let's reclassify vegetation below water to 0
values(modis_shrub)[values(modis_shrub) > 100] <- 0
values(modis_noveg)[values(modis_noveg) > 100] <- 100
values(modis_trees)[values(modis_trees) > 100] <- 0

# Visualize
plot(stack(modis_shrub, modis_noveg, modis_trees))

# Load averaged water layer
water <- raster("03_Data/02_CleanData/01_LandCover_Water_Averaged.tif")

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
plot(stack(modis_shrub, modis_noveg, modis_trees))

# Make sure that the layers range from 0 to 1 rather than from 0 to 100
modis_shrub <- modis_shrub / 100
modis_noveg <- modis_noveg / 100
modis_trees <- modis_trees / 100

# Visualize layers again
plot(stack(modis_shrub, modis_noveg, modis_trees))

# Put the stacks into a list
modis <- list(modis_shrub, modis_noveg, modis_trees)

# Prepare filenames
names <- c(
    "03_Data/02_CleanData/01_LandCover_NonTreeVegetation_Averaged.tif"
  , "03_Data/02_CleanData/01_LandCover_NonVegetated_Averaged.tif"
  , "03_Data/02_CleanData/01_LandCover_TreeCover_Averaged.tif"
)

# Store the rasterstacks
for (i in 1:length(names)){
  writeRaster(modis[[i]], names[i], overwrite = TRUE)
}

# Terminate the cluster
endCluster()
