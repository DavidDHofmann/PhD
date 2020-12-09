################################################################################
#### Distance to Water
################################################################################
# Description: Prepare layer that depicts the distance to the closest water cell

# clear r's brain
rm(list = ls())

# Define the working directories
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(spatstat)     # To calculate distances quickly
library(maptools)     # To calculate distances quickly
library(tidyverse)    # For data wrangling
library(raster)       # To handle spatial data
library(terra)        # To handle spatial data
library(davidoff)     # Custom functions

# Load water layers
water <- rast("03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd")

# We want to create an "average" water map. For this we calculate how often
# each cell was covered by water. If this is more than a desired threshold, we
# will use the cell for our "averaged" map. Let's calculate the threshold
nlayers <- nlyr(water)
threshold <- nlayers * 0.1

# Sum up all layers to get the number of times each cell was flooded
water <- sum(water)

# Every cell that was flooded more than x% of the times will be kept for our
# averaged map. Let's prepare the reclassification table for this operation
rcl <- data.frame(
    oldfrom = c(-Inf, threshold)
  , oldto   = c(threshold, Inf)
  , new     = c(0, 1)
)

# Apply the reclassification
water <- classify(water, rcl)

# Visualize
plot(water)

# Store the raster to file
writeRaster(
    water
  , "03_Data/02_CleanData/01_LandCover_WaterCoverAveraged_MERGED.tif"
  , overwrite = TRUE
)

# Use the averaged layer to calculate a distance to water layer
distance <- distanceTo(raster(water), value = 1)

# Plot the resulting layer
plot(distance)

# Store the layer to file
writeRaster(
    distance
  , "03_Data/02_CleanData/01_LandCover_DistanceToWaterAveraged_MERGED.tif"
  , overwrite = TRUE
)
