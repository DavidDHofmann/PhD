############################################################
#### Distance to Water
############################################################
# Description: Averaged Water map

# clear r's brain
rm(list = ls())

# Define the working directories
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(spatstat)
library(maptools)
library(tidyverse)
library(raster)
library(terra)

# Make use of multiple cores
beginCluster()

# Load water layers
water <- rast("03_Data/02_CleanData/01_LandCover_Water(Merged).tif")

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
    oldfrom = c(0, threshold)
  , oldto   = c(threshold, nlayers)
  , new     = c(0, 1)
)

# Apply the reclassification
water <- reclassify(raster(water), rcl)

# Visualize
plot(water)

# Store the raster to file
writeRaster(
    water
  , "03_Data/02_CleanData/01_LandCover_Water_Averaged.tif"
  , overwrite = TRUE
)

# Use the averaged layer to calculate a distance to water layer
waterppp <- water %>%
  rasterToPoints(., fun = function(x){x == 1}, spatial = TRUE) %>%
  spTransform(., CRS("+init=epsg:32734")) %>%
  as(., "ppp")

# Load the reference raster to prepare a distance to water
distance <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Now replace the values of the raster with the distances to the nearest water
# covered cell
values(distance) <- distance %>%
  as(., "SpatialPoints") %>%
  spTransform(., CRS("+init=epsg:32734")) %>%
  as(., "ppp") %>%
  nncross(., waterppp) %>%
  .[["dist"]]

# Plot the resulting layer
plot(distance)

# Store the layer to file
writeRaster(
    distance
  , "03_Data/02_CleanData/01_LandCover_Water_DistanceToWater.tif"
  , overwrite = TRUE
)

# Terminate the cluster
endCluster()
