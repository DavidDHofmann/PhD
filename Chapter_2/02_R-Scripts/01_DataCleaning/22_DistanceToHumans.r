################################################################################
#### Calculate Distance to Humans
################################################################################
# Description: Here, I calculate layers that indicate the distance to different
# anthropogenic features (humans, human influence, villages, etc.)

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)       # To manipulate spatial rasters
library(terra)        # To manipulate spatial rasters
library(spatstat)     # To calculate distances
library(maptools)     # To convert to psp
library(viridis)      # For nice colors
library(davidoff)     # Custom functions
library(parallel)     # For multicore abilities
library(tidyverse)    # For data wrangling

################################################################################
#### Distance To Humans
################################################################################
# Let's load the layers for which we want to know the distance
files <- c(
    "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence_FACEBOOK.tif"
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence_WORLDPOP.tif"
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_Villages_FACEBOOK.tif"
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_Villages_WORLDPOP.tif"
)

# Load them
humans <- stack(files)

# Make the maps binary
humans <- humans > 0

# Calculate distance on each of the layers
distance <- mclapply(1:nlayers(humans), mc.cores = detectCores() - 1, function(x){
  dist <- distanceTo(humans[[x]], value = 1)
  dist <- writeRaster(dist, tempfile())
  return(dist)
}) %>% stack()

# Plot the result
plot(distance, col = rev(magma(50)))

# Store the rasters to file
writeRaster(
    x         = distance[[1]]
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToHumans_FACEBOOK.tif"
  , overwrite = TRUE
)
writeRaster(
    x         = distance[[2]]
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToHumans_WORLDPOP.tif"
  , overwrite = TRUE
)
writeRaster(
    x         = distance[[3]]
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToVillages_FACEBOOK.tif"
  , overwrite = TRUE
)
writeRaster(
    x         = distance[[4]]
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToVillages_WORLDPOP.tif"
  , overwrite = TRUE
)
