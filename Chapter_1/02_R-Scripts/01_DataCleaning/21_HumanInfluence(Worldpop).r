################################################################################
#### Human Influence Layer
################################################################################
# Description: Here I combine the layers from Worldpop, Geofabrik, Croplands,
# Gabriele and Globelands to create two layers. First, a binary human influence
# layer that depicts in each cell whether there are humans or not. Second, I
# create a human density map, indicating whether human influence is strong or
# weak.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)       # To manipulate spatial rasters
library(terra)        # To manipulate spatial rasters
library(rgeos)        # To manipulate vector data
library(rgdal)        # To handle vector data
library(viridis)      # For nice colors
library(pbmcapply)    # For multicore with progress bar
library(tidyverse)    # For data wrangling

################################################################################
#### Human Base and Density Layers
################################################################################
# Import data representing human influence
population  <- rast("03_Data/02_CleanData/04_AnthropogenicFeatures_HumanDensity_WORLDPOP.tif")
farms_gabs  <- rast("03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_GABRIELE.tif")
farms_crops <- rast("03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_CROPLANDS.tif")
farms_glob  <- rast("03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_GLOBELAND.tif")
roads       <- rast("03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.tif")

# Combine all the layers (make sure that the farms are not double counted)
merged <- population + roads + max(farms_crops, farms_gabs, farms_glob)

# Let's check how the resulting values are distributed (we only care about
# values larger than 0)
summary(values(merged)[values(merged) > 0])
hist(values(merged)[values(merged) > 0])

# Apparently there are some heavy outliers. Let's remove values beyond 50 (set
# them to 50)
values(merged)[values(merged) > 50] <- 50

# Look at the histogram again
hist(values(merged)[values(merged) > 0])

# Load the shapefiles for the areas in which we know there are only camps but no
# human influence otherwise
delete <- vect("03_Data/01_RawData/DAVID/AnthropogenicInfluence(Delete).shp")

# Use the shapefile to delete the buildings that are inside the polygon
merged <- mask(merged, delete, inverse = TRUE, updatevalue = 0)

# Coerce to regular raster
merged <- raster(merged)

# Write the layer to file
writeRaster(
    x         = merged
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence_WORLDPOP.tif"
  , overwrite = TRUE
)

################################################################################
#### Buffered Human Density
################################################################################
# Define radii of the buffers that we are going to apply
radii <- seq(2.5, 20, by = 2.5)

# Create buffered rasters for the different radii
humans_buff <- pbmclapply(
    X                  = radii
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(x){
    w <- focalWeight(merged, d = x / 111, type = "circle")
    buffered <- focal(merged, w = w / max(w), FUN = sum, pad = T, padValue = 0)
    names(buffered) <- paste0("Buffer_", x * 1000)
    buffered <- writeRaster(buffered, tempfile())
    return(buffered)
}) %>% do.call(stack, .)

# Add an unbuffered layer to the stack (so that we get 0-20 km buffers)
humans_buff <- stack(merged, humans_buff)
names(humans_buff)[1] <- "Buffer_0000"

# Take the log of the layer (add one to assure that we only get positive
# numbers)
humans_buff <- log(humans_buff + 1)

# Plot the result
plot(humans_buff[[1:8]], col = viridis(50))

# Write the buffered layers to file
writeRaster(
    x         = humans_buff
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluenceBuff_WORLDPOP.tif"
  , overwrite = TRUE
)