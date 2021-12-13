################################################################################
#### Combining Different Water Layers
################################################################################
# Description: In this script I combine the water layers from Globeland, ORI,
# OSM, David, and MERIT.

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# load packages
library(tidyverse)  # For data wrangling
library(terra)      # To handle spatial data

# Load the layers we want to merge
flood <- "03_Data/02_CleanData/00_Floodmaps" %>%
  dir(path = ., pattern = ".tif$", full.names = T) %>%
  rast()
water <- rast("03_Data/02_CleanData/01_LandCover_WaterCover.tif")
river <- rast("03_Data/02_CleanData/03_LandscapeFeatures_Rivers.tif")

# Remove cloud cover (value = 2) from the floodmaps and call it dryland
rcl <- data.frame(old = c(0, 1, 2), new = c(0, 1, 0))
flood <- classify(flood, rcl)

# From the globeland land cover dataset, remove any water within the extent for
# which we have dynamic floodmaps
p <- as.polygons(ext(flood[[1]]))
water <- mask(water, p, inverse = T, updatevalue = 0)

# Combine the layer with the with river data
water <- max(water, river)

test <- extend(flood[[1]], water)
filled <- mask(water, test, maskvalue = 1, updatevalue = 1)

# Coerce to raster
water <- raster(water)
flood <- stack(flood)

# "Expand" floodmaps and fill values with values from static watermap
water <- suppressMessages(
  pbmclapply(
      X                  = 1:nlayers(flood)
    , mc.cores           = detectCores() - 1
    , ignore.interactive = T
    , FUN                = function(x){
      extended <- extend(flood[[x]], water, value = NA)
      filled <- mask(water, extended, maskvalue = 1, updatevalue = 1)
      filled <- writeRaster(filled, tempfile())
      return(filled)
  }) %>% stack()
)

# Let's also transfer the layernames
names(water) <- names(flood)

# Visualize some maps
plot(water[[1:4]], col = c("white", "blue"))

# Save the result to file. We'll store them uncompressed which allows faster
# reading times
writeRaster(
    x         = water
  , filename  = "03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd"
  , overwrite = TRUE
  , options   = c("COMPRESSION=NONE")
)
