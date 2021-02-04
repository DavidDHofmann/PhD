############################################################
#### Preparation Flood Animation
############################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_1")

# Set a seed
set.seed(12345)

# Load required packages
library(terra)
library(raster)
library(tidyverse)
library(davidoff)
library(viridis)
library(animation)

################################################################################
#### Animation of the Flood
################################################################################
# Load required data
maps <- dir(
    path       = "03_Data/02_CleanData/00_Floodmaps/01_Original"
  , pattern    = "*.tif$"
  , full.names = T
) %>% rast(., bands = 1)

# Get the dates from the map descriptions
dates <- as.Date(names(maps))

# We only want to plot relatively cloud free images. Identify the cloud cover in
# each image for this
clouds <- c()
for (i in 1:nlyr(maps)){
  clouds[i] <- sum(values(maps[[i]]) == 127)
  clouds[i] <- clouds[i] / ncell(maps[[i]])
}

# Let's identify the number of maps with a cloud coverage above 0.1
sum(clouds > 0.1)

# Keep only those images with a cloud coverage below 10%
indices <- which(clouds < 0.1)
maps <- maps[[indices]]
dates <- dates[indices]

# Load mask to cover any erronous pixels
mask <- vect("03_Data/01_RawData/MODIS/MCD43A4/02_Masks/MaskNoWater.shp")

# Reclassify the pixels
rcl <- data.frame(old = c(0, 127, 255), new = c(1, 0, 0))
maps <- classify(maps, rcl)

# Mask pixels below the nowater mask
maps <- mask(maps, mask, updatevalue = 0, inverse = T)

# Convert to regular raster
maps <- stack(maps)

# Subset
maps2 <- maps[[seq(1, 100, length.out = 5)]]
for (i in 1:5){
  png(paste0(i, ".png"), width = 1980, height = 1980, bg = "transparent")
  plot(maps2[[i]]
    , col    = c("gray25", "cornflowerblue")[c(1, 2)]
    , legend = FALSE
    , axes   = FALSE
    , box    = FALSE
  )
  dev.off()
}

# Prepare the animation
ani.options(interval = .1, ani.width = 1920, ani.height = 1080)
saveVideo({
  for (i in 1:nlayers(maps)){
    plot(maps[[i]]
      , col    = c("black", "cornflowerblue")[c(1, 2)]
      , legend = FALSE
      , axes   = FALSE
      , box    = FALSE
    )
    text(23.8, -20.5, dates[i], col = "white", cex = 4)
  }
}, video.name = "99_Fluctuations.mp4")
