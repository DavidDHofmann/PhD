############################################################
#### Preparation Flood Animation
############################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
input   <- "/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/00_LargeExtent"
output  <- "/home/david/ownCloud/University/14. FS 19/Masterarbeit/05_Presentation"
setwd(input)

# Set a seed
set.seed(12345)

# Load required packages
library(lubridate)
library(raster)
library(tmap)
library(tidyverse)
library(viridis)
library(rosm)
library(tmaptools)
library(RColorBrewer)
library(animation)

# Load custom functions
source("/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/Functions.r")

############################################################
#### Animation of the Flood
############################################################
# Load required data
setwd("/home/david/Downloads/MODIS/MCD43A4/Output")
maps <- dir(pattern = "*.tif$") %>% stack(., bands = 1)
setwd(input)

# Get the dates from the map descriptions
dates <- substr(names(maps), start = 2, stop = 11) %>% ymd(.)

# We only want to plot relatively cloud free images. Identify the cloud cover in
# each image for this
clouds <- c()
for (i in 1:nlayers(maps)){
  clouds[i] <- sum(values(maps[[i]]) == 127)
  clouds[i] <- clouds[i] / ncell(maps[[i]])
}

############################################################
#### IMPLEMENT SOON
############################################################
# rather than using this loop you can use the freq() function
clouds <- freq(maps, value == 127)

# Let's identify the number of maps with a cloud coverage above 0.1
sum(clouds > 0.1)

# Keep only those images with a cloud coverage below 10%
indices <- which(clouds < 0.1)
maps <- maps[[indices]]
dates <- dates[indices]

# Load our mask to cover any erronous pixels
setwd("/home/david/Downloads/MODIS/MCD43A4/Masks")
mask <- shapefile("MaskNoWater.shp")

# Crop the mask to the extent of the ORI images
mask <- crop(mask, extent(maps))

# Reclassify the pixels
rcl <- data.frame(old = c(0, 127, 255), new = c(1, 0, 0))
maps <- reclassify(maps, rcl)

# Prepare the animation
setwd(output)
ani.options(interval = .1, ani.width = 1920, ani.height = 1080)
saveVideo({
  for (i in 1:nlayers(maps)){
    plot(maps[[i]]
      , col     = viridis(2, end = 0.8)[c(1, 2)]
      , legend  = FALSE
      , axes    = FALSE
      , box     = FALSE
    )
    plot(mask
      , add     = TRUE
      , col     = viridis(2, end = 0.8)[1]
      , border  = viridis(2, end = 0.8)[1]
    )
    text(23.8, -20.5, dates[i], col = "white", cex = 4)
  }
}, video.name = "99_Fluctuations.mp4")
