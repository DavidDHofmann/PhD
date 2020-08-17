############################################################
#### Verification of the Floodmap Algorithm
############################################################
# Description: In this script I validate our classification algorithm by
# comparing images that we classified ourselfes, with pictures classified by ORI

# Clear R's brain
rm(list = ls())

# Specify working directories
input1 <- "/home/david/Schreibtisch/15. PhD/Chapter_0/03_Data/02_CleanData/00_Floodmaps/01_Original"
input2 <- "/home/david/Schreibtisch/15. PhD/Chapter_0/03_Data/02_CleanData/00_Floodmaps/03_Validation"

# Load required packages.
library(tidyverse)
library(raster)
library(lubridate)
library(terra)

# Identify the files of both input locations
files_ori <- dir(path = input1, pattern = "*.tif$")
files_mod <- dir(path = input2, pattern = "*.tif$")
dates_ori <- gsub(files_ori, pattern = "-", replacement = ".")
dates_mod <- gsub(files_mod, pattern = "-", replacement = ".")

# Identify the files which both locations have in common
ori <- files_ori[which(dates_ori %in% dates_mod)]
mod <- files_mod[which(dates_mod %in% dates_ori)]

# Make sure that there are 48 files in total
length(ori)
length(mod)

# Load two stacks, one for the ori files, one for the modis files
setwd(input1)
ori <- rast(ori)
setwd(input2)
mod <- rast(mod)

# We only care about water and dryland, but not about clouds.
# Let's reclassify the images so that water becomes 1, dryland 0 and clouds 0 as
# well
rcl <- data.frame(old = c(0, 127, 255), new = c(1, 0, 0))
ori <- classify(ori, rcl)
mod <- classify(mod, rcl)

# Prepare a function that produces a plot to compare the ori and modis images
compImg <- function(x, y, z){
  par(mfrow = c(1, 2), mar = c(1, 1, 1, 1))
  plot(x[[z]]
    , col = c("white", "blue")
    , legend = FALSE
    , axes = FALSE
    , main = "Prediction"
  )
  plot(y[[z]]
    , col = c("white", "blue")
    , legend = FALSE
    , axes = FALSE
    , main = "Comparison"
  )
}

# Run the function on a desired layer
compImg(mod, ori, 48)

# Prepare a function that produces a false positive or false negative plot
falseMap <- function(x, y, z){
  x[[z]] - y[[z]]
}

# Run the function on some desired layers and plot the result
wrong <- falseMap(mod, ori, 48)

# Plot the differences
plot(wrong
  , col = c("red", "white", "blue")
  , legend = FALSE
  , main = "False Classifications"
)
legend("bottomright"
  , c("False Negatives", "False Positives")
  , pch = c(15, 15)
  , col = c("red", "blue")
)

# Let's identify the number of wrongly classified pixels for each image
false <- c()
for (i in 1:nlyr(ori)){
  n_wrong   <- falseMap(mod, ori, i) %>% abs() %>% global(., sum)
  n_cell    <- ncell(wrong)
  false[i]  <- n_wrong / n_cell
}
false <- do.call(c, false)

# Let's see how many pixels we get correct on average
1 - mean(false)
