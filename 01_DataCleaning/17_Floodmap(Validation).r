############################################################
#### Verification of the Floodmap Algorithm
############################################################
# Description: In this script I validate our classification algorithm by
# comparing images that we classified ourselfes, with pictures classified by ORI

# Clear R's brain
rm(list = ls())

# Specify working directories
input1 <- "/home/david/Downloads/MODIS/MCD43A4/Output"
input2 <- "/home/david/Downloads/MODIS/MCD43A4/Validation"

# Load required packages.
library(tidyverse)
library(raster)
library(lubridate)

# Make use of multicore abilities
beginCluster()

# Identify the files of both input locations
setwd(input1)
ori <- dir(pattern = "*.tif$")
setwd(input2)
mod <- dir(pattern = "*.tif$")

# Identify the files which both locations have in common
files <- intersect(ori, mod)

# Make sure that there are 48 files in total
length(files)

# Load two stacks, one for the ori files, one for the modis files
setwd(input1)
ori <- stack(files)
setwd(input2)
mod <- stack(files)

# Make sure that the ori images have the same origin and extent.
# Let's resample the ORI images to the modis images
# ori <- resample(ori, mod, method = "ngb") %>% stack()

# We only care about water and dryland, but not about clouds.
# Let's reclassify the images so that water becomes 1, dryland 0 and clouds 0 as
# well
rcl <- data.frame(old = c(0, 127, 255), new = c(1, 0, 0))
ori <- reclassify(ori, rcl)
mod <- reclassify(mod, rcl)

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
for (i in 1:nlayers(ori)){
  n_wrong   <- falseMap(mod, ori, i) %>% abs() %>% cellStats(., sum)
  n_cell    <- ncell(wrong)
  false[i]  <- n_wrong / n_cell
}

# Let's see how many pixels we get correct on average
1 - mean(false)

# Terminate the cluster
endCluster()
