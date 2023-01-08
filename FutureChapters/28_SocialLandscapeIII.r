############################################################
#### Preparing the Social Landscape
############################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(data.table)
library(viridis)
library(adehabitatHR)
library(raster)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(rasterVis)
library(rgdal)
library(rgeos)
library(parallel)

# Load custom functions
source("Functions.r")

# Make use of multicore abilities
beginCluster()

# Reload data from Part II
hr_final <- readRDS("99_SocialLandscape(PartII).rds")

############################################################
#### Rasterize Homeranges
############################################################
# We want to rasterize all homeranges. Let's load the previously prepared social
# landscape for this
r <- "03_Data/02_CleanData/05_SocialFeatures_SocialLandscape(UtilisationDistributions)" %>%
  stack() %>%
  .[[1]]

# Let's also prepare layernames
names <- paste(hr_final$Date, hr_final$Pack, sep = "_")

# Split the hr polygons into a list so we can use mclapply
hr_ras <- list()
for (i in 1:length(hr_final)){
  hr_ras[[i]] <- hr_final[i, ]
}

# Rasterize all hrs onto the reference raster and stack them
hr_ras <- mclapply(hr_ras, mc.cores = (detectCores() - 1), function(x){
  y <- rasterize(x, r)
  y <- writeRaster(y, tempfile())
  gc()
  return(y)
}) %>% stack()

# Make sure that the NAs are replaced with 0s
hr_ras <- calc(hr_ras, function(x){
  x[is.na(x)] <- 0
  return(x)
})

# Finally we assign some nice layer names
names(hr_ras) <- names

# Write the stack to file
writeRaster(
    hr_ras
  , filename  = "03_Data/02_CleanData/05_SocialFeatures_SocialLandscape(HomeRanges).tif"
  , format    = "raster"
  , overwrite = TRUE
  , options   = c("INTERLEAVE = BAND", "COMPRESS = LZW")
)

# Remove the rds object from part II
file.remove("99_SocialLandscape(PartII).rds")

# End the cluster
endCluster()
