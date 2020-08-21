############################################################
#### Preparing Population Density Data
############################################################
# Description: Preparation of the tiles that were downloaded from facebook
# Includes stitching, aggregating and reprojecting the tiles

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(raster)
library(gdalUtils)
library(rgdal)
library(davidoff)

# Make use of multicore abilities
beginCluster()

############################################################
#### Merging and Cropping the Files
############################################################
# Identify all Facebook files
files <- dir(
    path    = "03_Data/01_RawData/FACEBOOK"
  , pattern = "population.*tif$"
  , full.names = T
)

# Load them into a list
pop <- lapply(files, raster)

# Merge all tiles together
pop <- do.call(merge, pop)

# Load the reference shapefile
s <- shapefile("03_Data/02_CleanData/00_General_Shapefile")

# Crop the merged file to our extent
pop <- crop(pop, s)

# Store the result to file
writeRaster(
    pop
  , "03_Data/01_RawData/FACEBOOK/PopulationDensity.tif"
  , overwrite = TRUE
)

# Aggregate the layer to 250m
coarse <- aggregate(pop, fact = round(250 / 30), fun = sum)

# Finally we need to resample the layer. Lets load the reference raster
r250 <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Resample the population density layer to the reference raster
coarse <- resample(coarse, r250, "bilinear")

# Replace NAs with 0s
coarse <- reclassify(coarse, rcl = c(NA, NA, 0))

# Store the result
writeRaster(coarse
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanDensity_Facebook.tif"
  , overwrite = TRUE
)

# End cluster
endCluster()
