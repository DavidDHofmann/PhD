################################################################################
#### Cleaning and Preparation of the Geofabrik Roads Data
################################################################################
# Description: In this script I clean the road files that I downloaded from
# Geofabrik (www.geofabrik.de), which provides ready to download shapefiles
# originating from the open stree maps project. I also cut down the number of
# roads to include only large tar roads.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)    # To handle raster data
library(terra)     # To handle raster data
library(tidyverse) # For data wrangling
library(davidoff)  # Access to custom functions

# Load road data
roads <- vect("03_Data/01_RawData/GEOFABRIK/Roads.shp")

# Check out the description of the classes as derived from the OSM webpage
# (https://wiki.openstreetmap.org/wiki/Key:highway).
legend <- read_csv2("03_Data/01_RawData/GEOFABRIK/RoadsDescription.csv")

# Keep only the largest roads (1-4) and their links (9-12)
roads <- subset(roads, roads$fclass %in% legend$Value[c(1:4, 9:12)])

# Crop the shapefile to our study extent. I therefore import the blank shapefile
# I created earlier
s <- vect("03_Data/02_CleanData/00_General_Shapefile.shp")
roads_crop <- crop(roads, s)

# Plot the remaining roads
plot(roads)

# Save the cropped shapefile
writeVector(
    x         = roads
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.shp"
  , overwrite = T
)

################################################################################
#### Rasterize Roads
################################################################################
# Let's load the reference raster so that we can rasterize our roads
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")

# Rasterize the roads
roads_r <- rasterize(roads, r, field = 1, background = 0)

# Store the layer
writeRaster(
    x         = raster(roads_r)
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.tif"
  , overwrite = T
)

################################################################################
#### Distance to Roads
################################################################################
# To calculate the distance to roads we need to replace 0s with NAs
roads_r <- classify(roads_r, matrix(c(0, NA), ncol = 2))
roads_r <- raster(roads_r)

# Calculate distance to roads (using our custom function)
distance2roads <- distanceTo(roads_r, value = 1)

# Visualize
plot(distance2roads)

# Store the rasters to file
writeRaster(
    x         = distance2roads
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToRoads_GEOFABRIK.tif"
  , overwrite = TRUE
)
