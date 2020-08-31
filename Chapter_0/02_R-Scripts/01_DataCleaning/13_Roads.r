############################################################
#### Cleaning and Preparation of the Geofabrik Roads Data
############################################################
# Description: In this script I clean the road files that I downloaded from
# Geofabrik (www.geofabrik.de), which provides ready to download shapefiles
# originating from the open stree maps project. I also cut down the number of
# roads to include only large tar roads.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(tidyverse)
library(raster)
library(rgdal)
library(spatstat)
library(maptools)
library(gdalUtils)
library(rgeos)

# Make use of multicores
beginCluster()

# Load road data
roads <- readOGR("03_Data/01_RawData/GEOFABRIK/Roads.shp")

# Crop the shapefile to our study extent. I therefore import the blank shapefile
# I created earlier
s <- readOGR("03_Data/02_CleanData/00_General_Shapefile.shp")
roads_crop <- crop(roads, s)

# Get the names of the columns
names(roads_crop)

# The column that indicates the road type or the road's size is called fclass.
# This is the column I will use to select bigger roads. Let's get the unique
# road classes first
unique(roads_crop$fclass)

# Check out the description of the classes as derived from the OSM webpage
# (https://wiki.openstreetmap.org/wiki/Key:highway). I created an excel from
# their descriptions that can be loaded for this purpose
legend <- read_csv2("03_Data/01_RawData/GEOFABRIK/RoadsDescription.csv")
legend

# Keep only the largest roads (1-4) and their links (9-12)
roads <- subset(roads_crop, roads_crop$fclass %in% legend$Value[c(1:4, 9:12)])

# Plot the selection
plot(roads)

# Save the cropped shapefile
writeOGR(roads
  , dsn       = "03_Data/02_CleanData"
  , layer     = "04_AnthropogenicFeatures_Roads_GEOFABRIK"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

############################################################
#### Rasterize Roads
############################################################
# For the human influence layer we also prepare a rasterized layer for all
# roads. Let's load the reference raster in order to burn the roads into it
r250 <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Replace all values with 0s
values(r250) <- 0

# Write the raster to file so we can burn the protected areas into it
writeRaster(r250
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.tif"
  , overwrite = TRUE
)

# Rasterize the roads using gdal
gdal_rasterize(
    src_datasource  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.shp"
  , dst_filename    = "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.tif"
  , burn            = 1
  , at              = TRUE
)

############################################################
#### Distance to Roads
############################################################
# Merge roads into single line object
roads <- gLineMerge(roads)

# We want to calculate distances in meters
roads <- spTransform(roads, CRS("+init=epsg:32734"))

# Coerce them to a ppp object
roads <- as(roads, "psp")

# We also need a reference raster from which we prepare the distance to roads
# raster
DistanceToRoads <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Replace raster values with distances
values(DistanceToRoads) <- DistanceToRoads %>%
  as(., "SpatialPoints") %>%
  spTransform(., CRS("+init=epsg:32734")) %>%
  as(., "ppp") %>%
  nncross(., roads) %>%
  .[["dist"]]

# Plot the result
plot(DistanceToRoads)

# Store the rasters to file
writeRaster(DistanceToRoads
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToRoads.tif"
  , overwrite = TRUE
)

# End cluster
endCluster()
