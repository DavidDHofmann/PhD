############################################################
#### Preparation of a Reference Shapefile
############################################################
# Description: Here I simply create a reference shapefile according to which I
# will resample and resize all the other shapefiles. This is comparable to the
# reference raster. I create the reference shapefile according to the reference
# raster due to limited flexibility in adjusting the raster (because of
# resolution).

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(raster)
library(rgdal)

# I will use the raster created in the previous script for this
r <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Create extent from the raster
extent <- extent(r)
s <- as(extent, "SpatialPolygons")
crs(s) <- crs(r)

# Plot the shapefile to make sure that it properly frames the reference raster
plot(r)
plot(s, border = "red", add = T)

# Now it is still only a spatial polygon and to save it we need a spatial
# polygon data frame. To do so we need to add some information. I will simply
# add an ID
ss <- SpatialPolygonsDataFrame(s, data = as.data.frame(1))

# Save the final object
writeOGR(
    ss
  , "03_Data/02_CleanData"
  , "00_General_Shapefile"
  , driver = "ESRI Shapefile"
  , overwrite = TRUE
)
