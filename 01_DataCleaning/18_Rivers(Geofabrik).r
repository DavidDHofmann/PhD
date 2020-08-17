############################################################
#### Preparation of the River Shapefiles (from Geofrabrik)
############################################################
# Description: Here I cut the rivers and water areas shapefiles downloaded from
# Geofabrik.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load packages
library(raster)
library(rgdal)

# Make use of multicore abilities
beginCluster()

############################################################
#### Rivers
############################################################
# Import file
riv <- shapefile("03_Data/01_RawData/GEOFABRIK/Rivers")

# Crop the data according to the reference shapefile
s <- shapefile("03_Data/02_CleanData/00_General_Shapefile")
riv <- crop(riv, s)

# We only want to keep bodies that are classified as rivers
riv <- subset(riv, fclass == "river")

# Save the file
writeOGR(riv
  , dsn       = "03_Data/02_CleanData"
  , layer     = "03_LandscapeFeatures_Rivers_GEOFABRIK"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

# End cluster
endCluster()
