############################################################
#### Preparation of the River Shapefile (from Dominik)
############################################################
# Description: Cleaing and projecting of rivers that I received from Dominik

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

# Import file
riv <- shapefile("03_Data/01_RawData/DOMINIK/Rivers")

# Save the file
writeOGR(riv
  , dsn       = "03_Data/02_CleanData"
  , layer     = "03_LandscapeFeatures_Rivers_DOMINIK"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

# End cluster
endCluster()
