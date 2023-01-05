############################################################
#### Preparation of the Africa Shapefile
############################################################
# Description: For plotting purposes it will be nice to have an Africa map as
# background and to show the extent of our study. The file was downloaded from
# http://www.maplibrary.org/library/stacks/Africa/index.htm

# Clean environment
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# load packages
library(rgdal)
library(raster)
library(rgeos)

# Load the Africa shapefile
africa <- readOGR("03_Data/01_RawData/ESRI/Africa.shp")

# Assign correct projection
crs(africa) <- CRS("+init=epsg:4326")

# Check if the geometry is valid
gIsValid(africa)

# Make it valid then
africa <- gBuffer(africa, width = 0, byid = TRUE)

# Check again
gIsValid(africa)

# Now save the projected file to our data folder
writeOGR(africa
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_Africa"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
