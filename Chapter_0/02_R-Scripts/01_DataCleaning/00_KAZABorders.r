################################################################################
#### Preparation of the Kaza Borders
################################################################################
# Description: In this script I import the files as they were downloaded from
# the website of the peace parks foundation
# (https://maps.ppf.org.za/KAZA_ME/public/index.html). The KAZA Shapefile is
# assigned a correct crs and stored properly.

# Clean environment
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(rgdal)
library(raster)

# Show the layers stored in the gdb file
layers <- ogrListLayers("03_Data/01_RawData/KAZA/data.gdb")

# Read the layer
kaza <- readOGR(
    dsn   = "03_Data/01_RawData/KAZA/data.gdb"
  , layer = "kztfca_A1_TransFrontierConservationArea"
)

# Plot the imported layer
plot(kaza)

# The projection is not yet set to WGS84. I therefore reproject it to WGS84
kaza <- spTransform(kaza, CRS("+init=epsg:4326"))

# Let's check the area of the kaza (in km2)
area(kaza) / 1e6

# Store the cleaned shapefile
writeOGR(
    kaza
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_KAZA_KAZA"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
