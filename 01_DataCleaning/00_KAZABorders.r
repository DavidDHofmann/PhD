############################################################
#### Preparation of the Kaza Borders
############################################################
# Description: Import and cleaning of the KAZA borders as downloaded from
# https://maps.ppf.org.za/KAZA_ME/public/index.html.

# Clean environment
rm(list = ls())

# Set working directory
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(rgdal)    # For loading spatial data
library(raster)   # For manipulating spatial data

# Show the layers stored in the gdb file
layers <- ogrListLayers("03_Data/01_RawData/KAZA/data.gdb")

# Read the layer
kaza <- readOGR(
    dsn   = "03_Data/01_RawData/KAZA/data.gdb"
  , layer = "kztfca_A1_TransFrontierConservationArea"
)

# Plot the imported layer
plot(KAZA)

# The projection is not yet set to WGS84. I therefore reproject it to WGS84
KAZA <- spTransform(KAZA, CRS("+init=epsg:4326"))

# Let's check the area of the kaza (in km2)
area(KAZA) / 1000000

# Store the cleaned shapefile
writeOGR(
    KAZA
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_KAZA_KAZA"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
