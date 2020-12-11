################################################################################
#### Identify Villages from Facebook Data
################################################################################
# Description: In this script I'll use the worldpop human density data to
# identify villages

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)     # To handle raster data
library(terra)      # To handle raster data
library(rgdal)      # To handle vector data
library(rgeos)      # To manipulate vector data

# Load the human density estimates from worldpop
dens <- rast("03_Data/01_RawData/FACEBOOK/HumanDensity.tif")

# Aggregate to a coarser resolution (should be equal to worldpop)
dens <- aggregate(dens, fact = 3, fun = sum)

# Remove pixels that are inhabited by less than two people
dens_bin <- dens > 2

# Aggregate to even coarser resolution
dens_bin <- aggregate(dens_bin, fact = 3, fun = max)

# Apply focal filter
dens_bin <- focal(dens_bin, w = 3, fun = modal)

# Remove 0s
dens_bin[dens_bin == 0] <- NA

# Coerce raster to polygon
pols <- as.polygons(dens_bin)

# Store polygons it to a file and then reload
file <- tempfile(fileext = ".shp")
writeVector(pols, file, overwrite = T)
pols <- readOGR(file)

# Apply small buffer
pols <- gBuffer(pols, width = 0.5 / 111)

# Rasterize them
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")
pols <- rasterize(vect(pols), r, field = 1, background = 0)

# Store to file
writeRaster(
    x         = raster(pols)
  , filename  = "03_Data/02_CleanData/04_AnthropogenicFeatures_Villages_FACEBOOK.tif"
  , overwrite = T
)
