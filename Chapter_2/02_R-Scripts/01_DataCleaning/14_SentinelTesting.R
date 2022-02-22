################################################################################
#### Making a Few Predictions
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)
library(sen2r)
library(raster)
library(lubridate)
library(sf)
library(terra)

# Function to compute the normalized difference (nd) index of two bands
nd <- function(img, band_x, band_y) {
  x <- img[[band_x]]
  y <- img[[band_y]]
  nd <- (x - y) / (x + y)
  return(nd)
}

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Generate two areas of interest (one for the GOMOTI, one for DAVID'S KINGFOM)
aoi_gomoti <- as(extent(c(23.497701, 23.579514, -19.544147, -19.501840)), "SpatialPolygons")
aoi_davids <- as(extent(c(23.229300, 23.288424, -19.313093, -19.279414)), "SpatialPolygons")
crs(aoi_gomoti) <- "+init=epsg:4326"
crs(aoi_davids) <- "+init=epsg:4326"

# Let's check with which sentinel tiles the two extents overlap
tile1 <- tiles_intersects(st_as_sf(aoi_gomoti))
tile2 <- tiles_intersects(st_as_sf(aoi_davids))
tiles <- c(tile1, tile2)

# Specify the directories to the processed sentinel files
dir_l2a <- "/media/david/Elements/L2A"

# Identify all processed sentinel 2 tiles
files <- dir(path = dir_l2a, include.dirs = T, full.names = T)

# Read metadata from filenames
meta <- safe_getMetadata(files, info = "fileinfo")
meta$Year <- year(meta$sensing_datetime)
meta$Month <- month(meta$sensing_datetime)

# Keep only September and October 2021 for the desired tile
index1 <- meta$Year == 2021
index2 <- meta$Month %in% c(9, 10)
index3 <- meta$id_tile %in% tiles
index <- index1 & index2 & index3
files_sub <- meta[index, ]

# Let's just keep one of the files
files_sub <- files_sub[1, ]

# Translate it
file <- file.path(dir_l2a, files_sub$name)
imag <- s2_translate(file, outdir = tempfile())
mask <- s2_translate(file, prod_type = "SCL", outdir = tempfile())

# Mask cloud cover
rf <- rast(imag)
rm <- rast(mask)

# Need to disaggregate the mask to match the other file
rm <- disagg(rm, fact = 2)
compareGeom(rf, rm)

# Mask has the following classes (check here: https://sen2r.ranghetti.info/articles/outstructure#accessory-layers)
final <- mask(rf, rm, maskvalue = c(0, 3, 8, 9, 10), updatevalue = NA)

# Rename the bands
names(final) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B11", "B12")

# Compute indices of interest
ndvi <- nd(final, "B8", "B4")
ndwi <- nd(final, "B3", "B8")
ndmi <- nd(final, "B8", "B11")
ndsi <- nd(final, "B3", "B11")
best <- nd(final, "B3", "B12")

# Put all together into a single stack
all <- c(final, ndvi, ndwi, ndmi, ndsi, best)
names(all) <- c(names(final), "ndvi", "ndwi", "ndmi", "ndsi", "best")

# Crop to the two areas of interest
plot(gomoti)
gomoti <- crop(all, project(vect(aoi_gomoti), crs(all)), snap = "out")
davids <- crop(all, project(vect(aoi_davids), crs(all)), snap = "out")

# Load the random forest model
library(randomForest)
model <- read_rds("03_Data/03_Results/99_PanMapping.rds")$ModelObject[[4]]
print(model)

# Make prediction
pred_gomoti <- predict(gomoti, model)
pred_davids <- predict(davids, model)

# Visualize
plot(pred_gomoti)
plot(pred_davids)

# Write to file
writeRaster(pred_gomoti, "03_Data/03_Results/99_GomotiPred.tif")
writeRaster(pred_davids, "03_Data/03_Results/99_MbomaPred.tif")
