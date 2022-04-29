################################################################################
#### A Few Pan Predictions
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)    # To wrangle data
library(sen2r)        # To handle sentinel data
library(raster)       # To handle spatial data
library(sf)           # To handle spatial data
library(terra)        # To handle spatial data
library(lubridate)    # To handle dates
library(randomForest) # To handle random forest models

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

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
files <- list.dirs(path = dir_l2a, full.names = T, recursive = F)

# Read metadata from filenames
meta <- safe_getMetadata(files, info = "nameinfo")
meta$filepath <- files

# Create proper timestamps
meta <- meta %>%
  mutate(
      Timestamp = ymd_hms(sensing_datetime)
    , Year      = year(Timestamp)
    , Month     = month(Timestamp)
  )

# Keep only September and October 2021 files
files_sub <- subset(meta, Year == 2021 & Month %in% c(9, 10) & id_tile %in% tiles)

# Let's just keep one of the files
files_sub <- files_sub$filepath[1]

# Translate it
imag <- s2_translate(files_sub, outdir = tempfile())
mask <- s2_translate(files_sub, prod_type = "SCL", outdir = tempfile())

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
gomoti <- crop(all, project(vect(aoi_gomoti), crs(all)), snap = "out")
davids <- crop(all, project(vect(aoi_davids), crs(all)), snap = "out")

# Load the random forest model
model <- read_rds("03_Data/03_Results/99_PanMapping.rds")
model <- model$ModelObject[[4]]

# Make prediction
pred_gomoti <- predict(gomoti, model)
pred_davids <- predict(davids, model)

# Visualize
plot(pred_gomoti)
plot(pred_davids)

# Write to file
writeRaster(pred_gomoti, "03_Data/03_Results/99_GomotiPred.tif")
writeRaster(pred_davids, "03_Data/03_Results/99_MbomaPred.tif")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/15_SentinelTesting.rds")
cat("Done :)\n")
