################################################################################
#### Cleaning Copernicus Land Cover Data
################################################################################
# Description: In this script I'll clean the land cover data for KAZA downloaded
# from here: https://lcviewer.vito.be/download

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load packages
library(terra)      # To handle rasters
library(raster)     # To handle rasters
library(tidyverse)  # For data wrangling

# Load the data
dat <- dir(
    path       = "03_Data/01_RawData/COPERNICUS"
  , full.names = T
  , pattern    = ".tif$"
)
dat <- lapply(dat, rast)

# Merge the tiles
dat <- do.call(merge, dat)

# Crop the data to our study area
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")
dat <- crop(dat, r, snap = "out")

# Store the merged object to file
writeRaster(
    x         = raster(dat)
  , filename  = "03_Data/01_RawData/COPERNICUS/Copernicus.tif"
  , overwrite = T
)

# Load land cover classes
info <- read_csv("03_Data/01_RawData/COPERNICUS/LandCoverClasses.csv")
info <- arrange(info, Code)

# Not all classes are represented in the cropped study area
vals <- freq(dat)
info <- subset(info, Code %in% vals[, 2])

# Prepare reclassification table
info$ClassNew <- info$Class
info$ClassNew[info$Code >= 111 & info$Code <= 126] <- "Forest"
info$ClassNew[info$Code == 0] <- "NA"
info$ClassNew[info$Code == 20] <- "Shrubs"
info$ClassNew[info$Code == 30] <- "Grassland"
info$ClassNew[info$Code == 40] <- "Cropland"
info$ClassNew[info$Code == 50] <- "Urban"
info$ClassNew[info$Code == 60] <- "Bare"
info$ClassNew[info$Code == 80] <- "Water"
info$ClassNew[info$Code == 90] <- "Water"

# Create new codes
info$CodeNew <- NA
info$CodeNew[info$ClassNew == "NA"] <- 0
info$CodeNew[info$ClassNew == "Water"] <- 1
info$CodeNew[info$ClassNew == "Urban"] <- 2
info$CodeNew[info$ClassNew == "Cropland"] <- 3
info$CodeNew[info$ClassNew == "Forest"] <- 4
info$CodeNew[info$ClassNew == "Shrubs"] <- 5
info$CodeNew[info$ClassNew == "Grassland"] <- 6
info$CodeNew[info$ClassNew == "Bare"] <- 7

# Arrange
info <- arrange(info, CodeNew)

# Assign a color to the new classes
info$Color[info$CodeNew == 0] <- "transparent"
info$Color[info$CodeNew == 1] <- "blue"
info$Color[info$CodeNew == 2] <- "red"
info$Color[info$CodeNew == 3] <- "pink"
info$Color[info$CodeNew == 4] <- "darkgreen"
info$Color[info$CodeNew == 5] <- "orange"
info$Color[info$CodeNew == 6] <- "beige"
info$Color[info$CodeNew == 7] <- "grey"

# Reclassify raster
rcl <- dplyr::select(info, c(Code, CodeNew))
dat <- classify(dat, rcl)

# Aggregate to 250 meters
coarse <- aggregate(dat, fact = round(250 / 110), fun = modal)

# Resample to reference raster
coarse <- resample(coarse, r, method = "near")

# Visualize it
plot(raster(coarse), col = unique(info$Color), breaks = 0:8 - 1)

# Store the raster
writeRaster(
    x         = raster(coarse)
  , filename  = "03_Data/02_CleanData/01_LandCover_LandCover_COPERNICUS.tif"
  , overwrite = T
)

# Also store the water-cover layer seperately
water <- coarse == 1
writeRaster(
    x         = raster(water)
  , filename  = "03_Data/02_CleanData/01_LandCover_WaterCover_COPERNICUS.tif"
  , overwrite = T
)

# Store the information table
info %>%
  dplyr::select(, Class = ClassNew, Code = CodeNew, Color) %>%
  distinct() %>%
  write_csv("03_Data/02_CleanData/01_LandCover_LandCover_COPERNICUS.csv")
