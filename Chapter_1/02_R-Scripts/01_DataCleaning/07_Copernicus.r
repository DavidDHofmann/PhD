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
library(terra)
library(raster)
library(sp)
library(tidyverse)

# Load the data
dat <- dir(
    path       = "03_Data/01_RawData/COPERNICUS"
  , full.names = T
  , pattern    = ".tif$"
)
dat <- lapply(dat, rast)

# Merge the tiles
dat <- merge(dat[[1]], dat[[2]], dat[[3]], dat[[4]])

# Crop the data to our study area
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")
dat <- crop(dat, r, snap = "out")

# Load land cover classes
info <- read_csv("03_Data/01_RawData/COPERNICUS/LandCoverClasses.csv")
info <- arrange(info, Code)

# Not all classes are represented in the cropped study area
vals <- freq(dat)
info <- subset(info, Code %in% vals[, 2])

# Prepare reclassification table
info$ClassNew <- info$Class
info$ClassNew[info$Code >= 111 & info$Code <= 126] <- "Forest"
info$ClassNew[info$Code == 90] <- "Permanent water bodies"

# Create new codes
info$CodeNew <- NA
info$CodeNew[info$ClassNew == "No input data available"] <- 0
info$CodeNew[info$ClassNew == "Permanent water bodies"] <- 1
info$CodeNew[info$ClassNew == "Urban / built up"] <- 2
info$CodeNew[info$ClassNew == "Cultivated and managed vegetation/agriculture (cropland)"] <- 3
info$CodeNew[info$ClassNew == "Forest"] <- 4
info$CodeNew[info$ClassNew == "Shrubs"] <- 5
info$CodeNew[info$ClassNew == "Herbaceous vegetation"] <- 6
info$CodeNew[info$ClassNew == "Bare / sparse vegetation"] <- 7

# Arrange
info <- arrange(info, CodeNew)

# Assign a color to the new classes
info$Color[info$CodeNew == 0] <- "transparent"
info$Color[info$CodeNew == 1] <- "blue"
info$Color[info$CodeNew == 2] <- "red"
info$Color[info$CodeNew == 3] <- "beige"
info$Color[info$CodeNew == 4] <- "darkgreen"
info$Color[info$CodeNew == 5] <- "orange"
info$Color[info$CodeNew == 6] <- "grey"
info$Color[info$CodeNew == 7] <- "white"

# Reclassify raster
rcl <- dplyr::select(info, c(Code, CodeNew))
dat <- classify(dat, rcl)

# Visualize it
plot(raster(dat), col = unique(info$Color), breaks = unique(info$CodeNew) - 1)

# Store the raster
writeRaster(
    x         = dat
  , filename  = "03_Data/02_CleanData/03_LandscapeFeatures_LandCover_COPERNICUS.tif"
  , overwrite = T
)

# Store the information table
info %>%
  dplyr::select(, Class = ClassNew, Code = CodeNew, Color) %>%
  distinct() %>%
  write_csv("03_Data/02_CleanData/03_LandscapeFeatures_LandCover_COPERNICUS.csv")
