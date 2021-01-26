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

# Make use of multiple cores
beginCluster()

################################################################################
#### Stitching the Tiles
################################################################################
# Load the data
dat <- dir(
    path       = "03_Data/01_RawData/COPERNICUS"
  , full.names = T
  , pattern    = ".tif$"
)
dat <- lapply(dat, raster)

# Merge the tiles
merged <- do.call(merge, dat)

################################################################################
#### Cropping, Aggregating, and Simplifying the Stitched Raster
################################################################################
# Load the reference shapefile and raster
s <- shapefile("03_Data/02_CleanData/00_General_Shapefile.shp")
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Crop the merged tiles to our reference shapefile
merged <- crop(merged, s)

# Store the merged object to file
writeRaster(
    x         = merged
  , filename  = "03_Data/01_RawData/COPERNICUS/Copernicus.tif"
  , overwrite = T
)

# Aggregate to coarser resolution
fact <- res(r)[1] / res(merged)[1]
coarse <- aggregate(merged, fact = round(fact), fun = modal)

# Check out the distribution of values
freq(coarse, useNA = "ifany")

# Load land cover classes
info <- read_csv("03_Data/01_RawData/COPERNICUS/LandCoverClasses.csv")
info <- arrange(info, Code)

# Not all classes are represented in the cropped study area
vals <- freq(coarse)
info <- subset(info, Code %in% vals[, 1])

# Prepare reclassification table. Note that I'll replace the NAs with grassland.
info$ClassNew <- info$Class
info$ClassNew[info$Code >= 111 & info$Code <= 126] <- "Forest"
info$ClassNew[info$Code == 0] <- "Grassland"
info$ClassNew[info$Code == 20] <- "Shrubs"
info$ClassNew[info$Code == 30] <- "Grassland"
info$ClassNew[info$Code == 40] <- "Cropland"
info$ClassNew[info$Code == 50] <- "Urban"
info$ClassNew[info$Code == 60] <- "Bare"
info$ClassNew[info$Code == 80] <- "Water"
info$ClassNew[info$Code == 90] <- "Water"

# Create new codes
info$CodeNew <- NA
info$CodeNew[info$ClassNew == "NA"] <- 6
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
info$Color[info$CodeNew == 1] <- "blue"
info$Color[info$CodeNew == 2] <- "red"
info$Color[info$CodeNew == 3] <- "pink"
info$Color[info$CodeNew == 4] <- "darkgreen"
info$Color[info$CodeNew == 5] <- "orange"
info$Color[info$CodeNew == 6] <- "beige"
info$Color[info$CodeNew == 7] <- "grey"

# Reclassify raster
rcl <- dplyr::select(info, c(Code, CodeNew))
new <- reclassify(coarse, rcl)

# Resample to reference raster
new <- resample(new, r, method = "ngb")

# Check out the frequency of different values
freq(new, useNA = "ifany")
sum(is.na(values(new)))

# Visualize it
plot(new, col = unique(info$Color), breaks = 0:6)

################################################################################
#### Store Final Raster
################################################################################
# Store the raster
writeRaster(
    x         = new
  , filename  = "03_Data/02_CleanData/01_LandCover_LandCover_COPERNICUS.tif"
  , overwrite = T
)

# Also store the water-cover layer seperately
water <- new == 1
writeRaster(
    x         = water
  , filename  = "03_Data/02_CleanData/01_LandCover_WaterCover_COPERNICUS.tif"
  , overwrite = T
)

# Store the information table
info %>%
  dplyr::select(, Class = ClassNew, Code = CodeNew, Color) %>%
  distinct() %>%
  write_csv("03_Data/02_CleanData/01_LandCover_LandCover_COPERNICUS.csv")

# End cluster
endCluster()
