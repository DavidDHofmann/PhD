################################################################################
#### Combining Different Water Layers
################################################################################
# Description: In this script I combine the water layers from Globeland, ORI,
# OSM, David, and MERIT.

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# load packages
library(tidyverse)  # For data wrangling
library(terra)      # To handle spatial data
library(raster)     # To handle spatial data
library(lubridate)  # To handle dates
library(pbmcapply)  # For multicore use with progress bar

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Merge Layers
################################################################################
# Load the layers we want to merge
water <- rast("03_Data/02_CleanData/LandCover.tif") == 1
river <- rast("03_Data/02_CleanData/Rivers.tif")

# Extract dates
flood_dates <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled" %>%
  dir(path = ., pattern = ".tif$") %>%
  ymd()

# From the floodmaps, we only want to keep those that are closest to some
# dispersal date. So identify unique dates of dispersal
disp_dates <- "03_Data/02_CleanData/Dispersers.csv" %>%
  read_csv() %>%
  subset(State == "Disperser") %>%
  pull(Timestamp) %>%
  as.Date() %>%
  unique()

# Find closest floodmap for each dispersal date
cat("Identifying floodmaps closest to dispersal data... \n")
closest <- lapply(disp_dates, function(x) {
  closest1 <- flood_dates[which(abs(x - flood_dates) == min(abs(x - flood_dates)))][1]
  closest2 <- flood_dates[which(abs(x - flood_dates) == min(abs(x - flood_dates)))][2]
  close <- c(closest1, closest2)
  return(close)
}) %>% do.call(c, .) %>% unique() %>% na.omit()

# Subset to only those floodmaps
flood <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled" %>%
  dir(pattern = ".tif$", full.names  = T) %>%
  data.frame(Filename = ., Date = ymd(basename(.)), stringsAsFactors = F) %>%
  subset(Date %in% closest) %>%
  pull(Filename) %>%
  rast()

# Get their dates
flood_dates <- ymd(names(flood))

# Remove cloud cover (value = 2) from the floodmaps and call it dryland
# This is not too much of an issue since we only require data around the actualy
# dispersal location
flood <- subst(flood, 2, 0)

# Expand maps to match the extent of the study area
cat("Extending floodmaps to the extent of the study area... \n")
flood <- extend(flood, water)

# From the globeland land cover dataset, remove any water within the extent for
# which we have dynamic floodmaps
p <- as.polygons(ext(trim(flood[[1]])))
water <- mask(water, p, inverse = T, updatevalue = 0, touches = F)

# Combine the layer with the with river data
dynamic <- max(water, river)

# Add dynamic floodmaps
cat("Creating dynamic water masks data... \n")
dynamic <- mask(dynamic, flood, maskvalue = 1, updatevalue = 1)

# Assign map dates again
names(dynamic) <- flood_dates

# Plot some of the maps
plot(dynamic[[sample(nlyr(dynamic), size = 4)]], col = c("white", "blue"))

# Save the result to file. We'll store them uncompressed which allows faster
# reading times
writeRaster(dynamic, "03_Data/02_CleanData/WaterCoverDynamic.grd", overwrite = TRUE)

################################################################################
#### Distance To Water
################################################################################
# Reload dynamic watermaps
water <- stack("03_Data/02_CleanData/WaterCoverDynamic.grd")

# Compute distance to water for each watermap
cat("Comptuing distance to dynamic water layers...\n")
distances <- pbmclapply(1:nlayers(water), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  dist <- distanceTo(water[[x]], value = 1)
  dist <- writeRaster(dist, tempfile(fileext = ".tif"))
  names(dist) <- names(water[[x]])
  return(dist)
})

# Convert back to a terra raster
distances <- rast(stack(distances))

# Store them
writeRaster(distances, "03_Data/02_CleanData/DistanceToWaterDynamic.grd", overwrite = T)

################################################################################
#### Create Averaged Watermap
################################################################################
# We also want to create a static watermap. This map basically resembles the
# type of data most people would consider for their analysis. For this, we'll
# create an "average representation of the flood" across the Okavango delta.
water <- rast("03_Data/02_CleanData/LandCover.tif") == 1
river <- rast("03_Data/02_CleanData/Rivers.tif")
flood <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled" %>%
  dir(pattern = ".tif$", full.names  = T) %>%
  rast()

# Reclassify all floodmaps (remove cloud cover in them)
cat("Reclassifying all floodmaps so we can create an averaged watermap...\n")
flood <- subst(flood, 2, 0)

# Sum them
summed <- sum(flood)

# Keep everything that is inundated most of the time
summed <- summed > nlyr(flood) * 0.1

# Store the static floodmap to file
writeRaster(summed, "03_Data/02_CleanData/00_Floodmaps/FloodMapStatic.tif", overwrite = T)

# Extend to the main study area
summed <- extend(summed, water)

# Merge with the other layers
static <- max(water, river)
static <- mask(static, summed, maskvalue = 1, updatevalue = 1)

# Compute distance to water
cat("Comptuing distance to static water...\n")
distance <- distanceTo(raster(static), value = 1)
distance <- rast(distance)

# Visualize the map
plot(static, col = c("white", "cornflowerblue"))
plot(distance)

# Store the file
writeRaster(static, "03_Data/02_CleanData/WaterCoverStatic.tif", overwrite = TRUE)
writeRaster(distance, "03_Data/02_CleanData/DistanceToWaterStatic.tif", overwrite = TRUE)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/06_MergeWater.rds")
cat("Done :)\n")
