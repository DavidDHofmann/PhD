################################################################################
#### Land Cover
################################################################################
# Description: Generate floodmaps that represent different seasons

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_10")

# Load required packages
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Water
################################################################################
# Print to terminal
cat("Computing flood summaries... \n")

# Find all floodmaps and their associated date
files <- tibble(
    File  = dir(path = "03_Data/01_RawData/FLOODMAPS", pattern = ".tif$", full.names = T)
  , Date  = ymd(basename(File))
)

# Load one of the files for reference
ref <- rast(files$File[1])

# Let's determine the flood extent in each image
flood_summary <- files$File %>%
  rast() %>%
  freq(bylayer = T) %>%
  as.data.frame() %>%
  mutate(count = count / ncell(ref)) %>%
  pivot_wider(
    , id_cols     = layer
    , names_from  = value
    , values_from = count
    , values_fill = 0
  ) %>%
  rename(Flood = "0", Dryland = "255", Cloud = "127") %>%
  mutate(Date = files$Date)

# Put data into a nice tibble for plotting
flood_summary <- tibble(
    Date    = flood_summary$Date
  , Flood   = flood_summary$Flood
  , Dryland = flood_summary$Dryland
  , Cloud   = flood_summary$Cloud
  , Year    = year(Date)
  , Month   = month(Date)
  , Week    = week(Date)
  , Day     = yday(Date)
)

# Join data
files <- left_join(files, flood_summary, by = "Date")
rm(flood_summary, ref)

# Keep only floodmaps with cloud cover below 5%
files <- subset(files, Cloud < 0.05)

# Let's create an averaged map
cat("Computing averaged floodmap...\n")
avg <- files %>%
  pull(File) %>%
  rast() %>%
  classify(cbind(c(0, 127, 255), c(1, 0, 0))) %>%
  mean()

# Threshold
flood <- avg > 0.5

# Combine the maps with the globeland dataset
globe <- rast("03_Data/01_RawData/GLOBELAND/Water.tif")

# Get polygon of the extent
p <- as.polygons(ext(trim(flood)))

# Mask globeland layer
globe <- mask(globe, p, inverse = T, updatevalue = 0, touches = F)

# Equalize extent, resolution and origin of the floodmaps and the globeland
# layer
flood <- extend(flood, globe)
flood <- disagg(flood, fact = 2)
flood <- resample(flood, globe, method = "near")

# Put them together
water <- mask(globe, flood, maskvalue = 1, updatevalue = 1)

################################################################################
#### Vegetation
################################################################################
# Load vegetation data
cat("Preparing vegetation layers... \n")
trees <- rast("03_Data/01_RawData/MODIS/Trees.tif")
shrub <- rast("03_Data/01_RawData/MODIS/Shrubs.tif")

# Reclassify water and scale to 0, 1
trees <- subst(trees, from = 200, to = 0) / 100
shrub <- subst(shrub, from = 200, to = 0) / 100

# Use the merged watermaps (without the rivers yet) to set vegetation to 0
trees <- mask(trees, water, maskvalue = 1, updatevalue = 0)
shrub <- mask(shrub, water, maskvalue = 1, updatevalue = 0)

################################################################################
#### Distance To Water
################################################################################
# To finalize our water layer, we will add some rivers when the flood is at its
# mean or maximum extent. We'll not add them when the flood is at a minimum
# level
cat("Computing distance to water... \n")
rivers1 <- rast("03_Data/01_RawData/MERIT/Rivers.tif")
rivers2 <- vect("03_Data/02_CleanData/MajorRivers.shp")
rivers2 <- rasterize(rivers2, rivers1)
rivers2 <- subst(rivers2, NA, 0)
rivers  <- max(rivers1, rivers2)
water   <- max(water, rivers)

# We will also calculate the (sqrt) distance to water
dist_water <- distanceTo(water, value = 1)
dist_water <- rast(dist_water, crs = crs(dist_water))
dist_water <- sqrt(dist_water)
names(dist_water) <- "SqrtDistanceToWater"

################################################################################
#### Human Influence
################################################################################
# We also want to add the human influence layer
humans <- rast("03_Data/02_CleanData/HumanInfluence.tif")

################################################################################
#### Store Layers
################################################################################
# Put all layers into a single stack
layers <- c(water, dist_water, trees, shrub, humans)
names(layers) <- c("Water", "SqrtDistanceToWater", "Trees", "Shrubs", "HumansBuff5000")
plot(layers, legend = F)

# Store all
terra::writeRaster(layers
  , "03_Data/02_CleanData/SpatialLayers.tif"
  , overwrite = T
)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/02_LandCover.rds")

# Print to terminal
cat("Done :) \n")
