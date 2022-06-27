################################################################################
#### Simulation
################################################################################
# Description: Simulation of dispersal through seasonal landscapes

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(glmmTMB)        # To handle the movement model

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Load Required Data
################################################################################
# Load the calibrated movement model
mod <- read_rds("03_Data/02_CleanData/MovementModel.rds")
mod <- prepareModel(mod)

# Load the fitted gamma distribution (for sampling step lengths)
sl_dist <- read_rds("03_Data/02_CleanData/GammaDistribution.rds")

# Load the habitat layers
water <- rast("03_Data/02_CleanData/WaterCover.tif")
trees <- rast("03_Data/02_CleanData/TreeCover.tif")
shrub <- rast("03_Data/02_CleanData/ShrubCover.tif")
human <- rast("03_Data/02_CleanData/HumanInfluence.tif")
prote <- vect("03_Data/02_CleanData/Protected.shp")
areas <- vect("03_Data/02_CleanData/SourceAreas.shp")

################################################################################
#### Simulate Dispersal
################################################################################
# Simulation parameters

# Source areas

################################################################################
#### Flood
################################################################################
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

# Visualize
ggplot(files, aes(x = Day, y = Flood, col = Year, group = Year)) +
  geom_line() +
  theme_minimal() +
  scale_color_viridis_c() +
  ggtitle("Flood Extent")

# How many files do we have for each month?
plot(table(month(files$Date)))

# Let's pull the 50 maps with the lowest flood extent and the 50 maps with the
# highest flood extent
min <- files %>%
  arrange(Flood) %>%
  slice(1:100) %>%
  pull(File) %>%
  rast() %>%
  classify(cbind(c(0, 127, 255), c(1, 0, 0))) %>%
  mean()
max <- files %>%
  arrange(desc(Flood)) %>%
  slice(1:100) %>%
  pull(File) %>%
  rast() %>%
  classify(cbind(c(0, 127, 255), c(1, 0, 0))) %>%
  mean()

# Let's also create an averaged map
avg <- files %>%
  pull(File) %>%
  rast() %>%
  classify(cbind(c(0, 127, 255), c(1, 0, 0))) %>%
  mean()

# Threshold
min <- min > 0.5
max <- max > 0.5
avg <- avg > 0.5

# Put the maps into a stack
flood <- c(min, avg, max)
names(flood) <- c("min", "mean", "max")

# Visualize
plot(flood, col = c("white", "cornflowerblue"))

# Combine the maps with the globeland dataset
globe <- rast("03_Data/01_RawData/GLOBELAND/Water.tif")

# Get polygon of the extent
p <- as.polygons(ext(trim(flood[[1]])))

# Mask globeland layer
globe <- mask(globe, p, inverse = T, updatevalue = 0, touches = F)

# Equalize extent, resolution and origin of the floodmaps and the globeland
# layer
flood <- extend(flood, globe)
flood <- disagg(flood, fact = 2)
flood <- resample(flood, globe, method = "near")

# Put them together
water <- mask(globe, flood, maskvalue = 1, updatevalue = 1)
names(water) <- c("min", "mean", "max")

# Visualize
plot(water, col = c("white", "cornflowerblue"))

################################################################################
#### Vegetation
################################################################################
# Load vegetation data
trees <- rast("03_Data/01_RawData/MODIS/Trees.tif")
shrub <- rast("03_Data/01_RawData/MODIS/Shrubs.tif")

# Reclassify water and scale to 0, 1
trees <- subst(trees, from = 200, to = 0) / 100
shrub <- subst(shrub, from = 200, to = 0) / 100

# Use the merged watermaps (without the rivers yet) to set vegetation to 0
trees <- mask(trees, globe, maskvalue = 1, updatevalue = 0)
shrub <- mask(shrub, globe, maskvalue = 1, updatevalue = 0)
names(trees) <- names(water)
names(shrub) <- names(water)

# Now add the Merit rivers
merit <- rast("03_Data/01_RawData/MERIT/Rivers.tif")
water <- max(water, merit)

# Store all
writeRaster(trees, "03_Data/02_CleanData/TreeCover.tif")
writeRaster(shrub, "03_Data/02_CleanData/ShrubCover.tif")
writeRaster(water, "03_Data/02_CleanData/WaterCover.tif")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/03_Simulation.rds")
