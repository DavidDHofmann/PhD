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
library(lubridate)  # To handle dates

################################################################################
#### Merge Layers
################################################################################
# Load the layers we want to merge
water <- rast("03_Data/02_CleanData/01_LandCover_LandCover.tif") == 1
river <- rast("03_Data/02_CleanData/03_LandscapeFeatures_Rivers.tif")
plot(water, col = c("white", "cornflowerblue"))

# Extract dates
flood_dates <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled" %>%
  dir(path = ., pattern = ".tif$") %>%
  ymd()

# From the floodmaps, we only want to keep those that are closest to some
# dispersal date. So identify unique dates of dispersal
disp_dates <- "03_Data/02_CleanData/00_General_Dispersers.csv" %>%
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
writeRaster(
    x         = dynamic
  , filename  = "03_Data/02_CleanData/01_LandCover_WaterCoverDynamic.grd"
  , overwrite = TRUE
)

################################################################################
#### Create Averaged Watermap
################################################################################
# We also want to create a static watermap. This map basically resembles the
# type of data most people would consider for their analysis. For this, we'll
# create an "average representation of the flood" across the Okavango delta.
flood <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled" %>%
  dir(pattern = ".tif$", full.names  = T) %>%
  rast()

# Reclassify all floodmaps (remove cloud cover in them)
cat("Reclassifying all floodmaps...\n")
flood <- subst(flood, 2, 0)

# Sum them
summed <- sum(flood)

# Keep everything that is inundated most of the time
summed <- summed > nlyr(flood) * 0.1

# Extend to the main study area
summed <- extend(summed, water)

# Merge with the other layers
static <- max(water, river)
static <- mask(static, summed, maskvalue = 1, updatevalue = 1)

# Visualize the map
plot(static, col = c("white", "cornflowerblue"))

# Store the file
writeRaster(
    x         = static
  , filename  = "03_Data/02_CleanData/01_LandCover_WaterCoverStatic.tif"
  , overwrite = TRUE
)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/06_MergeWater.rds")
cat("Done :)\n")
