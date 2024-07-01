################################################################################
#### Combining Different Water Layers
################################################################################
# Description: In this script I combine the water layers from Globeland, ORI,
# OSM, David, and MERIT.

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# load packages
library(tidyverse)  # For data wrangling
library(terra)      # To handle spatial data
library(raster)     # To handle spatial data
library(lubridate)  # To handle dates
library(pbmcapply)  # For multicore use with progress bar
library(spatstat)   # To compute distances

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Merge Layers
################################################################################
# Only do this if the file does not already exist
if (!file.exists("03_Data/02_CleanData/WaterDynamic.tif")) {

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
  # plot(dynamic[[sample(nlyr(dynamic), size = 4)]], col = c("white", "blue"))

  # Save the result to file. We'll store them uncompressed which allows faster
  # reading times
  writeRaster(dynamic
    , filename  = "03_Data/02_CleanData/WaterDynamic.tif"
    , overwrite = TRUE
    , gdal      = "COMPRESS=NONE"
  )

}

################################################################################
#### Distance To Water
################################################################################
# Only do this if the layer does not exist
if (!file.exists("03_Data/02_CleanData/DistanceToWaterDynamic.tif")) {

  # Reload dynamic watermaps
  water <- rast("03_Data/02_CleanData/WaterDynamic.tif")

  # Compute distance to water for each watermap
  # cat("Computing distance to dynamic water layers...\n")
  # distances <- pbmclapply(1:nlayers(water), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  #   dist <- distanceTo(water[[x]], value = 1)
  #   dist <- writeRaster(dist, tempfile(fileext = ".tif"))
  #   names(dist) <- names(water[[x]])
  #   return(dist)
  # })

  cat("Computing distance to dynamic water layers...\n")
  pb <- txtProgressBar(min = 0, max = nlyr(water), style = 3)
  distances <- lapply(1:nlyr(water), function(x) {
    water_proj  <- project(water[[x]], "epsg:32734", method = "near")
    water_proj  <- subst(water_proj, 0, NA)
    dist        <- distance(water_proj)
    dist        <- project(dist, "epsg:4326", method = "near")
    dist        <- crop(dist, water[[x]])
    dist        <- writeRaster(dist, tempfile(fileext = ".tif"))
    names(dist) <- names(water[[x]])
    setTxtProgressBar(pb, x)
    return(dist)
  })

  # Convert back to a terra raster and store
  distances <- rast(distances)
  writeRaster(distances
    , filename  = "03_Data/02_CleanData/DistanceToWaterDynamic.tif"
    , overwrite = T
    , gdal      = "COMPRESS=NONE"
  )

}

################################################################################
#### Create Averaged Watermap
################################################################################
# We also want to create a static watermap. This map basically resembles the
# type of data most people would consider for their analysis. For this, we'll
# create an "average representation of the flood" across the Okavango delta.

# Load required layers
water <- rast("03_Data/02_CleanData/LandCover.tif") == 1
river <- rast("03_Data/02_CleanData/Rivers.tif")
flood <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled" %>%
  dir(pattern = ".tif$", full.names  = T) %>%
  rast()

# Only run this if the layer doesn't exist yet
if (!file.exists("03_Data/02_CleanData/00_Floodmaps/FloodMapStatic.tif")) {

  # Reclassify all floodmaps (remove cloud cover in them)
  cat("Reclassifying all floodmaps so we can create an averaged watermap...\n")
  flood <- subst(flood, 2, 0)

  # Sum them (and write this to file as it takes ages)
  summed <- sum(flood)
  writeRaster(summed, "03_Data/02_CleanData/00_Floodmaps/FloodMapSummed.tif")

  # Keep everything that is inundated most of the time
  summed <- summed > nlyr(flood) * 0.5

  # Store the static floodmap to file
  writeRaster(summed
    , "03_Data/02_CleanData/00_Floodmaps/FloodMapStatic.tif"
    , overwrite = T
    , gdal      = "COMPRESS=NONE"
  )

}

# Check again and skip if all files already exist
if (!file.exists("03_Data/02_CleanData/WaterStatic.tif") | !file.exists("03_Data/02_CleanData/DistanceToWaterStatic.tif")) {

  # Reload static floodmap
  summed <- rast("03_Data/02_CleanData/00_Floodmaps/FloodMapStatic.tif")

  # Extend to the main study area
  summed <- extend(summed, water)

  # Merge with the other layers
  static <- max(water, river)
  static <- mask(static, summed, maskvalue = 1, updatevalue = 1)

  # Compute distance to water
  cat("Computing distance to static water...\n")
  dist <- distanceTo(static, value = 1)

  # Visualize the map
  # plot(static, col = c("white", "cornflowerblue"))
  # plot(dist)

  # Store the file
  writeRaster(static
    , "03_Data/02_CleanData/WaterStatic.tif"
    , overwrite = TRUE
    , gdal      = "COMPRESS=NONE"
  )
  writeRaster(dist
    , "03_Data/02_CleanData/DistanceToWaterStatic.tif"
    , overwrite = TRUE
    , gdal      = "COMPRESS=NONE"
  )

}

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/06_MergeWater.rds")
cat("Done :)\n")
