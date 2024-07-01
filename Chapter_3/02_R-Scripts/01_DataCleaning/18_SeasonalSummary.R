################################################################################
#### Seasonal Summary
################################################################################
# Description: Here, we prepare the landscape layers that we need for our
# seasonal simulation. Specifically, we are going to generate, for each
# covariate, a sequence of layers that represents a typical year. For this,
# we'll use all available layers and aggregate them across years.

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(terra)        # To handle spatial data
library(raster)       # To handle spatial data
library(tidyverse)    # To wrangle data
library(lubridate)    # To handle dates
library(hms)          # To handle timestamps
library(pbmcapply)    # For multicore lapply with progress bar

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Function to convert yday to year 2000 (this is really just necessary as it's
# much easier to work with "day of the year", rather than with the proper date,
# as this alleviates issues with products that are generated for the same day,
# but different dates)
ydate <- function(x) {
  as.Date(x - 1, origin = "2000-01-01")
}

################################################################################
#### NDVI, Precipitation, Temperature, Humans, ForestCover
################################################################################
# Load data and extract timestamps from layernames
covs <- tibble(
    Covariate      = c("NDVI", "Precipitation", "Temperature")
  , InputFilename  = paste0("03_Data/02_CleanData/", Covariate, "Dynamic.tif")
  , OutputFilename = paste0("03_Data/02_CleanData/", Covariate, "DynamicAggregated.tif")
)

# Go through each covariate and aggregate data
covs$Layerdates <- lapply(1:nrow(covs), function(i) {

  # Extract relevant values
  covariate <- covs$Covariate[i]
  infile    <- covs$InputFilename[i]
  outfile   <- covs$OutputFilename[i]

  # Load data
  dat <- rast(infile)

  # Extract dates from layernames (note that we will transform some of the dates
  # to day of the year, and then backtransform to a proper date, yet for the
  # year 2000. This will lead to minor offsets by 1 day for some dates, but will
  # substantially reduce the number of layers, as e.g. 17/18 are not treated as
  # separate days anymore, particularly for MODIS products where the yday
  # matters, not the true date)
  if (covariate %in% c("NDVI")) {
      covariate_dates <- dat %>%
        names() %>%
        substr(start = nchar(.) - 9, stop = nchar(.)) %>%
        ymd(tz = "UTC") %>%
        yday() %>%
        ydate()
    } else if (covariate %in% c("Precipitation", "Temperature")) {
      covariate_dates <- dat %>%
        names() %>%
        substr(start = nchar(.) - 18, stop = nchar(.)) %>%
        ymd_hms(tz = "UTC") %>%
        tibble(Date = as_date(.), Time = as_hms(.)) %>%
        mutate(Date = yday(Date)) %>%
        mutate(Date = ydate(Date)) %>%
        mutate(Timestamp = ymd_hms(paste(Date, Time), tz = "UTC")) %>%
        pull(Timestamp)
    } else {
      covariate_dates <- now() # Date really doesn't matter
  }

  # Put into a nice tibble
  info <- tibble(
        Layer     = 1:nlyr(dat)
      , Timestamp = covariate_dates
    ) %>%
    arrange(Timestamp) %>%
    nest(Layers = -Timestamp)

  # Aggregate across years
  if (!file.exists(outfile)) {
    cat("Aggregating", covariate, "data across years... \n")
    pb <- txtProgressBar(min = 0, max = nrow(info), style = 3)
    info$Raster <- lapply(1:nrow(info), function(j) {
      ind <- info$Layers[[j]]$Layer
      avg <- mean(dat[[ind]])
      setTxtProgressBar(pb, value = j)
      return(avg)
    })

    # Load the data into a single stack and assign timestamps as layernames
    data <- rast(info$Raster)
    names(data) <- info$Timestamp

    # Store the resulting layers to file
    writeRaster(data
      , filename  = outfile
      , overwrite = T
    )

  }

  # Keep track of the layerdates
  return(info$Timestamp)

})

################################################################################
#### Floodmaps
################################################################################
# Load data
covariate <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled/" %>%
  dir(pattern = ".tif", full.names = T) %>%
  rast()

# Extract layerdates
covariate_dates <- names(covariate) %>%
      substr(start = nchar(.) - 9, stop = nchar(.)) %>%
      ymd(tz = "UTC") %>%
      yday() %>%
      ydate()

# Put into a nice tibble
info <- tibble(
      Layer     = 1:nlyr(covariate)
    , Timestamp = covariate_dates
  ) %>%
  arrange(Timestamp) %>%
  nest(Layers = -Timestamp)

# Aggregate by the finest possible temporal resolution
cat("Aggregating floodmaps across years \n")
pb <- txtProgressBar(min = 0, max = nrow(info), style = 3)
info$Raster <- lapply(1:nrow(info), function(x) {
  ind    <- info$Layers[[x]]$Layer
  maps   <- covariate[[ind]]
  maps   <- subst(maps, 2, 0)
  nmaps  <- nlyr(maps)
  thresh <- nmaps / 2
  maps   <- sum(maps)
  map    <- maps > thresh
  setTxtProgressBar(pb, value = x)
  return(map)
})

# Load the static water and river layers
water <- rast("03_Data/02_CleanData/LandCover.tif") == 1
river <- rast("03_Data/02_CleanData/Rivers.tif")

# Expand maps to match the extent of the study area
cat("Extending floodmaps to the extent of the study area... \n")
flood <- rast(info$Raster)
flood <- extend(flood, water)

# From the globeland land cover dataset, remove any water within the extent for
# which we have dynamic floodmaps
p <- as.polygons(ext(trim(flood[[1]])), crs = "+init=epsg:4326")
water <- mask(water, p, inverse = T, updatevalue = 0, touches = F)

# Combine the layer with the with river data
dynamic <- max(water, river)

# Add dynamic floodmaps
cat("Creating dynamic water masks data... \n")
dynamic <- mask(dynamic, flood, maskvalue = 1, updatevalue = 1)

# Assign map dates again
names(dynamic) <- info$Timestamp

# Store the resulting layers to file
writeRaster(dynamic
  , filename  = "03_Data/02_CleanData/WaterDynamicAggregated.tif"
  , overwrite = T
  , gdal      = "COMPRESS=NONE"
)

################################################################################
#### Distance To Water
################################################################################
# Reload dynamic watermaps
water <- rast("03_Data/02_CleanData/WaterDynamicAggregated.tif")

# Compute distance to water for each watermap
cat("Computing distance to dynamic water layers...\n")
distances <- pbmclapply(1:nlyr(water), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  dist <- distanceTo(water[[x]], value = 1)
  dist <- raster(dist)
  dist <- writeRaster(dist, tempfile(fileext = ".tif"))
  names(dist) <- names(water[[x]])
  return(dist)
})

# Put into stack
distances <- stack(distances)
distances <- rast(distances)
names(distances) <- names(dynamic)

# Store them
writeRaster(distances
  , filename  = "03_Data/02_CleanData/DistanceToWaterDynamicAggregated.tif"
  , overwrite = T
  , gdal      = "COMPRESS=NONE"
)

################################################################################
#### Vegetation
################################################################################
# Load data
covariate <- "03_Data/02_CleanData/00_Vegmaps/" %>%
  dir(pattern = ".tif", full.names = T) %>%
  rast()

# Extract layerdates
covariate_dates <- names(covariate) %>%
      substr(start = nchar(.) - 9, stop = nchar(.)) %>%
      ymd(tz = "UTC") %>%
      yday() %>%
      ydate()

# Put into a nice tibble
info <- tibble(
      Layer     = 1:nlyr(covariate)
    , Timestamp = covariate_dates
    , Category  = substr(names(covariate), start = 1, stop = nchar(names(covariate)) - 11)
  ) %>%
  arrange(Timestamp) %>%
  nest(Layers = -c(Timestamp, Category))

# Aggregate by the finest possible temporal resolution
cat("Aggregating Vegetation data \n")
pb <- txtProgressBar(min = 0, max = nrow(info), style = 3)
info$Raster <- lapply(1:nrow(info), function(x) {
  ind <- info$Layers[[x]]$Layer
  avg <- mean(covariate[[ind]])
  setTxtProgressBar(pb, value = x)
  return(avg)
})

# Put into a stack
veget <- rast(info$Raster)
names(veget) <- info$Category

# Separate tree cover from shrub cover
trees_map <- veget[[grepl(names(veget), pattern = "Tree")]]
shrub_map <- veget[[grepl(names(veget), pattern = "Shrub")]]

# Merge with floodmaps
trees_map <- mask(trees_map, dynamic, maskvalue = 1, updatevalue = 0)
shrub_map <- mask(shrub_map, dynamic, maskvalue = 1, updatevalue = 0)

# Take over the layernames
names(trees_map) <- names(dynamic)
names(shrub_map) <- names(dynamic)

# Store them
writeRaster(trees_map
  , filename  = "03_Data/02_CleanData/TreesDynamicAggregated.tif"
  , overwrite = T
  , gdal      = "COMPRESS=NONE"
)
writeRaster(shrub_map
  , filename  = "03_Data/02_CleanData/ShrubsDynamicAggregated.tif"
  , overwrite = T
  , gdal      = "COMPRESS=NONE"
)

# # Add info to the tibble
# covs_add <- tibble(
#     Type       = "Dynamic"
#   , Covariate  = c("Water", "DistanceToWater", "Trees", "Shrubs")
#   , Filename   = c(
#       "03_Data/02_CleanData/00_SeasonalAggregates/Water.tif"
#     , "03_Data/02_CleanData/00_SeasonalAggregates/DistanceToWater.tif"
#     , "03_Data/02_CleanData/00_SeasonalAggregates/Trees.tif"
#     , "03_Data/02_CleanData/00_SeasonalAggregates/Shrubs.tif"
#   )
#   , Dates = list(
#       time(dynamic)
#     , time(dynamic)
#     , time(dynamic)
#     , time(dynamic)
#   )
# )
#
# # Put all together
# covs <- covs %>%
#   select(Type, Covariate, Filename = OutputFilename, Dates) %>%
#   rbind(., covs_add)
#
# # Store the timestamps to file
# write_rds(covs, "03_Data/02_CleanData/00_SeasonalAggregates/Covariates.rds")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/18_SeasonalSummary.rds")
cat("Done :)\n")
