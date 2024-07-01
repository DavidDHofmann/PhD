################################################################################
#### Computing the Distance to Pans
################################################################################
# Description: In this script, I'll first resample all pan maps to a common
# reference raster. Then, I'll use the resampled maps to generate a composite
# "static" pan map. This map will later allow us to fill some of the gaps in the
# dynamic maps. I'll then generate dynamic "distance to pan".

# Clear R's brain
rm(list = ls())

# Load required packages
library(terra)        # To handle spatial data
library(tidyverse)    # For data-wrangling
library(lubridate)    # To handle dates
library(parallel)     # To work in parallel
library(gdalUtils)    # Some raster manipulation utilities

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Resampling and Creating Static Pan Map
################################################################################
# Create directory if needed
if (!dir.exists("03_Data/02_CleanData/00_Panmaps/Resampled")) {
  dir.create("03_Data/02_CleanData/00_Panmaps/Resampled")
}

# Identify all maps we want to loop through and prepare their output filenames
files <- tibble(
    Filename = dir(path = "03_Data/02_CleanData/00_Panmaps", full.names = T, pattern = ".tif$")
  , Outname  = file.path("03_Data/02_CleanData/00_Panmaps/Resampled", basename(Filename))
  , Date     = substr(basename(Filename), start = 1, stop = nchar(basename(Filename)) - 4)
)

# Let's create a reference raster to which we'll resample the maps
r <- read_rds("03_Data/02_CleanData/Windows.rds") %>%
  pull(Window) %>%
  do.call(rbind, .) %>%
  vect() %>%
  aggregate() %>%
  rast(res = 10 / 111000)

# Loop through the files and resample them to the reference raster
cat("Resampling pan maps to reference raster...\n")
pb <- txtProgressBar(min = 0, max = nrow(files), style = 3)
for (i in seq_len(nrow(files))) {
  if (file.exists(files$Outname[i])) {
    setTxtProgressBar(pb, i)
    next
  }
  map <- rast(files$Filename[i])
  map <- resample(map, r, method = "near", threads = T)
  map <- trim(map) # To save space
  writeRaster(map, files$Outname[i], overwrite = T)
  setTxtProgressBar(pb, i)
}

################################################################################
#### Composite "Static" Pan Map
################################################################################
# Create a composite panmap across all seasons.
if (!file.exists("03_Data/02_CleanData/PansStatic.tif")) {
  cat("Generating composite pan-map...\n")
  pb <- txtProgressBar(min = 0, max = nrow(files), style = 3)

  # Load maps and extend them to the study area. I know, it appears a bit silly
  # # to first trim the maps, only to later expand them, but it allows me to
  # safe a # substantial amount of storage space.
  pans <- lapply(1:nrow(files), function(i) {
    map <- rast(files$Outname[i])
    map <- extend(map, r)
    setTxtProgressBar(pb, i)
    return(map)
  })

  # Compute a composite map
  static        <- rast(pans)
  static        <- mean(static, na.rm = T)
  static_thresh <- static >= 0.5
  static_thresh <- as.numeric(static_thresh)

  # Store it to file
  writeRaster(static_thresh
    , filename  = "03_Data/02_CleanData/PansStatic.tif"
    , overwrite = T
    , gdal      = "COMPRESS=NONE"
  )
}

################################################################################
#### Static "Distance-To-Pan" Map
################################################################################
# The static pan map is rather big, which is why I can't compute the spatial
# distances using terra. Instead, I'll use the gdal_proximity command, which has
# to be executed through the command line. Importantly, the static maps needs
# to projected for this
if (!file.exists("03_Data/02_CleanData/DistanceToPansStatic.tif")) {

  # Load static pan map and derive a mask from the NA values
  # directory)
  pans  <- rast("03_Data/02_CleanData/PansStatic.tif")
  maski <- is.na(pans)

  # Project the map
  file_proj <- tempfile(fileext = ".tif")
  pans_proj <- project(pans, "+init=epsg:32734"
    , method   = "near"
    , filename = file_proj
    , threads  = T
  )

  # Run the gdal command to compute distances
  file_dist <- tempfile(fileext = ".tif")
  command <- paste(
    "gdal_proximity.py -srcband 1 -distunits GEO -values 1 -nodata 0.0 -ot Float32 -of GTiff"
    , file_proj
    , file_dist
  )
  system(command)

  # Now load the distances, resample them, and aggregate the raster to a coarser
  # resolution
  dist <- rast(file_dist)
  dist <- project(dist, pans, method = "near", threads = T)
  dist <- mask(dist, maski, maskvalue = T)
  dist <- aggregate(dist, fact = 10, fun = mean)

  # Write the result to file
  writeRaster(dist
    , filename  = "03_Data/02_CleanData/DistanceToPansStatic.tif"
    , overwrite = T
    , gdal      = "COMPRESS=NONE"
  )

  # Remove temporary files
  file.remove(file_proj)
  file.remove(file_dist)
}

# # The pan map is a massive layer, prohibiting to compute distances in a single
# # go. Hence, I'd like to split it into smaller tiles before computing distances.
# # However, we need to be careful to avoid edge effects, which is why I will
# # buffer the tiles so they overlap. Furthermore, each tile MUST contain some
# # pans, otherwise, it can't compute a distance.
# if (!file.exists("03_Data/02_CleanData/DistanceToPansStatic.tif")) {
#
#   # Split the pan map into tiles (note that the buffer is in pixels)
#   cat("Calculating distance to pans from static layer...\n")
#   buffer <- 20000 / 10
#   pans_split <- splitRas("03_Data/02_CleanData/PansStatic.tif", s = 5, buffer = buffer)
#
#   # Find all pan-maps that are empty (as the respective tile is in a corner) and
#   # remove them
#   empty <- sapply(pans_split, function(x) {
#     r <- rast(x)
#     r <- freq(r)
#     return(nrow(r) == 1)
#   })
#   pans_split <- pans_split[!empty]
#
#   # Now loop through them, compute the distance, and merge them again
#   pb <- txtProgressBar(min = 0, max = length(pans_split), style = 3)
#   dists <- lapply(seq_len(length(pans_split)), function(i) {
#
#     # Get the tile
#     pans_crop <- rast(pans_split[[i]])
#
#     # Compute distances
#     dist <- distanceTo(pans_crop, value = 1, retainNAs = T)
#
#     # Store file
#     outfile <- tempfile(fileext = ".tif")
#     dist <- writeRaster(dist, outfile)
#     setTxtProgressBar(pb, i)
#     return(outfile)
#
#   })
#
#   # The data still has a resolution of 10 meters. I didn't want to coarsen the
#   # data earlier to avoid dilluting the distances. However, at this stage we
#   # can make the data much coarser and reduce the resultion to about 100
#   # meters
#   dist <- lapply(dists, rast) %>% sprc() %>% mosaic(fun = "mean")
#   dist <- aggregate(dist, fact = 10, fun = mean)
#
#   # Write the result to file
#   writeRaster(dist
#     , filename  = "03_Data/02_CleanData/DistanceToPansStatic.tif"
#     , overwrite = T
#     , gdal      = "COMPRESS=NONE"
#   )
#
# }

################################################################################
#### Dynamic "Distance-To-Pan" Maps
################################################################################
# We'll store the distance to maps into a subfolder
if (!dir.exists("03_Data/02_CleanData/00_Panmaps/DistanceTo")) {
  dir.create("03_Data/02_CleanData/00_Panmaps/DistanceTo")
}

# Identify all maps we want to loop through and prepare their output filenames
files <- tibble(
    Filename = dir(path = "03_Data/02_CleanData/00_Panmaps/Resampled", full.names = T, pattern = ".tif$")
  , Outname  = file.path("03_Data/02_CleanData/00_Panmaps/DistanceTo", basename(Filename))
  , Date     = substr(basename(Filename), start = 1, stop = nchar(basename(Filename)) - 4)
)

# Load the moving windows
files <- read_rds("03_Data/02_CleanData/Windows.rds") %>%
  mutate(Date = paste0(Year, "_", sprintf("%02d", Month))) %>%
  ungroup() %>%
  dplyr::select(Window, Date) %>%
  left_join(files, ., by = "Date")

# Also load the "static" pan map. We'll use it to fill NAs
static <- rast("03_Data/02_CleanData/PansStatic.tif")

# Loop through the files and compute distance to the nearest pan
cat("Computing distance to nearest pan...\n")
pb <- txtProgressBar(min = 0, max = nrow(files), style = 3)
for (i in seq_len(nrow(files))) {

  # Extract loop information
  inname  <- files$Filename[i]
  outname <- files$Outname[i]
  window  <- files$Window[[i]]

  # If the file already exists, skip an iteration
  if (file.exists(outname)) {
    setTxtProgressBar(pb, value = i)
    next
  }

  # Load the input map and try to fill any NAs
  pans <- rast(inname)
  pans <- cover(pans, crop(static, pans))

  # Loop through the moving windows (there may be more than just one and we
  # don't want to compute distances for the areas in between them)
  dists <- lapply(seq_len(nrow(window)), function(j) {

    # Crop to the respective window and reproject
    pans_crop <- crop(pans, vect(window[j, ]))

    # Project so we can compute spatial distances
    pans_back <- pans_crop # Keep a backup so we can project stuff back later
    pans_crop <- project(pans_crop, "+init=epsg:32734", method = "near")

    # Get an NA mask for later and then remove NAs
    # pans_mask <- is.na(pans_crop)
    pans_crop <- subst(pans_crop, 0, NA)

    # Compute distance and add back the NAs
    dist <- distance(pans_crop)
    # dist <- mask(dist, pans_mask, maskvalues = T)

    # Reproject
    dist <- project(dist, pans_back, method = "near", threads = T)
    dist <- mask(dist, vect(window[j, ]))

    # The data still has a resolution of 10 meters. I didn't want to coarsen the
    # data earlier to avoid dilluting the distances. However, at this stage we
    # can make the data much coarser and reduce the resultion to about 100
    # meters
    dist <- aggregate(dist, fact = 10, fun = mean)
    dist <- writeRaster(dist, tempfile(fileext = ".tif"))
    return(dist)
  })

  # If there is more than one raster (because there are multiple windows),
  # combine them again
  if (length(dists) > 1) {

      # Get the layersources and define a name for the temporary virtual raster
      tempfiles <- sapply(dists, sources)
      tempvrt   <- tempfile(fileext = ".vrt")

      # Mosaic them
      vrt(tempfiles, tempvrt)
      dists <- rast(tempvrt)
    } else {
      dists <- rast(dists)
  }

  # Store raster to file
  writeRaster(dists, outname, overwrite = T)

  # Remove temporary files and purge other stuff from the tempdir
  remo <- dir(path = tempdir(), pattern = ".tif|.vrt", full.names = T)
  file.remove(remo)

  # Update progress
  setTxtProgressBar(pb, value = i)

}

################################################################################
#### Merge them
################################################################################
# Load all "distance to pan" maps and combine them into a stack. Because
# previously we aggregated them (and they all have a different extent), we need
# to resample them as well so they all match the same reference raster
if (!file.exists("03_Data/02_CleanData/DistanceToPansDynamic.tif")) {
  cat("Merging distance to pan maps...\n")

  # Load all panmaps
  files <- "03_Data/02_CleanData/00_Panmaps/DistanceTo" %>%
    dir(full.names = T, pattern = ".tif$")
  dists <- files %>%
    lapply(rast) %>%
    sprc()

  # Extract layernames
  dates <- substr(basename(files), start = 1, stop = 7)

  # Get an extent and prepare a reference raster
  ext <- ext(dists)
  ref <- rast(ext, res = res(dists[1]))

  # Loop through all files, extend them, and match their sampling to the reference
  pb <- txtProgressBar(min = 0, max = length(dists), style = 3)
  dists <- lapply(1:length(dists), function(i) {
    tmp <- tempfile(fileext = ".tif")
    r <- extend(dists[i], ref)
    r <- resample(r, ref, method = "bilinear", filename = tmp)
    gc()
    setTxtProgressBar(pb, i)
    return(r)
  }) %>% rast()

  # Assign dates as names
  names(dists) <- dates

  # Store resulting layer to file
  writeRaster(dists
    , filename  = "03_Data/02_CleanData/DistanceToPansDynamic.tif"
    , overwrite = T
    , gdal      = "COMPRESS=NONE"
  )
}

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/16_FinalizePanmaps.rds")
cat("Done :)\n")
