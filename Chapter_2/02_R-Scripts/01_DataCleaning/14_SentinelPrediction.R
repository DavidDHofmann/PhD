################################################################################
#### Predicting Pans
################################################################################
# Description: Make predictions on the processed L2A products and merge the
# layers of each month into a single layer

# Clear R's brain
rm(list = ls())

# Load required packages
library(sen2r)         # To translate SAFE format
library(parallel)      # To check for the number of cores
library(tidyverse)     # To wrangle data
library(lubridate)     # To handle dates
library(randomForest)  # For classification
library(pbmcapply)     # To run stuff in parallel
library(sf)            # To handle spatial data
library(terra)         # To handle spatial data
library(gdalUtils)     # To handle spatial data

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

# Load all information on the Sentinel tiles
sent <- read_rds("03_Data/03_Results/99_SentinelJoined.rds")
print(sent)

# Create a temporary directory
tempi <- "/media/david/Elements/Temporary"
dir.create(tempi, showWarnings = F)

# Function to split a raster
splitRas <- function(file, outdir = tempdir(), s = 2, cores = 1) {
  # file <- sources(rb)
  # s <- 5
  # cores <- 8
  # outdir <- tempi

  # Prepare output names
  outnames <- substr(basename(file), start = 1, stop = nchar(basename(file)) - 4)
  outnames <- paste0(outnames, "_Tile_", sprintf("%02d", 1:(s**2)), ".tif")
  outnames <- file.path(outdir, outnames)

  # Get size of each side
  meta <- gdalinfo(file)
  size <- meta[grepl("Size is .*", x = meta)]
  x <- as.numeric(gsub("[^0-9]", "", unlist(strsplit(size, ", "))))[1]
  y <- as.numeric(gsub("[^0-9]", "", unlist(strsplit(size, ", "))))[2]

  # Create windows
  t <- s - 1
  srcwin <- expand.grid(x = 0:t * x / s, y = 0:t * y / s)
  srcwin <- cbind(srcwin, x / s, y / s)
  srcwin <- as.matrix(srcwin)

  # t is nr. of iterations per side
  invisible(mclapply(1:nrow(srcwin), mc.cores = cores, function(x) {
    gdal_translate(
        src_dataset = file
      , dst_dataset = outnames[x]
      , srcwin      = srcwin[x, ]
    )
  }))

  # Return the filenames
  return(outnames)
}

################################################################################
#### Merge Files
################################################################################
# General process for merging the tiles:
# - Load tile
# - Mask problematic pixels
# - Make prediction on that tile (reduces the tile to a single band)
# - Store the tile
# - Pool tiles from the same month and merge them

# Load the random forest classifier that we will use to make pan-predictions.
# Note that we'll use the "Random Forest" classifier for the "Sentinel" data
model <- read_rds("03_Data/03_Results/99_PanMapping.rds")
model <- model$ModelObject[[4]]

# Specify scl values -> Quality control -> Some of these we want to mask
scl <- data.frame(
    Value    = c(0:11)
  , Category = c(
      "NoData"
    , "Saturated or Defective"
    , "Dark Area Pixels"
    , "Cloud Shadows"
    , "Vegetation"
    , "Bare Soils"
    , "Water"
    , "Clouds Low Probability / Unclassified"
    , "Clouds Medium Probability"
    , "Clouds High Probability"
    , "Cirrus"
    , "Snow / Ice"
  )
)

# Specify the values that we want to mask
scl$ValueNew <- ifelse(scl$Category %in%
    c("NoData", "Cloud Shadows", "Clouds Medium Probability", "Clouds High Probability", "Cirrus")
  , yes = 1
  , no  = NA
)

# # Let's only retain the values that we want to mask later
# tomask <- scl %>%
#   subset(Category %in% c("NoData", "Cloud Shadows", "Clouds Medium Probability", "Clouds High Probability", "Cirrus")) %>%
#   pull(Value)

# Let's take a look at the rows that we need to go through. Each row represents
# a collection of files that need to be predicted and merged.
print(sent)

# Let's generate a filename each of the final (i.e. predicted and merged) layers
sent <- mutate(sent, Filename = paste0(
  "03_Data/02_CleanData/00_Panmaps/"
  , Year
  , "_"
  , sprintf("%02d", Month)
  , ".tif"
))

# We might need to create the associated directory
dir.create("03_Data/02_CleanData/00_Panmaps", showWarnings = F)

# Subset to the files that have not been created yet
exists <- file.exists(sent$Filename)
cat(sum(exists), "of", nrow(sent), "Layers already exist(s)...\n")
sent <- subset(sent, !exists)

# Function to clean a file (mask bad quality pixels), crop it to the moving
# window, and make predictions
prepareFile <- function(file, outname, ext, n = 5, ncores = 1) {
  # file <- "/media/david/Elements/L2A/S2A_MSIL2A_20151218T082342_N0201_R121_T34KGD_20151218T084413.SAFE"
  # file <- files_sub$Filename[1]
  # outname <- "/home/david/Schreibtisch"
  # ext <- vect(sent$Window[[1]])
  # n <- 5

  # Translate the SAFE files into virtual rasters. One containing all desired
  # bands, one containing the quality layers
  vrt_band <- suppressMessages(s2_translate(file, outdir = tempi))
  vrt_mask <- suppressMessages(s2_translate(file, prod_type = "SCL", outdir = tempi))

  # Load them (rm for "raster bands", rf for "raster mask")
  rb <- rast(vrt_band)
  rm <- rast(vrt_mask)

  # Reclassify mask values (we want to mask anything that is == 1)
  rm <- classify(rm, as.matrix(scl[, c("Value", "ValueNew")]))

  # Check if the resolution of the two files is the same. If it's not, we need
  # to disaggregate the mask
  if (!all(res(rb) == res(rm))) {
    rm <- disagg(rm, fact = 2, filename = tempfile(fileext = ".tif"))
  }

  # Unless the file lies fully within the moving window, crop it
  ext_proj <- project(ext, rb)
  # within <- as.vector(relate(ext(rb), ext_proj, "within"))
  # if (!within) {
  #   rb <- crop(rb, ext_proj, snap = "out", filename = tempfile(fileext = ".tif"))
  #   rm <- crop(rm, ext_proj, snap = "out", filename = tempfile(fileext = ".tif"))
  # }
  #
  # Split all rasters (bands + mask) into multiple tiles
  rb_tiled <- splitRas(sources(rb), outdir = tempi, s = n, cores = detectCores() - 1)
  rm_tiled <- splitRas(sources(rm), outdir = tempi, s = n, cores = detectCores() - 1)

  # This process generates many xml files which I want to remove. Also remove
  # the temporary files we created earlier.
  file.remove(dir(tempi, ".aux.xml", full.names = T))
  file.remove(sources(rb))
  file.remove(sources(rm))
  rm(rm, rb)
  gc()

  # Run through the tiles, clean them, and run predictions
  predictions <- mclapply(1:length(rb_tiled), mc.cores = ncores, function(z) {

    # Load tile and associated mask
    band_tiled <- rast(rb_tiled[z])
    mask_tiled <- rast(rm_tiled[z])

    # Buffer the mask and apply it
    mask_tiled <- as.polygons(mask_tiled)
    if (length(mask_tiled) > 0) {
        mask_tiled <- buffer(mask_tiled, width = 250)
        masked <- mask(band_tiled, mask_tiled, updatevalue = NA, inverse = T)
      } else {
        masked <- band_tiled
    }

    # Give the bands proper names
    names(masked) <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B11", "B12")

    # Compute indices of interest
    ndvi <- nd(masked, "B8", "B4")
    ndwi <- nd(masked, "B3", "B8")
    ndmi <- nd(masked, "B8", "B11")
    ndsi <- nd(masked, "B3", "B11")
    best <- nd(masked, "B3", "B12")

    # Put all together into a single stack
    all <- c(masked, ndvi, ndwi, ndmi, ndsi, best)
    names(all) <- c(names(masked), "ndvi", "ndwi", "ndmi", "ndsi", "best")
    rm(masked, ndvi, ndwi, ndmi, ndsi, best)
    gc()

    # Make prediction
    prd <- predict(all, model)

    # Reclassify to water & dryland
    prd <- classify(prd, rbind(c(1, 0), c(2, 1), c(3, 1)))

    # Mask everything outside the moving window
    # prd <- mask(prd, ext_proj)

    # Store to file
    filename <- substr(outname, start = 1, stop = nchar(outname) - 4)
    filename <- paste0(filename, "_", sprintf("%02d", z), ".tif")
    writeRaster(prd, filename, overwrite = T)
    rm(all)
    gc()

    # Return the filename
    return(filename)
  })

  # Merge the tiles again
  pred <- lapply(predictions, rast)
  pred <- sprc(pred)
  pred <- mosaic(pred)

  # Remove temporary files
  file.remove(unlist(predictions))
  if (file.exists(vrt_band)) {
    file.remove(vrt_band)
  }
  if (file.exists(vrt_mask)) {
    file.remove(vrt_mask)
  }
  file.remove(rb_tiled)
  file.remove(rm_tiled)

  # Remove ".vrt" folder
  if (dir.exists(file.path(tempi, ".vrt"))) {
    unlink(file.path(tempi, ".vrt"), recursive = T)
  }

  # Store the file
  writeRaster(pred, outname, overwrite = T)
  rm(pred)
  gc()

  # Return the filename
  return(outname)
}

# Check files
print(sent, n = 100)

# Run it on all files
final <- lapply(47:nrow(sent), function(x) {

  # Print some info on the iteration
  cat("Merging files from month", x, "out of", nrow(sent), "...\n")

  # Get the window/extent of the current month
  ext <- suppressWarnings(vect(st_as_sf(sent$Window[[x]])))

  # Get the filenames of the current month and prepare filenames for the cleaned
  # and predicted files
  files <- tibble(
        Filename = sent$data[[x]]$filepath
      , Year     = sent$Year[x]
      , Month    = sent$Month[x]
      , Tile     = sent$data[[x]]$id_tile
    ) %>%
    mutate(
        FileNo        = sprintf("%02d", 1:n())
      , Month         = sprintf("%02d", Month)
      , FilenameFinal = paste0(tempi, "/", Year, "_", Month, "_", FileNo, ".tif")
  )

  # Subset to the files that do not exist yet
  exists <- file.exists(files$FilenameFinal)
  files_sub <- subset(files, !exists)

  # Loop through the files and run cleaning and prediction
  if (nrow(files_sub) > 0) {
    pb <- txtProgressBar(min = 0, max = nrow(files_sub), style = 3)
    predictions <- lapply(1:nrow(files_sub), function(y) {

      # Prepare the file
      pred <- prepareFile(files_sub$Filename[y], files_sub$FilenameFinal[y]
        , ext    = ext
        , n      = 5
        , ncores = 5
      )

      # Print progress
      setTxtProgressBar(pb, y)

      # Return the filename
      return(pred)
    })
    return(predictions)
  }

  # # Generate a reference raster
  # r <- rast(ext, resolution = 10 / 111000)
  # r[] <- NA
  #
  # # Loop through the files and put values onto the reference raster
  # for (i in files$FilenameFinal) {
  #   b <- rast(i)
  #   b <- mask(b, project(ext, crs(b)), updatevalue = NA)
  #   b <- resample(b, r, method = "near")
  #   r <- min(b, r, na.rm = T)
  # }
  #
  # # Store it to file
  # writeRaster(r, sent$Filename[x])

})

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/14_SentinelPrediction.rds")
cat("Done :)\n")
