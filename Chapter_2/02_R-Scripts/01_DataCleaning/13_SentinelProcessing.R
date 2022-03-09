################################################################################
#### Sentinel Data Processing L1C to L2A
################################################################################
# Description: Processing the Sentinel 2 L1C data to achieve L2A products.

# Clear R's brain
rm(list = ls())

# Load required packages
library(sen2r)         # To automate the correction of Sentinel 2 data
library(parallel)      # To check for the number of cores
library(tidyverse)     # To wrangle data
library(lubridate)     # To handle dates
library(terra)         # To handle spatial data
library(randomForest)  # For classification
library(parallel)      # To run stuff in parallel
library(pbmcapply)     # To run stuff in parallel

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Specify the directories to the sentinel files
dir_l1c <- "/media/david/Elements/L1C"
dir_l2a <- "/media/david/Elements/L2A"

# Create a temporary directory
tempi <- "/media/david/Elements/Temporary"
dir.create(tempi, showWarnings = F)

# Function to compute the normalized difference (nd) index of two bands
nd <- function(img, band_x, band_y) {
  x <- img[[band_x]]
  y <- img[[band_y]]
  nd <- (x - y) / (x + y)
  return(nd)
}

################################################################################
#### Process TOA to BOA
################################################################################
# Check out all data that needs to be adjusted (all from L1C)
files <- dir(
    path         = dir_l1c
  , include.dirs = T
  , full.names   = T
  , pattern      = ".SAFE$"
)

# Put files into groups of 15
group <- rep(1:ceiling(length(files) / 15), length.out = length(files))

# Function to determine the name of the corrected product
correctedName <- function(x) {
  corr <- gsub(basename(x), pattern = "MSIL1C", replacement = "MSIL2A")
  return(corr)
}

# Go through the groups, move the files to the computer and run the correction
lapply(group, function(x) {

  # Copy the files
  cat("Copying files to computer...\n")
  files_sub <- files[group == x]
  files_new <- file.path(tempdir(), basename(files_sub))
  file.copy(files_sub, tempdir(), recursive = T)

  # Run the correction
  cat("Running the correction...\n")
  sen2cor(files_new
    , outdir    = dir_l2a
    , parallel  = T
    , overwrite = F
    , use_dem   = F
  )

  # Make sure that all files have been correctly converted
  success <- all(file.exists(file.path(dir_l2a, correctedName(files_sub))))

  # If all files have been correctly converted, remove the originals as well as
  # the copies in the tempdir
  unlink(files_new, recursive = T)
  unlink(files_sub, recursive = T)

})

################################################################################
#### Getting Metadata
################################################################################
# Load paths to all downloaded and processed files
files <- list.dirs(path = dir_l2a, full.names = T, recursive = F)

# We now want to collect metainformation on each of the downloaded files. Let's
# write a function that can do this on multiple cores
getMeta <- function(files) {

  # Get metadata
  meta <- pbmclapply(files
    , mc.cores           = detectCores() - 1
    , ignore.interactive = T
    , FUN                = function(x) {
      meta <- safe_getMetadata(x, abort = F)
      meta$filepath <- x
      meta <- select(meta, filepath, everything())
      return(meta)
    }
  )

  # Do some cleaning
  meta <- meta %>%
    do.call(rbind, .) %>%
    as_tibble() %>%
    mutate(
        Timestamp = ymd_hms(sensing_datetime)
      , Year      = year(Timestamp)
      , Month     = month(Timestamp)
      , clouds    = as.numeric(clouds)
    )

  # Return metainformation
  return(meta)
}

# Note that the function will also check for erronous files. Thus, we can use it
# for a "Quality Control". However, getting the metadata on all files takes
# quite some time and I do not necessarily want to run it over and over again.
# Hence, I'll store the results to file so that we can reload them.
if (file.exists("03_Data/03_Results/99_SentinelMetadata.rds")) {
    meta <- read_rds("03_Data/03_Results/99_SentinelMetadata.rds")
    toadd <- !(files %in% meta$filepath)
    if (length(files[toadd]) > 0) {
      toadd <- getMeta(files[toadd])
      meta <- rbind(meta, toadd)
    }
  } else {
    meta <- getMeta(files)
}

# Store the metadata to file
write_rds(meta, "03_Data/03_Results/99_SentinelMetadata.rds")
meta <- read_rds("03_Data/03_Results/99_SentinelMetadata.rds")

# Ensure that all files are valid
table(meta$validname)

# Nest the data by year and month so that we can join it with the moving windows
# for each month
meta <- meta %>%
  group_by(Year, Month) %>%
  nest() %>%
  mutate(NumberFiles = map_int(data, nrow))

# Load moving windows that we generated based on the GPS data. We will later use
# these to crop the sentinel layers
sent <- "03_Data/03_Results/99_SentinelResults.rds" %>%
  read_rds() %>%
  select(-c(data, From, To, Files, FilesMetaData, NumberFiles))

# Join the metadata with all moving windows. This will allow us to crop each
# tile to the moving window that we want to retain for that month
sent <- left_join(sent, meta, by = c("Year", "Month"))

# Store this to file
write_rds(sent, "03_Data/03_Results/99_SentinelJoined.rds")
sent <- read_rds("03_Data/03_Results/99_SentinelJoined.rds")

# Cleanup
rm(correctedName, dir_l1c, dir_l2a, meta, getMeta, group, toadd, files)
gc()

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
gc()

# Specify scl values -> Quality control -> Some of these we want to mask
scl <- data.frame(
    Value    = c(0:11)
  , Category = c("NoData", "Saturated or Defective", "Dark Area Pixels", "Cloud Shadows", "Vegetation", "Bare Soils", "Water", "Clouds Low Probability / Unclassified", "Clouds Medium Probability", "Clouds High Probability", "Cirrus", "Snow / Ice")
)

# Let's only retain the values that we want to mask later
tomask <- scl %>%
  subset(Category %in% c("NoData", "Cloud Shadows", "Clouds Medium Probability", "Clouds High Probability", "Cirrus")) %>%
  pull(Value)

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
sent <- subset(sent, !exists)
print(sent)

# Function to clean a file and make predictions
prepareFile <- function(file, outname, ext) {
  # file <- file_in
  # outname <- file_out

  # Translate the SAFE files into virtual rasters. One containing all desired
  # bands, one containing the quality layers
  band <- suppressMessages(s2_translate(file, outdir = tempdir()))
  mask <- suppressMessages(s2_translate(file, prod_type = "SCL", outdir = tempdir()))

  # Load them (rm for "raster bands", rf for "raster mask")
  rb <- rast(band)
  rm <- rast(mask)

  # Check if the resolution of the two files is the same. If it's not, we need
  # to disaggregate the mask
  same <- compareGeom(rb, rm, stopOnError = F)
  if (!same) {
    rm <- disagg(rm, fact = 2)
  }

  # If tiles intersect with the boundary of the moving window, crop them to the
  # extent of the moving window of the current month
  ext_proj <- project(ext, rb)
  crosses <- relate(ext(rb), ext_proj, "crosses") %>% as.vector()
  if (crosses) {
    rb <- crop(rb, ext_proj, snap = "out")
    rm <- crop(rm, ext_proj, snap = "out")
  }

  # Generate a grid according to which we can split the tiles
  tl <- rast(nrows = 5, ncols = 5, extent = ext(rb), crs = crs(rb))

  # Generate filenames for the split rasters (we will put them into a temporary
  # folder)
  names_rb <- lapply(1:ncell(tl), function(x) {tempfile(fileext = ".tif")})
  names_rm <- lapply(1:ncell(tl), function(x) {tempfile(fileext = ".tif")})

  # Split all rasters accordingly and store them under said filenames
  rb_tiled <- makeTiles(rb, tl, filename = names_rb, overwrite = T)
  rm_tiled <- makeTiles(rm, tl, filename = names_rm, overwrite = T)
  rm(rb, rm)
  gc()

  # Run through the tiles, clean them, and run predictions
  predictions <- mclapply(1:length(rb_tiled), mc.cores = detectCores() / 2, function(z) {

    # Load tile and associated mask
    band_tiled <- rast(rb_tiled[z])
    mask_tiled <- rast(rm_tiled[z])

    # Mask has the following classes (check here:
    # https://sen2r.ranghetti.info/articles/outstructure#accessory-layers)
    masked <- mask(band_tiled, mask_tiled, maskvalue = tomask, updatevalue = NA)

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
    rm(ndvi, ndwi, ndmi, ndsi, best)
    gc()

    # Make prediction
    prd <- predict(all, model)

    # Store to file
    filename <- substr(outname, start = 1, stop = nchar(outname) - 4)
    filename <- paste0(filename, "_", sprintf("%02d", z), ".tif")
    writeRaster(prd, filename, overwrite = T)

    # Return the filename
    return(filename)
  })

  # Merge the tiles again
  pred <- lapply(predictions, rast)
  pred <- sprc(pred)
  pred <- mosaic(pred)

  # Remove temporary files
  file.remove(unlist(predictions))
  file.remove(band)
  file.remove(mask)
  file.remove(rb_tiled)
  file.remove(rm_tiled)

  # Store the file
  writeRaster(pred, outname, overwrite = T)
  return(outname)
}

# Run it on all files
final <- lapply(1:nrow(sent), function(x) {

  # Print some info on the iteration
  cat("Merging files from month", x, "out of", nrow(sent), "...\n")

  # Get the window/extent of the current month
  ext <- suppressWarnings(vect(sent$Window[[x]]))

  # Get the filenames of the current month and prepare filenames for the cleaned
  # and predicted files
  files <- tibble(
        Filename = sent$data[[x]]$filepath
      , Year     = sent$Year[x]
      , Month    = sent$Month[x]
    ) %>% mutate(
        FileNo        = sprintf("%02d", 1:n())
      , Month         = sprintf("%02d", Month)
      , FilenameFinal = paste0(tempi, "/", Year, "_", Month, "_", FileNo, ".tif")
  )

  # Subset to the files that do not exist yet
  exists <- file.exists(files$FilenameFinal)
  files <- subset(files, !exists)

  # Loop through the files and run cleaning and prediction
  pb <- txtProgressBar(min = 0, max = nrow(files), style = 3)
  predictions <- lapply(1:nrow(files), function(y) {

    # Prepare the file
    file_in  <- files$Filename[y]
    file_out <- files$FilenameFinal[y]
    pred <- prepareFile(file_in, file_out, ext)

    # Print progress
    setTxtProgressBar(pb, y)

    # Return the filename
    return(pred)
  })
  return(predictions)
})
