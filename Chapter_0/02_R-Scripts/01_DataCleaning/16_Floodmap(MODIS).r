################################################################################
#### Download and Classify MODIS MCD43A4 Imagery
################################################################################
# Description: This script allows to download, process and classify MODIS
# imagery automatically.

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)
library(terra)
library(tidyverse)
library(lubridate)
library(splitstackshape)
library(floodmapr)

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Are you downloading data for validation?
validation <- F

################################################################################
#### Specify Dates
################################################################################
# Specify dates for which you want inundation maps
dates <- c(
    "2012-01-09"
  # , "2012-01-17"
  # , "2012-01-25"
  # , "2012-02-02"
  # , "2012-02-10"
  # , "2015-11-17"
  # , "2015-11-25"
  # , "2015-12-03"
  # , "2015-12-11"
  # , "2015-12-19"
  # , "2015-12-27"
  # , "2016-01-01"
  # , "2016-01-09"
  # , "2016-01-17"
  # , "2016-11-08"
  # , "2016-11-16"
  # , "2016-11-24"
  # , "2016-12-02"
  # , "2017-03-14"
  # , "2017-03-22"
  # , "2017-11-01"
  # , "2017-11-09"
  # , "2017-11-17"
  # , "2017-11-25"
  # , "2017-12-03"
  # , "2017-12-11"
  # , "2017-12-19"
  # , "2017-12-27"
  # , "2018-01-01"
  # , "2018-01-09"
  # , "2018-01-17"
  # , "2018-03-06"
  # , "2018-03-14"
  # , "2018-03-22"
  # , "2018-03-30"
  # , "2018-04-07"
  # , "2018-04-15"
  # , "2018-04-23"
  # , "2018-05-01"
  # , "2018-05-09"
  # , "2018-05-17"
  # , "2018-05-25"
  # , "2018-08-29"
  # , "2018-09-06"
  # , "2018-09-14"
  # , "2018-09-22"
  # , "2018-09-30"
  # , "2018-10-08"
  # , "2018-10-16"
  # , "2018-10-24"
  # , "2018-11-01"
  # , "2018-11-09"
  # , "2018-11-17"
  # , "2018-11-25"
  # , "2018-12-03"
  # , "2018-12-11"
  # , "2018-12-19"
  # , "2019-01-17"
  # , "2019-01-25"
  # , "2019-02-18"
  # , "2019-02-26"
  # , "2019-03-06"
  # , "2019-03-14"
  # , "2019-03-22"
  # , "2019-03-30"
  # , "2019-04-07"
  # , "2019-04-15"
  # , "2019-04-23"
  # , "2019-05-01"
  # , "2019-05-17"
  # , "2019-05-25"
  # , "2019-06-02"
  # , "2019-06-10"
  # , "2019-06-18"
  # , "2019-06-26"
  # , "2019-07-04"
  # , "2019-07-12"
  # , "2019-07-20"
  # , "2019-07-28"
  # , "2019-08-05"
  # , "2019-08-13"
  # , "2019-08-21"
  # , "2019-08-29"
  # , "2019-09-06"
  # , "2019-09-14"
  # , "2019-09-22"
)

################################################################################
#### For Validation Purposes Only
################################################################################
# # This part contains some code that allows to compare the results of our own
# # classification algorithm against the results of the original ORI algorithm
# set.seed(1234)
# dates <- dir(
#     path        = "03_Data/01_RawData/ORI/01_Maps"
#   , pattern     = ".tif$"
# ) %>% substr(., start = 1, stop = 10) %>% as.Date(format = "%Y-%m-%d")
#
# # Extract the months
# dates <- data.frame(Date = dates, Month = as.factor(month(dates)))
#
# # Randomly draw 48 dates (4 x 12) and make sure that each month has equal
# # representation
# dates <- stratified(dates, "Month", 4)[[1]] %>% as.character()
#
# # Sort them
# dates <- sort(dates)
#
# # Check them
# dates
#
# # Download only dates that are not available yet
# avail <- dir(path = "/home/david/Schreibtisch/15. PhD/Chapter_0/03_Data/02_CleanData/00_Floodmaps/03_Validation")
# avail <- as.Date(substr(avail, start = 1, stop = 10), format = "%Y.%m.%d")
# dates <- dates[!as.Date(dates) %in% avail]

################################################################################
#### Download Data
################################################################################
# Load reference raster (floodmap reference raster)
r <- raster("03_Data/01_RawData/ORI/FloodmapReference.tif")

# Loop through the dates and download and classify floodmaps. Note: While it is
# possible to paralellize this process, I decided not to. The reason is that I
# will use dynamic floodmaps. Hence, later floodmaps will depend on earlier
# classifications. Therefore, I don't want to run classifications in paralell.
for (i in 1:length(dates)){

  # Download modis image
  downloaded <- modis_download(
      dates     = dates[i]
    , outdir    = "03_Data/01_RawData/MODIS/MCD43A4/01_Maps"
    , tmpdir    = tempdir()
    , username  = "DoDx9"
    , password  = "EarthData99"
    , overwrite = T
  )

  # Print update
  cat("File downloaded and extracted. Proceeding with classification...\n")

  # Load it
  loaded <- modis_load(downloaded)

  # Identify filenames of all already classified watermaps
  filenames <- dir(
      path        = "03_Data/02_CleanData/00_Floodmaps/01_Original"
    , pattern     = ".tif$"
    , full.names  = T
  )

  # Identify dates of all already classified watermaps
  filedates <- filenames %>%
    basename() %>%
    substr(start = 1, stop = 10) %>%
    as.Date(format = "%Y.%m.%d")

  # Create a dynamic watermask
  watermask <- modis_watermask(
      date      = as.Date(dates[i])
    , filenames = filenames
    , filedates = filedates
    , years     = 5
    , threshold = 0.99
  )

  # Check for bimodality in the modis image
  is_bimodal <- modis_bimodal(
      x         = loaded
    , watermask = watermask
    , drymask   = NULL
  )

  # In case the file is not bimodal, we can't classify it and need to skip
  if (!is_bimodal){
    cat("Image is not bimodal. Can't classify watermap. Skipping to next image\n")
    next
  }

  # Classify modis image. We provide a dynamic watermask here
  classified <- modis_classify(
      x                 = loaded
    , watermask         = watermask
    , drymask           = NULL
    , ignore.bimodality = T
  )

  # Print update
  cat("Image classified. Resampling and Storing now...\n")

  # We also want to resample the classified image, so that it matches the
  # reference raster. We'll use the nearest neighbor method here
  classified <- resample(
      x       = classified
    , y       = r
    , method  = "ngb"
  )

  # We need to make sure that the numbers align with the ORI classification.
  # That is, we want dryland = 255, water = 0, clouds = 127
  rcl <- data.frame(old = c(0, 1, NA), new = c(255, 0, 127))
  reclassified <- reclassify(classified, rcl)

  # Specify final filename
  name <- gsub(dates[i], pattern = "-", replacement = ".")
  name <- paste0("03_Data/02_CleanData/00_Floodmaps/01_Original/", name, ".tif")

  # In case we are validating
  if (validation){
    dir.create("03_Data/02_CleanData/00_Floodmaps/03_Validation")
    name <- gsub(dates[i], pattern = "-", replacement = ".")
    name <- paste0("03_Data/02_CleanData/00_Floodmaps/03_Validation/", name, ".tif")
  }

  # For whatever reason the files are much smaller when stored using "raster",
  # rather than using "terra".
  # reclassified <- raster(reclassified)

  # Store raster
  writeRaster(reclassified, name, overwrite = TRUE)

  # Print update
  cat("Image", i, "out of", length(dates), "finished...\n")
}
