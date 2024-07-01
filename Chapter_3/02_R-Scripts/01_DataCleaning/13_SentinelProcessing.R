################################################################################
#### Sentinel Data Processing L1C to L2A
################################################################################
# Description: Processing the Sentinel 2 L1C data to achieve L2A products.

# Clear R's brain
rm(list = ls())

# Load required packages
library(sen2r)         # To automate the correction of Sentinel 2 data
library(tidyverse)     # To wrangle data
library(lubridate)     # To handle dates
library(terra)         # To handle spatial data
library(pbmcapply)     # To run stuff in parallel

# Need to ensure the following versions are installed!
if (packageVersion("sen2r") != "1.5.1") {
  stop("Please install sen2r version 1.5.1")
}
# devtools::install_version("sen2r", version = "1.5.1", repos = "http://cran.us.r-project.org")

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

# Specify the directories to the sentinel files
dir_l1c <- "/media/david/Elements/L1C"
dir_l2a <- "/media/david/Elements/L2A"

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

# Go through the groups, move the files to the computer and run the correction
lapply(group, function(x) {

  # Copy the files
  cat("Computing BOA reflectance from TOA images...\n")
  files_sub <- files[group == x]
  files_new <- file.path(tempdir(), basename(files_sub))
  file.copy(files_sub, tempdir(), recursive = T)

  # Run the correction
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
#### Getting Metadata & Quality Control
################################################################################
# Load paths to all downloaded and processed files
files <- list.dirs(path = dir_l2a, full.names = T, recursive = F)

# Check files for validity
nameval <- pbmclapply(files, safe_isvalid
  , ignore.interactive = T
  , mc.cores = detectCores() - 1
) %>% do.call(c, .) %>% unname()

# Some of the files fail when trying to convert them into a vrt. Let's test and
# see for which ones it works
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
convertable <- lapply(1:length(files), function(x) {
  suppressMessages({
    result <- tryCatch({
      s2_translate(files[x], outdir = tempdir())
    }, warning = function(w){return(w)}, error = function(e){return(e)})
  })
  result <- ifelse(length(result) == 0, "invalid", result)
  setTxtProgressBar(pb, x)
  if (grepl(result, pattern = ".vrt")) {
    return(T)
  } else {
    return(F)
  }
}) %>% do.call(c, .)

# Put all into a nice table and identify invalid products
valid   <- tibble(files, convertable, nameval)
invalid <- subset(valid, !convertable | !nameval)

# Remove the invalid folders
cat(paste0("\n", nrow(invalid), " folders are invalid and removed\n"))
unlink(invalid$files, recursive = T, force = T)

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
# Hence, I'll store the results to file so that we can reload them and only need
# to append metadata if necessary
if (file.exists("03_Data/03_Results/SentinelMetadata.rds")) {
    meta <- read_rds("03_Data/03_Results/SentinelMetadata.rds")
    toadd <- !(files %in% meta$filepath)
    if (length(files[toadd]) > 0) {
      cat("Collecting metadata for Sentinel tiles...\n")
      toadd <- getMeta(files[toadd])
      meta <- rbind(meta, toadd)
    }
  } else {
    meta <- getMeta(files)
}

# Store the metadata to file
write_rds(meta, "03_Data/03_Results/SentinelMetadata.rds")
meta <- read_rds("03_Data/03_Results/SentinelMetadata.rds")

# Ensure that all files are valid
table(meta$validname)

# In previous runs I have downloaded tiles that are not strictly needed here
# (i.e. they are not in the "todownload" file anymore). However, I don't want to
# delete them yet, in case we can use them at a later stage. However, I still
# want to remove their entries from the metadata just for now in order to avoid
# that we will further process all of them. Thus, let's only keep metadata of
# the files that we actually needed to download according to the "todownload"
# file.
todo <- "/media/david/Elements/Todownload.rds" %>%
  read_rds() %>%
  pull(name) %>%
  correctedName() %>%
  substr(start = 1, stop = nchar(.) - 21) %>%
  unique()

# Let's subset to those entries in the metadata
meta <- subset(meta, substr(name, start = 1, stop = nchar(name) - 21) %in% todo)

# Nest the data by year and month so that we can join it with the moving windows
# for each month
meta <- meta %>%
  group_by(Year, Month) %>%
  nest() %>%
  mutate(NumberFiles = map_int(data, nrow))

# Load moving windows that we generated based on the GPS data. We will later use
# these to crop the sentinel layers
sent <- "03_Data/03_Results/SentinelResults.rds" %>%
  read_rds() %>%
  select(-c(data, From, To, Files, FilesMetaData))

# Join the metadata with all moving windows. This will allow us to crop each
# tile to the moving window that we want to retain for that month
sent <- left_join(sent, meta, by = c("Year", "Month"))

# Store this to file
write_rds(sent, "03_Data/03_Results/SentinelJoined.rds")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/13_SentinelProcessing.rds")
cat("Done :)\n")
