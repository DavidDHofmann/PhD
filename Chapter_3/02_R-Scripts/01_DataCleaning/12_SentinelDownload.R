################################################################################
#### Sentinel Data Download
################################################################################
# Description: Download of the files specified in the previous script

# Clear R's brain
rm(list = ls())

# Load required packages
library(sen2r)         # To automate the download of Sentinel 2 data
library(tidyverse)     # For data wrangling
library(pbmcapply)     # For multicore abilities with progress bar
library(sf)            # To parse wkt
library(lubridate)     # To handle dates

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

# Note output directories
outdir_l1c <- "/media/david/Elements/L1C"
outdir_l2a <- "/media/david/Elements/L2A"

# Create directory into which we will download the data
dir.create(outdir_l1c, showWarnings = F)
dir.create(outdir_l2a, showWarnings = F)

# Login to scihub
load("/home/david/ownCloud/Dokumente/Bibliothek/Wissen/R-Scripts/ScihubLogin.rds")
# username <- "username"
# password <- "password"
write_scihub_login(username, password)

# Read the file containing all data that we need to download
todownload <- read_rds("/media/david/Elements/Todownload.rds")

# The output directory of the file depends on whether we can download the l2a
# product or if we need to go with the l1c product
todownload$outdir <- ifelse(todownload$level == "1C", outdir_l1c, outdir_l2a)

# Let's generate the L2A product names, so that we also know the final names of
# the converted L1C products
todownload$L2A_name <- correctedName(todownload$name)

# Give each job a unique ID
todownload$JobID <- 1:nrow(todownload)

# Instead of requesting all downloads at once, I'd like to go through the files
# in groups. Thus, let's create groups with around 20 members in each group
n <- 20
ngroups <- round(nrow(todownload) / n, -1)
todownload$GroupID <- sort(rep(1:ngroups, length.out = nrow(todownload)))
print(todownload)

# Depending on the group we will use different login credentials
todownload$Username <- ifelse(todownload$GroupID %% 2 == 1, "dodx9", "dodx92")
todownload$Password <- password

# Let's write a function that we can use to check if a file already exists
# The function needs to tolerate minor mismatches in names as the ingestion date
# between the esa and google servers not always matches 100%
exists <- function(x) {

  # List all files that are already downloaded
  present_l1c <- list.dirs(outdir_l1c, recursive = F)
  present_l2a <- list.dirs(outdir_l2a, recursive = F)

  # Keep only base names
  present_l1c <- basename(present_l1c)
  present_l2a <- basename(present_l2a)

  # Remove ending characters
  present_l1c_short <- substr(present_l1c, start = 1, stop = nchar(present_l1c) - 21)
  present_l2a_short <- substr(present_l2a, start = 1, stop = nchar(present_l2a) - 21)

  # Prepare l2a names for all files
  y <- correctedName(x)

  # We remove the last few characters from the filenames as they are irrelevant
  x_short <- substr(x, start = 1, stop = nchar(x) - 21)
  y_short <- substr(y, start = 1, stop = nchar(y) - 21)

  # Check if the file "x" is already present, either through its l1c or l2a name
  exists_l1c <- x_short %in% present_l1c_short
  exists_l2a <- y_short %in% present_l2a_short
  ex <- exists_l1c | exists_l2a

  # Return the result
  return(ex)
}

# Check which files already exist
exists     <- exists(todownload$name)
todownload <- subset(todownload, !exists)
print(todownload)

# For some files, the download through sentinel hub doesn't work properly. For
# these, I'll use gcloud
if (nrow(todownload) > 0) {
  ext <- st_as_sf(todownload, wkt = "footprint")
  ext <- st_set_crs(ext, 4326)
  cat("Downloading from Google Cloud...\n")
  for (i in 1:nrow(ext)) {
    gcfile <- s2_list(
        # spatial_extent = ext$footprint[i]
        tile           = ext$id_tile[i]
      , orbit          = ext$id_orbit[i]
      , time_interval  = c(
            as.Date(ymd_hms(ext$sensing_datetime[i])) - days(1)
          , as.Date(ymd_hms(ext$sensing_datetime[i])) + days(1)
        )
      , server         = "gcloud"
    )
    exists <- dir.exists(file.path(ext$outdir[i], names(gcfile)))
    if (!exists) {
      s2_download(gcfile, outdir = ext$outdir[i], overwrite = T)
    }
  }
}

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/12_SentinelDownload.rds")
cat("Done :)\n")
