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

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
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
write_scihub_login("dodx9", "Scihubbuster69_")

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
todownload$Password <- "Scihubbuster69_"

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
exists <- exists(todownload$name)
todownload <- subset(todownload, !exists)

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

# ################################################################################
# #### LEGACY
# ################################################################################
# # Order jobs so that the newest dates are coming first -> not yet archvied in
# # the long-term archive
# # todownload <- arrange(todownload, desc(sensing_datetime))
#
# # Nest the data by its group
# todownload <- todownload %>%
#   # group_by(GroupID, Username, Password) %>%
#   group_by(Username, Password) %>%
#   nest()
#
# # Remove the first couple of entries as there appears to be an issue
# # todownload <- todownload[-(1:4), ]
#
# # Check out what we need to download
# print(todownload)
#
# # Specify the number of cores to use
# cores <- 4
#
# # Run through the groups and download data
# for (i in 1:nrow(todownload)) {
#
#   # Extract the files we want to download in that group
#   getit <- todownload$data[[i]]
#   if (nrow(todownload) > 1 & i != nrow(todownload)) {
#     getitnext <- todownload$data[[i + 1]]
#   }
#
#   # Check if the files are already available
#   cat("Checking if files are online...\n")
#   write_scihub_login(todownload$Username[i], todownload$Password[i])
#   online <- safe_is_online(getit$todownload)
#   if (nrow(todownload) > 1 & i != nrow(todownload)) {
#     write_scihub_login(todownload$Username[i + 1], todownload$Password[i + 1])
#     onlinenext <- safe_is_online(getitnext$todownload)
#   }
#
#   if (!all(online)) {
#     ordered <- tryCatch({
#       s2_order(getit$todownload[!online])
#     }, error = function(e) {return(NA)})
#
#     # Wait until some products are online
#     letuswait <- T
#     while (letuswait) {
#       online <- safe_is_online(getit$todownload)
#       letuswait <- ifelse(!any(online), T, F)
#       if (letuswait) {
#         cat("No  products available. Waiting for 5 Minutes...\n")
#         Sys.sleep(5 * 60)
#       }
#     }
#   }
#
#   # Order files from the next group if necessary
#   if (nrow(todownload) > 1 & i != nrow(todownload)) {
#     if (!all(onlinenext)) {
#       tryCatch({
#         s2_order(getitnext$todownload[!onlinenext])
#       }, error = function(e) {return(NA)})
#     }
#   }
#
#   # Login with correct account
#   write_scihub_login(todownload$Username[i], todownload$Password[i])
#
#   # Remove the products that are not available. Note that they will not be
#   # downloaded in this or any of the next iterations. The code will have to be
#   # rerun multiple times to ensure that these files are eventually downloaded as
#   # well. The reason I chose this approach is because some fails always fail to
#   # be dwonloaded, yet I want the code to continue running nevertheless.
#   getit <- getit[online, ]
#
#   # If there is more than one product for download, we want to download two
#   # images at a time
#   cat("Downloading Sentinel 2 imagery...\n")
#   if (nrow(getit) > 1) {
#
#       # Split data into two groups
#       getit$GroupID2 <- rep(1:cores, length.out = nrow(getit))
#
#       # Download data in parallel
#       pbmclapply(1:cores, mc.cores = cores, function(j) {
#         getit_sub <- subset(getit, GroupID2 == j)
#         for (k in 1:nrow(getit_sub)) {
#           suppressWarnings(
#             tryCatch({
#               k <- 3
#               s2_download(getit_sub$todownload[k], outdir = getit_sub$outdir[k])
#             }, error = function(e) {return(e)})
#           )
#         }
#       })
#
#     # If there is only one product for download, get it
#     } else {
#       success <- tryCatch({
#         s2_download(getit$todownload[1], outdir = getit$outdir[1])
#       }, error = function(e) {return(e)})
#   }
# }
