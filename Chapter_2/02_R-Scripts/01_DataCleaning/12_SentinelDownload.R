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

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Note output directories
outdir_l1c <- "/media/david/Elements/L1C"
outdir_l2a <- "/media/david/Elements/L2A"

# Create directory into which we will download the data
dir.create(outdir_l1c, showWarnings = F)
dir.create(outdir_l2a, showWarnings = F)

# Make sure the directories exist
dir.exists(outdir_l1c)
dir.exists(outdir_l2a)

# Login to scihub
write_scihub_login("dodx9", "Scihubbuster69_")

# Read the file containing all data that we need to download
todownload <- read_rds("/media/david/Elements/Todownload.rds")

# Depending on the availability of L1C or L2A products, we put the downloaded
# data into two different directories
todownload$outdir <- ifelse(todownload$level == "1C", outdir_l1c, outdir_l2a)

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

# We only need to download files that haven't been downloaded yet. Thus, let's
# check which files are already present
exists <-
  dir.exists(paste0(todownload$outdir, "/", todownload$name)) |
  dir.exists(paste0(todownload$outdir, "/", todownload$name))
todownload <- subset(todownload, !exists)

# Order jobs so that the newest dates are coming first -> not yet archvied in
# the long-term archive
todownload <- arrange(todownload, desc(sensing_datetime))

# Nest the data by its group
todownload <- todownload %>%
  group_by(GroupID, Username, Password) %>%
  nest()

# Remove the first couple of entries as there appears to be an issue
todownload <- todownload[-(1:14), ]

# Check out what we need to download
print(todownload)

# Specify the number of cores to use
cores <- 4

# Run through the groups and download data
for (i in 1:nrow(todownload)) {

  # Extract the files we want to download in that group
  getit <- todownload$data[[i]]
  getitnext <- todownload$data[[i + 1]]

  # Check if the files are already available
  cat("Checking if files are online...\n")
  write_scihub_login(todownload$Username[i], todownload$Password[i])
  online <- safe_is_online(getit$todownload)
  write_scihub_login(todownload$Username[i + 1], todownload$Password[i + 1])
  onlinenext <- safe_is_online(getitnext$todownload)

  if (!all(online)) {
    ordered <- tryCatch({
      s2_order(getit$todownload[!online])
    }, error = function(e) {return(NA)})

    # Wait until some products are online
    letuswait <- T
    while (letuswait) {
      online <- safe_is_online(getit$todownload)
      letuswait <- ifelse(!any(online), T, F)
      if (letuswait) {
        cat("No  products available. Waiting for 5 Minutes...\n")
        Sys.sleep(5 * 60)
      }
    }
  }

  # Order files from the next group if necessary
  if (!all(onlinenext)) {
    tryCatch({
      s2_order(getitnext$todownload[!onlinenext])
    }, error = function(e) {return(NA)})
  }

  # Login with correct account
  write_scihub_login(todownload$Username[i], todownload$Password[i])

  # Remove the products that are not available. Note that they will not be
  # downloaded in this or any of the next iterations. The code will have to be
  # rerun multiple times to ensure that these files are eventually downloaded as
  # well. The reason I chose this approach is because some fails always fail to
  # be dwonloaded, yet I want the code to continue running nevertheless.
  getit <- getit[online, ]

  # If there is more than one product for download, we want to download two
  # images at a time
  cat("Downloading Sentinel 2 imagery...\n")
  if (nrow(getit) > 1) {

      # Split data into two groups
      getit$GroupID2 <- rep(1:cores, length.out = nrow(getit))

      # Download data in parallel
      pbmclapply(1:cores, mc.cores = cores, function(j) {
        getit_sub <- subset(getit, GroupID2 == j)
        for (k in 1:nrow(getit_sub)) {
          suppressWarnings(
            tryCatch({
              s2_download(getit_sub$todownload[k], outdir = getit$outdir[k])
            }, error = function(e) {return(e)})
          )
        }
      })

    # If there is only one product for download, get it
    } else {
      success <- tryCatch({
        s2_download(getit$todownload[1], outdir = getit$outdir[1])
      }, error = function(e) {return(e)})
  }
}
