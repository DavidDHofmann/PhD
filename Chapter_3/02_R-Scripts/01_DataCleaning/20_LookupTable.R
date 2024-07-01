################################################################################
#### Lookup Table
################################################################################
# Description: Prepare a lookup table that we can use to generate the landscape
# for a desired date

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_3")

# Load required packages
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(hms)            # To handle tims
library(terra)          # To handle spatial data
library(pbmcapply)      # For multicore ability with progress bar

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Lookup for Model Fitting and Validation
################################################################################
# Load covariates and dispersers
covs <- read_rds("03_Data/02_CleanData/Covariates.rds")
disp <- read_csv("03_Data/02_CleanData/Dispersers.csv") %>% subset(State == "Disperser")

# Span a range of dates for which we might want to extract data
disp_dates <- tibble(
  Timestamp = seq(
        from = update(range(disp$Timestamp)[1], hour = 3, minutes = 0, seconds = 0)
      , to   = update(range(disp$Timestamp)[2], hour = 23, minutes = 0, seconds = 0)
      , by   = "4 hours"
    )
  ) %>% subset(hour(Timestamp) %in% c(3, 7, 15, 19, 23))
sim_dates <- tibble(
  Timestamp = seq(
        from = ymd_hms("2000-01-01 03:00:00")
      , to   = ymd_hms("2000-12-31 23:00:00")
      , by   = "4 hours"
    )
  ) %>% subset(hour(Timestamp) %in% c(3, 7, 15, 19, 23))

# Put the dates together
all_dates <- rbind(disp_dates, sim_dates)

# Prepare a lookup table that gives the index of the layer that is closest in
# date
lookup <- expand_grid(
    Timestamp = all_dates$Timestamp
  , covs[, c("Type", "Covariate")]
)

# Join with the covariate dates
lookup <- covs %>%
  left_join(lookup, ., by = c("Type", "Covariate"))

# Go through the timestamps for which we will later extract covariates (either
# for model fitting or simulation / prediction)
cat("Preparing lookup table...\n")
close <- pbmclapply(
    X                  = seq_len(nrow(lookup))
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(i) {

    # Extract loop information
    timestamp    <- lookup$Timestamp[i]
    layerdates   <- lookup$Dates[[i]]$Layerdate

    # If we're looking at pan data, we'll not extract at the closest date but
    # from the month within which a step falls (if possible)
    if (lookup$Covariate[i] == "DistanceToPans" & lookup$Type[i] == "Dynamic") {
      index <- which(
        paste0(year(layerdates), "-", month(layerdates)) ==
        paste0(year(timestamp), "-", month(timestamp))
      )
      if (length(index) == 0) {
        index <- which.min(abs(difftime(timestamp, layerdates, units = "hours")))
      }
    } else {
      index <- which.min(abs(difftime(timestamp, layerdates, units = "hours")))
    }
    closest <- lookup$Dates[[i]][index, ]
    return(closest)
}) %>% bind_rows()

# Bind this with the lookup table
lookup <- lookup %>%
  dplyr::select(-c(Dates, Filename)) %>%
  cbind(., close)

# Store the table to file
write_rds(lookup, "03_Data/02_CleanData/LookupTable.rds")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/20_LookupTable.rds")
cat("Done :)\n")
