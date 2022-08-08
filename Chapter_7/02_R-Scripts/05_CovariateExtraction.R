################################################################################
#### Extraction of Covariates
################################################################################
# Description: Extracting climate covariates for the activity data. Note that
# I'll use the GPS coordinates assigned to each activity fix to extract
# precipitation and temperature values

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
# wd <- "C:/Users/david/switchdrive/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(terra)       # To handle spatial data
library(tidyverse)   # For data wrangling
library(pbmcapply)   # For multicore abilities with progress bar
library(lubridate)   # To handle timestamps

# Load activity data
act <- read_csv("03_Data/02_CleanData/ActivityDataMoonphase.csv")

# Round timestamps to the nearest hour
act$TimestampRounded <- round_date(act$Timestamp, unit = "hour")

# The two layers from which we will extract
layers <- c("Precipitation", "Temperature", "CloudCover")

################################################################################
#### Covariate Extraction
################################################################################
# Loop through the different rasterlayers (actually stacks) and extract
# covariate values
for (i in layers) {

  # Load covariate layer
  covariate <- rast(paste0("03_Data/02_CleanData/", i, ".tif"))

  # Generate dates from layernames
  covariate_dates <- covariate %>%
    names() %>%
    substr(start = 1, stop = 21) %>%
    ymd_h()

  # I want to avoid repeatedly extracting covariate values at the same pixel for
  # the same date. However, this would happen if we simply loop through all the
  # activity data and extract from the closes location (as most fixes fall into
  # the same pixel). Instead, I want to identify at which pixels we need to
  # extract, then we can concentrated on those only.
  covariate_id <- covariate[[1]]
  covariate_id[] <- 1:ncell(covariate_id)

  # Identify for which pixels we need to extract values
  act$CovariatePixelID <- terra::extract(covariate_id, cbind(act$x, act$y))[, 1]

  # You can see that the number of pixels we need to extract from is fairly
  # small
  pixels <- unique(act$CovariatePixelID)
  length(pixels)

  # For each unique timestamp, find the layer that is closest in time to it
  dates <- unique(act$TimestampRounded)
  closest <- pbmclapply(
      X                  = dates
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x) {
      index <- which.min(abs(x - covariate_dates))[1]
      close <- data.frame(
          TimestampRounded    = x
        , CovariateLayerIndex = index
      )
      return(close)
  }) %>% do.call(rbind, .)

  # Join those to the activity data
  act <- left_join(act, closest, by = "TimestampRounded")

  # Run the extraction at the pixels we are interested in
  extr <- terra::extract(x = covariate, y = pixels)
  extr <- cbind(Pixel = pixels, extr)
  names(extr) <- c("CovariatePixelID", 1:length(covariate_dates))
  extr <- pivot_longer(extr, 2:ncol(extr), names_to = "CovariateLayerIndex", values_to = i)
  extr$CovariateLayerIndex <- as.numeric(extr$CovariateLayerIndex)

  # Left join to the activity data
  act <- left_join(act, extr, by = c("CovariatePixelID", "CovariateLayerIndex"))

  # Remove undesired columns
  act$CovariatePixelID    <- NULL
  act$CovariateLayerIndex <- NULL

}

# Organize the columns more nicely
test <- dplyr::select(act, c(
    DogID
  , CollarID
  , ActX
  , Date
  , Year
  , Month
  , Time
  , Hour
  , TimestampRounded
  , Timestamp
  , Lag
  , GPSTimestamp
  , x
  , y
  , DOP
  , State
  , SunAngle
  , ToD
  , ToD2
  , Cycle
  , minMoonlightIntensity
  , meanMoonlightIntensity
  , maxMoonlightIntensity
  , minMoonPhase
  , meanMoonPhase
  , maxMoonPhase
  , Temperature
  , Precipitation
  , CloudCover
))

# Store the final data to file
write_csv("03_Data/02_CleanData/ActivityDataCovariates.csv")
