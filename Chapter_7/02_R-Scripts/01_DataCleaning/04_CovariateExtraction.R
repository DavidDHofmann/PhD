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
setwd(wd)

# Load required packages
library(terra)       # To handle spatial data
library(tidyverse)   # For data wrangling
library(pbmcapply)   # For multicore abilities with progress bar
library(lubridate)   # To handle timestamps
library(runner)      # To have a moving window

# Load activity data
act <- read_csv("03_Data/02_CleanData/ActivityData.csv")

# Round timestamps to the nearest hour
act$TimestampRounded <- round_date(act$Timestamp, unit = "hour")

# The three layers from which we will extract
layers <- c("Precipitation", "Temperature", "CloudCover")

################################################################################
#### Covariate Extraction of Climate Data
################################################################################
# Loop through the different rasterlayers (actually stacks) and extract
# covariate values
for (i in layers) {

  # Load covariate layer
  cat("Extracting", i,"values...\n")
  covariate <- rast(paste0("03_Data/02_CleanData/", i, ".tif"))

  # Generate dates from layernames
  if (i == "CloudCover") {
      covariate_dates <- covariate %>%
        names() %>%
        substr(start = 1, stop = 10) %>%
        ymd()
    } else {
      covariate_dates <- covariate %>%
        names() %>%
        substr(start = 1, stop = 21) %>%
        ymd_h()
  }

  # Stop if there are any NAs
  if (any(is.na(covariate_dates))) {
    stop("At least one of the dates failed to parse\n")
  }

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
    , mc.cores           = (detectCores() - 1) / 2
    , FUN                = function(x) {
      if (i == "CloudCover") {
        index <- which.min(abs(as_date(x) - covariate_dates))[1]
      } else {
        index <- which.min(abs(x - covariate_dates))[1]
      }
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

# Ensure that there are no NAs
if (sum(is.na(act[, layers])) > 0) {
  stop("Some of the extracted values are NA...\n")
}

################################################################################
#### Moon Angle
################################################################################
# Load moonangle data
moon <- read_csv("03_Data/02_CleanData/MoonAngle.csv")
moon <- dplyr::select(moon, timestamp, moonPhase, moonAltDegrees)

# Join the data
act <- left_join(act, moon, by = c("Timestamp" = "timestamp"))

# Ensure that there are no NAs
if (sum(is.na(act$MoonAngle) | is.na(act$moonPhase)) > 0) {
  stop("Some of the extracted values are NA...\n")
}

################################################################################
#### Moonlight Statistics
################################################################################
# Load moonlight data
moon <- read_csv("03_Data/02_CleanData/Moonlight.csv")

# Remove undesired columns
moon <- dplyr::select(moon, -c(lon, lat))

# Compute the delay since sunset to maximum moonlight
moon$maxMoonDelay <- difftime(moon$maxMoonTime, moon$sunset, units = "mins")

# Specify which date we want to match. Before 12:00, we match the day before,
# after 12, the day after
act <- mutate(act, DateToMatch = if_else(hour(Timestamp) >= 12
  , Date
  , Date - days(1))
)

# Join the nightly statistics respectively
act <- left_join(act, moon, by = c("DateToMatch" = "date"))

################################################################################
#### Activity Metrics
################################################################################
# Combine the x and y axis into a single activity metric
act$Act <- act$ActX + act$ActY

# Function to apply a moving window to the data
windowActivity <- function(x, timestamps, window_size = "60_mins") {
  x_coarse <- runner(x
    , k   = window_size
    , idx = timestamps
    , f   = function(x) {mean(x)}
  )
  return(x_coarse)
}

# Loop through the individuals and compute activity within moving windows
act_nested <- act %>% nest(Data = -DogID)
act_nested$Data <- pbmclapply(
    X                  = act_nested$Data
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(x) {
  x <- arrange(x, Timestamp)
  x$Act15 <- windowActivity(x$Act, x$Timestamp, window_size = "15 mins")
  x$Act30 <- windowActivity(x$Act, x$Timestamp, window_size = "30 mins")
  x$Act45 <- windowActivity(x$Act, x$Timestamp, window_size = "45 mins")
  x$Act60 <- windowActivity(x$Act, x$Timestamp, window_size = "60 mins")
  return(x)
})

# Unnest again
act <- unnest(act_nested, Data)

################################################################################
#### Reorder and Store
################################################################################
# Reorder the columns slightly
act <- dplyr::select(act, c(
    DogID
  , CollarID
  , Act
  , Act15
  , Act30
  , Act45
  , Act60
  , Date
  , TimestampRounded
  , Timestamp
  , GPSTimestamp
  , x
  , y
  , DOP
  , State
  , SunAngle
  , MoonAngle   = moonAltDegrees
  , MoonPhase   = moonPhase
  , SunsetLast  = sunset_last
  , Sunrise     = sunrise
  , Sunset      = sunset
  , SunriseNext = sunrise_next
  , ToD
  , ToD2
  , Cycle
  , minMoonlightIntensity
  , meanMoonlightIntensity
  , maxMoonlightIntensity
  , minMoonPhase
  , meanMoonPhase
  , maxMoonPhase
  , maxMoonTime
  , maxMoonDelay
  , Temperature
  , Precipitation
  , CloudCover
))

# Function to distinguish seasons
getSeason <- function(x) {
    DS <- as.Date("2020-04-01", format = "%Y-%m-%d") # Start Dry
    WS <- as.Date("2020-11-01", format = "%Y-%m-%d") # Start Wet
    d <- as.Date(strftime(x, format = "2020-%m-%d"))
    season <- ifelse(d >= DS & d < WS, "Dry", "Wet")
    return(season)
}

# Apply it to separate the data into the seasons
act$Season <- getSeason(act$Date)

# Generate a binary variable indicating if it rained or not
act$Rain <- act$Precipitation > 0

# Store the final data to file
write_csv(act, "03_Data/02_CleanData/ActivityDataCovariates.csv")
