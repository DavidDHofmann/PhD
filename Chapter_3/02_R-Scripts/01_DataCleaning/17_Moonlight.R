################################################################################
#### Moonlight Statistics
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)    # For data-wrangling
library(lubridate)    # To handle dates
library(pbmcapply)    # For multicore abilities
library(suncalc)      # To identify sunrise and sunset times
library(hms)          # To work with times

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

# Define coordinates we want to get astronomical data for
lon <- 23.5
lat <- -19

# Span a vector of timestamps for which we want to identify moonlight and sun
# metrics
timestamps1 <- seq(
    from = ymd_hms("1999-12-31 23:00:00")
  , to   = ymd_hms("2001-01-01 03:00:00")
  , by   = "5 mins"
)
timestamps2 <- seq(
    from = ymd_hms("2011-01-01 03:00:00")
  , to   = ymd_hms("2023-12-31 23:00:00")
  , by   = "5 mins"
)
timestamps <- c(timestamps1, timestamps2)

# If you want to rerun
# file.remove("03_Data/02_CleanData/MoonlightDetailed.rds")
# file.remove("03_Data/02_CleanData/Moonlight.rds")

# Only run if the file has not been created yet
if (!file.exists("03_Data/02_CleanData/Moonlight.rds")) {

  # Compute moonlight statistics for all datapoints
  moon <- moonlightIntensity(
      date = timestamps
    , lon  = 23.5
    , lat  = -19
    , e    = 0.21
  )

  # Group by the hours we are interested in
  moon <- mutate(moon, date_new = case_when(
      hour(date) %in% c(3, 4, 5, 6) ~ update(date, hour = 3, minutes = 0, seconds = 0)
    , hour(date) %in% c(7, 8, 9, 10, 11, 12, 13, 14) ~ update(date, hour = 7, minutes = 0, seconds = 0)
    , hour(date) %in% c(15, 16, 17, 18) ~ update(date, hour = 15, minutes = 0, seconds = 0)
    , hour(date) %in% c(19, 20, 21, 22) ~ update(date, hour = 19, minutes = 0, seconds = 0)
    , hour(date) %in% c(23) ~ update(date, hour = 23, minutes = 0, seconds = 0)
    , hour(date) %in% c(0, 1, 2) ~ update(date, hour = 23, minutes = 0, seconds = 0) - days(1)
    , .default = date
  ))

  # Summarize
  moon_aggr <- moon %>%
    group_by(date_new) %>%
    summarize(
        meanMoonPhase    = mean(moonPhase)
      , meanMoonAlt      = mean(moonAltDegrees)
      , meanSunAlt       = mean(sunAltDegrees)
      , meanMoonlight    = mean(moonlightModel)
      , meanMoonlightLux = mean(moonlightModel) * 0.32
      , Night            = sum(sunAltDegrees < -18) / length(sunAltDegrees)
    ) %>%
    rename(Timestamp = date_new)

  # Split into dark and bright periods
  moon_aggr <- mutate(moon_aggr, LightType = if_else(
      Night > 0.5 & meanMoonlightLux < 0.02, "Dark", "Bright"
  ))

  # Store to file
  write_rds(moon_aggr, "03_Data/02_CleanData/Moonlight.rds")

}

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/17_Moonlight.rds")
cat("Done :)\n")
