################################################################################
#### Cleaning Collar Data
################################################################################
# Cleaning GPS and Activity data downloaded from dropbox

# Clear R's brain
rm(list = ls())

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_7")

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(pbmcapply)   # To allow multicore progress bars
library(suncalc)     # To compute the sun angle

################################################################################
#### GPS Data
################################################################################
# Identify all files containing GPS data
files <- dir(path = "03_Data/01_RawData/POPECOL/GPS", full.names = T)

# Loop through the GPS files and clean them
cat("Loading and cleaning GPS data...\n")
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
gps_dat <- lapply(1:length(files), function(i) {

  # Extract DogID from filename so that we can put it into the dataframe
  dog <- unlist(strsplit(basename(files[i]), split = "_"))[2]

  # Identify column separator
  sep <- readLines(files[i], n = 1)
  sep <- if_else(grepl(x = sep, pattern = ";", useBytes = T), ";", ",")

  # Load the file as plain text and remove funny characters
  gps_dat <- read_file(files[i], local = locale(encoding = "latin1"))
  if (sep == ",") {
      gps_dat <- gsub(gps_dat, pattern = ",", replacement = ";")
    } else {
      gps_dat <- gsub(gps_dat, pattern = ",", replacement = ".")
  }
  gps_dat <- gsub(gps_dat, pattern = "°", replacement = "")
  gps_dat <- gsub(gps_dat, pattern = "/", replacement = ".")
  gps_dat <- gsub(gps_dat, pattern = " \\[(.*?)\\]", replacement = "")

  # Load gps_dat into a data frame and do some more cleaning
  gps_dat <- suppressMessages(read_delim(gps_dat
      , local     = locale(encoding = "latin1")
      , delim     = ";"
      , col_types = cols()
      , col_names = T
    )) %>%

    # Add Dog Name
    mutate(DogID = dog) %>%

    # Retrieve timestamp
    mutate(Timestamp = as.POSIXct(
      paste(UTC_Date, UTC_Time), tz = "UTC", format = "%d.%m.%Y %H:%M:%S")
    ) %>%

    # Remove special characters like [°]
    setNames(gsub(names(.), pattern = " \\[*", replacement = "")) %>%
    setNames(gsub(names(.), pattern = "\\]", replacement = "")) %>%

    # Keep only desired columns
    dplyr::select(.
      , DogID     = DogID
      , CollarID  = CollarID
      , x         = Longitude
      , y         = Latitude
      , Timestamp = Timestamp
      , DOP       = DOP
    )

  # Return the filepath
  setTxtProgressBar(pb, i)
  return(gps_dat)

}) %>% do.call(rbind, .)

# Let's load the collar handling dates
handling <- read_csv("03_Data/01_RawData/POPECOL/CollarHandling.csv")

# Loop through the collar handling dates and only keep valid data, i.e. data
# that falls within the range of dates after the first but before the last date
cat("Subsetting to valid data based on collar handling dates...\n")
pb <- txtProgressBar(min = 0, max = nrow(handling), style = 3)
gps_dat <- lapply(1:nrow(handling), function(i) {
  start <- handling$StartDate_UTC[i]
  stop <- handling$EndDate_UTC[i]
  stop <- if_else(is.na(stop), ymd_hms(Sys.time(), tz = "UTC"), stop)
  sub <- subset(gps_dat
    , Timestamp >= start & Timestamp <= stop &
      DogID == handling$DogID[i] & CollarID == handling$CollarID[i]
  )
  setTxtProgressBar(pb, i)
  return(sub)
}) %>% do.call(rbind, .)

# Keep only data that is distinct
initial <- nrow(gps_dat)
gps_dat <- distinct(gps_dat)
cat("Removed", initial - nrow(gps_dat), "non-unique rows\n")

# We also want to know whether an individual was dispersing or resident. For
# this we need the dispersal dates.
dispersal <- read_csv("03_Data/01_RawData/POPECOL/DispersalDates.csv")

# Note that there are multiple dispersal periods for some individuals. Hence, we
# can't simply join the data by Dog. By default, let's assume that dogs are
# resident.
gps_dat$State <- "Resident"

# Go through the dogs and assign dispersal dates
names <- unique(dispersal$DogID)
for (i in seq_along(names)) {
  cutoff <- subset(dispersal, DogID == names[i])
  index <- which(gps_dat$DogID == names[i])
  for (h in 1:nrow(cutoff)) {
    gps_dat$State[index][gps_dat$Timestamp[index] >= cutoff$StartDate_UTC[h] &
    gps_dat$Timestamp[index] <= cutoff$EndDate_UTC[h]] <- "Disperser"
  }
}

# Check if there are any duplicates
dups_indices <- gps_dat %>%
  select(DogID, Timestamp) %>%
  duplicated()
table(dups_indices)

# We also need to remove any data where the coordinates are NA
gps_dat <- subset(gps_dat, !is.na(x) & !is.na(y))

################################################################################
#### Activity Data
################################################################################
# Identify all files containing activity data
files <- dir(path = "03_Data/01_RawData/POPECOL/ACT", full.names = T)

# Loop through the GPS files and clean them
cat("Loading and cleaning activity data...\n")
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
act_dat <- lapply(1:length(files), function(i) {

  # Extract DogID from filename so that we can put it into the dataframe
  dog <- unlist(strsplit(basename(files[i]), split = "_"))[2]

  # Identify column separator
  sep <- readLines(files[i], n = 1)
  sep <- if_else(grepl(x = sep, pattern = ";", useBytes = T), ";", ",")

  # Load the file as plain text and remove funny characters
  act_dat <- read_file(files[i], local = locale(encoding = "latin1"))
  if (sep == ",") {
      act_dat <- gsub(act_dat, pattern = ",", replacement = ";")
    } else {
      act_dat <- gsub(act_dat, pattern = ",", replacement = ".")
  }
  act_dat <- gsub(act_dat, pattern = "°", replacement = "")
  act_dat <- gsub(act_dat, pattern = "/", replacement = ".")
  act_dat <- gsub(act_dat, pattern = " \\[(.*?)\\]", replacement = "")

  # Load data into a data frame and do some more cleaning
  act_dat <- suppressMessages(read_delim(act_dat
      , local     = locale(encoding = "latin1")
      , delim     = ";"
      , col_types = cols()
      , col_names = T
    )) %>%

    # Add Dog Name
    mutate(DogID = dog) %>%

    # Retrieve timestamp
    mutate(Timestamp = as.POSIXct(
      paste(UTC_Date, UTC_Time), tz = "UTC", format = "%d.%m.%Y %H:%M:%S")
    ) %>%

    # Remove special characters like [°]
    setNames(gsub(names(.), pattern = " \\[*", replacement = "")) %>%
    setNames(gsub(names(.), pattern = "\\]", replacement = "")) %>%

    # Keep only desired columns
    dplyr::select(.
      , DogID     = DogID
      , CollarID  = CollarID
      , ActX      = ActivityX
      , ActY      = ActivityY
      , ActZ      = ActivityZ
      , Timestamp = Timestamp
    )

  # Return the filepath
  setTxtProgressBar(pb, i)
  return(act_dat)

}) %>% do.call(rbind, .)

# Let's load the collar handling dates
handling <- read_csv("03_Data/01_RawData/POPECOL/CollarHandling.csv")

# Loop through the collar handling dates and only keep valid data, i.e. data
# that falls within the range of dates after the first but before the last date
cat("Subsetting to valid data based on collar handling dates...\n")
pb <- txtProgressBar(min = 0, max = nrow(handling), style = 3)
act_dat <- lapply(1:nrow(handling), function(i) {
  start <- handling$StartDate_UTC[i]
  stop <- handling$EndDate_UTC[i]
  stop <- if_else(is.na(stop), ymd_hms(Sys.time(), tz = "UTC"), stop)
  sub <- subset(act_dat
    , Timestamp >= start & Timestamp <= stop &
      DogID == handling$DogID[i] & CollarID == handling$CollarID[i]
  )
  setTxtProgressBar(pb, i)
  return(sub)
}) %>% do.call(rbind, .)

# Keep only data that is distinct
initial <- nrow(act_dat)
act_dat <- distinct(act_dat)
cat("Removed", initial - nrow(act_dat), "non-unique rows\n")

# We also want to know whether an individual was dispersing or resident. For
# this we need the dispersal dates.
dispersal <- read_csv("03_Data/01_RawData/POPECOL/DispersalDates.csv")

# Note that there are multiple dispersal periods for some individuals. Hence, we
# can't simply join the data by Dog. By default, let's assume that dogs are
# resident.
act_dat$State <- "Resident"

# Go through the dogs and assign dispersal dates
names <- unique(dispersal$DogID)
for (i in seq_along(names)) {
  cutoff <- subset(dispersal, DogID == names[i])
  index <- which(act_dat$DogID == names[i])
  for (h in 1:nrow(cutoff)) {
    act_dat$State[index][act_dat$Timestamp[index] >= cutoff$StartDate_UTC[h] &
    act_dat$Timestamp[index] <= cutoff$EndDate_UTC[h]] <- "Disperser"
  }
}

# Check how many duplicates there are
dups_indices <- act_dat %>%
  select(DogID, Timestamp) %>%
  duplicated()
table(dups_indices)

# We also need to remove any data where the activity was NA
act_dat <- subset(act_dat, !is.na(ActX) & !is.na(ActY))

################################################################################
#### Merge Data
################################################################################
# Now we want to combine both the GPS and the activity data into a single
# dataframe. For this, we assign to every activity-fix the closest GPS-fix.
# Let's first make sure that the DogIDs and CollarIDs match.
gps_dat <- nest(gps_dat, GPS_Data = -c(DogID, CollarID))
act_dat <- nest(act_dat, ACT_Data = -c(DogID, CollarID))

# Merge them
dat <- left_join(gps_dat, act_dat, by = c("DogID", "CollarID"))

# Only keep entries where the activity data is not null
dat <- subset(dat, !sapply(ACT_Data, is.null))

# Within each row of the tibble, we now want to assign gps and activity data by
# their dates
cat("Joining activity and gps data...\n")
dat$MergedData <- pbmclapply(
    X                  = 1:nrow(dat)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {
  act <- dat$ACT_Data[[x]]
  gps <- dat$GPS_Data[[x]]
  indices <- sapply(act$Timestamp, function(y) {
    which.min(abs(y - gps$Timestamp))
  })
  gps_nearest <- gps[indices, ]
  gps_nearest <- select(gps_nearest, c(GPSTimestamp = Timestamp, x, y, DOP))
  merged <- cbind(act, gps_nearest)
  return(merged)
})

# Keep only the merged data and unnest
dat <- dat %>%
  select(c(DogID, CollarID, MergedData)) %>%
  unnest(MergedData)

# Compute the sunangle at each coordinate at the given time
sun <- getSunlightPosition(data = select(dat
  , date = Timestamp
  , lon = x
  , lat = y
))

# We want to work with the angle in degrees
dat$SunAngle <- sun$altitude * 180 / pi

# Let's use the angle to determine day and night periods
dat <- mutate(dat, ToD = case_when(
    SunAngle >= 0 ~ "Day"
  , SunAngle < 0 & SunAngle >= -6 ~ "CivilTwilight"
  , SunAngle < -6 & SunAngle >= -12 ~ "NauticalTwilight"
  , SunAngle < -12 & SunAngle >= -18 ~ "AstronomicalTwilight"
  , SunAngle < -18 ~ "Night"
))

# Let's also generate a binary indicator, assuming that crepuscular time is
# considered day
dat <- mutate(dat, ToD2 = ifelse(ToD == "Night", "Night", "Day"))

# Arrange such that the oldest data comes first
dat <- arrange(dat, DogID, CollarID, Timestamp)

# Add an indicator for the date
dat$Date <- date(dat$Timestamp)

# Nest the data by dog
dat <- dat %>% nest(data = -DogID)

# Identify tuples of night and day. That is, put together data collected during
# the day, plus the data of the following night. For this, we first grab the
# daily data and then add the data collected during the following night.
pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
dat$Cycles <- lapply(1:nrow(dat), function(i) {
  x <- dat$data[[i]]
  dates <- unique(x$Date)
  cycles <- lapply(dates, function(y) {

    # Get the data collected during the day on that date
    day <- x[x$Date == y & x$ToD2 == "Day", ]

    # Get the data collected on the following night (may be from the same or the
    # following day)
    night <- x[
      (x$Date == y & hour(x$Timestamp) > 12) |
      (x$Date == y + days(1) & hour(x$Timestamp) < 12), ]
    night <- night[night$ToD2 == "Night", ]

    # Put the data together and arrange by the timestamp
    cycle <- rbind(day, night)
    cycle <- arrange(cycle, Timestamp)

    # Return it
    setTxtProgressBar(pb, i)
    return(cycle)
  })

  # Put the cycles into a dataframe
  final <- tibble(Date = dates, Data = cycles, Cycle = 1:length(dates))
  final$Date <- NULL

  # Unnest the data
  final <- unnest(final, Data)

  # Return them
  return(final)
})

# Unnest the data and make cycle numbers unique across animals
dat <- dat %>%
  select(DogID, Cycles) %>%
  unnest(Cycles) %>%
  group_by(DogID, Cycle) %>%
  mutate(Cycle = cur_group_id()) %>%
  ungroup()

# Write the final data to file
write_csv(dat, "03_Data/02_CleanData/ActivityData.csv")
