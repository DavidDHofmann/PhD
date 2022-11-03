################################################################################
#### Time when Dogs Start Moving
################################################################################
# Description: Take the raw activity data to determine at what time in the
# afternoon start moving

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(hms)         # To handle times
library(pbmcapply)   # To run stuff in parallel

# Load cleaned activity data (note that the timestamps are all in UTC)
dat <- read_csv("03_Data/02_CleanData/ActivityDataCovariates.csv")

# Compute some useful time metrics
dat$Hour <- as_hms(dat$Timestamp)

# We deem an animal active if ActX >= 20
dat$Active <- dat$ActX > 20

# We are only interested in data that is collected after 12:00 (14:00 Botswana
# time)
dat <- subset(dat, Hour >= as_hms("12:00:00"))

# Nest the data by dog, collar, state, and date
dat <- dat %>% nest(Data = -c(DogID, CollarID, State, Date))

# Check the number of datapoints per day
dat$N <- sapply(dat$Data, function(x) {nrow(x)})

# Let's only keep days where we have at least 140 datapoints
dat <- subset(dat, N > 140)

# Compute temporal lag between fixes (the lag of the first fix doesn't matter,
# so set it to 5 minutes)
dat$Data <- pbmclapply(
    X                  = dat$Data
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(x) {
    x$Lag <- difftime(x$Timestamp, lag(x$Timestamp), units = "mins")
    x$Lag[1] <- 5
    return(x)
})

# We now create a 30 min. moving window to find the time at which time the dogs
# start moving. Specifically, the dogs are deemed to start moving once all fixes
# inside the moving window are "active". Similarly, the dogs are resting once
# all fixes become "inactive"
periods <- 6
dat$Data <- pbmclapply(
    X                  = dat$Data
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(x) {
    x$Status <- NA
    for (i in 1:nrow(x)) {
      if (i < periods | as.numeric(x$Lag[i]) > 5 | is.na(x$Lag[i])) {
        x$Status[i] <- NA
      } else {
        if (all(x$Active[(i - (periods - 1)):i])) {
          x$Status[i] <-"Moving"
        } else if (all(!x$Active[(i - (periods - 1)):i])) {
          x$Status[i] <- "Resting"
        } else {
          x$Status[i] <- x$Status[i - 1]
        }
      }
    }
  return(x)
})

# Each day, determine time of the day when the dogs became active first
dat$FirstMovement <- pbmclapply(
    X                  = dat$Data
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(x) {
  first <- subset(x, Status == "Moving")[1, ]
  first <- select(first, Timestamp)
  return(first)
})

# Unnest
first <- dat %>%
  select(-c(Data, N)) %>%
  unnest(FirstMovement)

# Generate a column indicating that the timing relates to the evening burst
first$ToD <- "Evening"
names(first)[names(first) == "Timestamp"] <- "StartMoving"

# Combine with the activity data that was aggregated by "Time of Day"
act <- read_csv("03_Data/02_CleanData/ActivityDataCovariatesAggregated.csv")
act <- left_join(act, first, by = c("DogID", "CollarID", "Date", "State", "ToD"))

# Store combined data to file
write_csv(act, "03_Data/02_CleanData/ActivityDataCovariatesAggregated.csv")
