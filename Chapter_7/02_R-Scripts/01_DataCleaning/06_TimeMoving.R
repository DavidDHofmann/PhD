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

# Function to assess activity depending on an activity threshold that has to be
# exceeded for a minimal number of periods
getActivity <- function(data, activity, periods) {
  # data <- dat
  # activity <- 20
  # periods <- 6
  data$Data <- pbmclapply(
      X                  = data$Data
    , mc.cores           = detectCores() - 1
    , ignore.interactive = T
    , FUN                = function(x) {
      x$Active <- x$ActX > activity
      x$Status <- NA
      for (i in periods:nrow(x)) {
        if (as.numeric(x$Lag[i]) > 5 | is.na(x$Lag[i])) {
          x$Status[i] <- NA
        } else {
          if (all(x$Active[(i - (periods - 1)):i])) {
            x$Status[(i - (periods - 1)):i] <-"Moving"
          } else if (all(!x$Active[(i - (periods - 1)):i])) {
            x$Status[(i - (periods - 1)):i] <- "Resting"
          } else {
            x$Status[i] <- x$Status[i - 1]
          }
        }
      }
    return(x)
  })

  # Each day, determine time of the day when the dogs became active first
  data$FirstMovement <- pbmclapply(
      X                  = data$Data
    , mc.cores           = detectCores() - 1
    , ignore.interactive = T
    , FUN                = function(x) {
    first <- subset(x, Status == "Moving")[1, ]
    first <- select(first, Timestamp)
    return(first)
  })

  # Unnest
  first <- data %>%
    select(-c(Data, N)) %>%
    unnest(FirstMovement)

  # Generate a column indicating that the timing relates to the evening burst
  first$ToD <- "Evening"
  names(first)[names(first) == "Timestamp"] <- "StartMoving"
  return(first)
}

# Try it for different values
design <- expand_grid(activity = c(20, 50, 100), periods = c(4, 6, 8))
design$Results <- lapply(1:nrow(design), function(x) {
  first <- getActivity(dat, activity = design$activity[x], periods = design$periods[x])
  return(first)
})

# Unnest and plot
design %>%
  unnest(Results) %>%
  ggplot(aes(x = yday(Date), y = as_hms(StartMoving))) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
    theme_minimal() +
    facet_wrap(~ activity + periods)

# What is the percentage of fixes where individuals are active before 1400?
design %>%
  unnest(Results) %>%
  mutate(Before1400 = hour(StartMoving) < 14) %>%
  group_by(activity, periods) %>%
  count(Before1400) %>%
  pivot_wider(names_from = Before1400, values_from = n) %>%
  setNames(c("activity", "period", "After1400", "Before1400", "NeverActive")) %>%
  mutate(PercentBefore1400 = Before1400 / sum(Before1400, After1400))

# # We now create a 30 min. moving window to find the time at which time the dogs
# # start moving. Specifically, the dogs are deemed to start moving once all fixes
# # inside the moving window are "active". Similarly, the dogs are resting once
# # all fixes become "inactive"
# activity <- 20
# periods  <- 6
# dat$Data <- pbmclapply(
#     X                  = dat$Data
#   , mc.cores           = detectCores() - 1
#   , ignore.interactive = T
#   , FUN                = function(x) {
#     x$Active <- x$ActX > activity
#     x$Status <- NA
#     for (i in periods:nrow(x)) {
#       if (as.numeric(x$Lag[i]) > 5 | is.na(x$Lag[i])) {
#         x$Status[i] <- NA
#       } else {
#         if (all(x$Active[(i - (periods - 1)):i])) {
#           x$Status[(i - (periods - 1)):i] <-"Moving"
#         } else if (all(!x$Active[(i - (periods - 1)):i])) {
#           x$Status[(i - (periods - 1)):i] <- "Resting"
#         } else {
#           x$Status[i] <- x$Status[i - 1]
#         }
#       }
#     }
#   return(x)
# })
#
# # Each day, determine time of the day when the dogs became active first
# dat$FirstMovement <- pbmclapply(
#     X                  = dat$Data
#   , mc.cores           = detectCores() - 1
#   , ignore.interactive = T
#   , FUN                = function(x) {
#   first <- subset(x, Status == "Moving")[1, ]
#   first <- select(first, Timestamp)
#   return(first)
# })
#
# # Unnest
# first <- dat %>%
#   select(-c(Data, N)) %>%
#   unnest(FirstMovement)
#
# # Generate a column indicating that the timing relates to the evening burst
# first$ToD <- "Evening"
# names(first)[names(first) == "Timestamp"] <- "StartMoving"

# How does the time of becoming active change throghout the year?
ggplot(first, aes(x = yday(Date), y = as_hms(StartMoving))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  theme_minimal()

# Combine with the activity data that was aggregated by "Time of Day"
act <- read_csv("03_Data/02_CleanData/ActivityDataCovariatesAggregated.csv")
act <- left_join(act, first, by = c("DogID", "CollarID", "Date", "State", "ToD"))

# Store combined data to file
write_csv(act, "03_Data/02_CleanData/ActivityDataCovariatesAggregated.csv")

################################################################################
#### CONTINUE HERE
################################################################################
# For this, we'll only work with "evening" data
evening <- subset(act, ToD == "Evening")

# Compute time of the day as regular seconds since midnight
evening$StartMovingNumeric <- as.numeric(as_hms(evening$StartMoving))

# Nest data by individual and tod
evening_nested <- evening %>% nest(Data = -c(DogID, ToD))

# Only work with rows for which there are at least 50 datapoints
evening_nested$N <- sapply(evening_nested$Data, function(x) {nrow(x)})
evening_nested <- subset(evening_nested, N >= 50)
print(evening_nested)

# How does the time of becoming active change throghout the year?
ggplot(evening, aes(x = yday(Date), y = as_hms(StartMoving))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  theme_minimal()
