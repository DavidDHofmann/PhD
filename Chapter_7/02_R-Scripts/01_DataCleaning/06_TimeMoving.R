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

################################################################################
#### Consecutive Activity Bouts
################################################################################
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

# Function that determines whether the dogs are moving or resting and returns
# the time of becoming active. It uses two parameters: "activity", which is the
# activity threshold. If an activity fix is above this threshold, the animals is
# said to be active. The second parameter is the "periods" parameter, which
# tells how many times animals need to be subsequently "active" to become
# "moving".
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
  return(data)
}

# Let's write another function that uses the function above, but then returns
# the first time at which the dogs become active in the evening
getFirstMovement <- function(data, activity, periods) {

  # Determine phases of activity and inactivity
  # activity <- 20
  # periods <- 6
  # data <- dat
  data <- getActivity(data, activity, periods)

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
design <- expand_grid(activity = seq(20, 160, 10), periods = c(2:6))
design$Results <- lapply(1:nrow(design), function(x) {
  first <- getFirstMovement(dat, activity = design$activity[x], periods = design$periods[x])
  return(first)
})

# # Unnest and plot
# design %>%
#   unnest(Results) %>%
#   ggplot(aes(x = yday(Date), y = as_hms(StartMoving))) +
#     geom_point(alpha = 0.5) +
#     geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
#     theme_minimal() +
#     facet_wrap(~ activity + periods)

# What is the percentage of fixes where individuals are active before 1400?
perc <- design %>%
  unnest(Results) %>%
  mutate(Before1400 = hour(StartMoving) < 14) %>%
  group_by(activity, periods) %>%
  count(Before1400) %>%
  pivot_wider(names_from = Before1400, values_from = n) %>%
  setNames(c("activity", "period", "After1400", "Before1400", "NeverActive")) %>%
  mutate(Total = sum(After1400, Before1400, NeverActive)) %>%
  mutate(
      After1400   = After1400 / Total
    , Before1400  = Before1400 / Total
    , NeverActive = NeverActive / Total
    , Minimize    = Before1400 + NeverActive
  )

# Visualize
ggplot(perc, aes(x = activity, y = Minimize, col = as.factor(period))) +
  geom_line() +
  theme_minimal()

# Find the value both values are minimal
perc %>% subset(Minimize == min(Minimize))

# Store them to file
write_csv(perc, "03_Data/03_Results/99_ActiveInactiveSeparation.csv")

# We'll go with 4 periods of over 50
first <- design[design$activity == 90 & design$periods == 3, ]
first <- unnest(first, Results)
first[, c("activity", "periods")] <- NULL

# How does the time of becoming active change throghout the year?
ggplot(first, aes(x = yday(Date), y = as_hms(StartMoving))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  theme_minimal()

# Combine with the activity data that was aggregated by "Time of Day"
act <- read_csv("03_Data/02_CleanData/ActivityDataCovariatesAggregated.csv")
act <- left_join(act, first, by = c("DogID", "CollarID", "Date", "State", "ToD"))

# Store combined data to file
write_csv(act, "03_Data/02_CleanData/ActivityDataCovariatesAggregatedTime.csv")
