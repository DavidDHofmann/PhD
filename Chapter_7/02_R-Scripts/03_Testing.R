################################################################################
####
################################################################################
# Determination of moonlight intensity

# Clear R's brain
rm(list = ls())

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_7")
# setwd("C:/Users/david/Switchdrive/University/15. PhD/Chapter_7")

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(hms)         # To handle times
library(pbmcapply)   # To run stuff in parallel
library(moonlit)     # For moonlight statistics

# Reload cleaned activity data (note that the timestamps are all in UTC)
dat <- read_csv("03_Data/02_CleanData/ActivityData.csv")
head(dat)

# Derive some time metrics
dat <- mutate(dat
  , Date = as_date(Timestamp)
  , Time = as_hms(Timestamp)
  , Hour = hour(Timestamp)
)

# We only want to work with the x-axis and remove the rest
dat$ActY <- dat$ActZ <- NULL

# Compute average 5 minute activity
means <- dat %>%
  group_by(Time) %>%
  summarize(
      Mean = mean(ActX)
    , SD   = sd(ActX)
  )

# Visualize them
ggplot(means, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  geom_hline(yintercept = 20, col = "gray30") +
  geom_vline(xintercept = as_hms("10:00:00"), col = "gray30") +
  geom_vline(xintercept = as_hms("12:00:00"), col = "gray30") +
  geom_vline(xintercept = as_hms("14:00:00"), col = "gray30") +
  geom_text(aes(x = as_hms("10:00:00"), label = "10:00:00", y = 150), angle = 90, size = 3, fontface = 3, nudge_x = -1000) +
  geom_text(aes(x = as_hms("12:00:00"), label = "12:00:00", y = 150), angle = 90, size = 3, fontface = 3, nudge_x = -1000) +
  geom_text(aes(x = as_hms("14:00:00"), label = "14:00:00", y = 150), angle = 90, size = 3, fontface = 3, nudge_x = -1000) +
  theme(axis.text.x = element_text(angle = 45))

# Visualize activity by hour and day
dat %>%
  group_by(Hour, Date) %>%
  summarize(ActX = mean(ActX)) %>%
  ggplot(aes(x = Hour, y = ActX, col = Date, group = Date)) +
    geom_line(lwd = 0.5, alpha = 0.5) +
    scale_color_viridis_c() +
    theme_minimal()

# # For now, focus on a couple of dogs
# dat <- subset(dat, DogID %in% sample(unique(dat$DogID), size = 5))

# Nest by dog
dat <- nest(dat, Data = -DogID)

# # Keep only a couple of entries per dog
# dat$Data <- lapply(dat$Data, function(x) {x[1:20000, ]})

# Check the time-lag between measurements
dat$Data <- pbmclapply(dat$Data, mc.cores = detectCores() - 1, ignore.interactive = T, function(x) {
  x$Lag <- difftime(x$Timestamp, lag(x$Timestamp), units = "mins")
  return(x)
})

# Unnest data
dat <- unnest(dat, Data)

# There should not be many datapoints with a time lag beyond 5 mins
dat %>%
  pull(Lag) %>%
  as.numeric() %>%
  summary()

# For each activity record, compute the moonlight intensity of the nearest night
pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
moonstats <- pbmclapply(1:nrow(dat), mc.cores = detectCores() - 1, ignore.interactive = T, function(x) {

  # Now calculate moonlight intensity of the following night
  pdf(file = NULL)
  moon <- calculateMoonlightStatistics(
      date     = dat$Timestamp[x]
    , lon      = dat$x[x]
    , lat      = dat$y[x]
    , e        = 0.21
    , t        = "15 mins"
    , timezone = "UTC"
  )
  dev.off()

  # Remove undesired data
  moon <- select(moon, -c(sunset, sunrise, date))

  # Return the data
  return(moon)

  # Update progress bar
  setTxtProgressBar(pb, value = x)

}) %>% do.call(rbind, .)

# Bind with original data
dat <- cbind(dat, moonstats)

# Write to file
write_csv(dat, "03_Data/02_CleanData/ActivityDataMoonphase.csv")

# Plot activity against moonlight intensity
ggplot(dat, aes(x = Time, y = maxMoonlightIntensity, z = ActX)) +
  stat_summary_2d(fun = "mean", bins = 100) +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal()

# Plot activity against date
ggplot(dat, aes(x = Time, y = Date, z = ActX)) +
  stat_summary_2d(fun = "mean", bins = 100) +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal()

# We deem an animal active whenever its activity is above 20
dat <- mutate(dat, Data = map(Data, function(x) {
  x$Active <- x$ActX > 20
  return(x)
}))

# Based on this, we now create a 30 min. moving window to find the time at which
# they start becoming active
periods <- 6
dat <- mutate(dat, Data = map(Data, function(x) {
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
        x$Status[i] <- x$Status[i-1]
      }
    }
  }
  return(x)
}))

# Find the time of first movement for every day
dat <- mutate(dat, FirstMovement = map(Data, function(x) {

  # Find all days through which we need to loop
  days <- unique(x$Date)
  firstmoves <- lapply(days, function(y) {
    first <- x[x$Date > y & x$Time >= hms::as_hms("14:00:00") & x$Status == "Moving", ]
    first <- first[1, ]
    return(first)
  }) %>% do.call(rbind, .)

  # Return the first moves
  firstmoves <- tibble(Day = days, FirstMovement = firstmoves$Timestamp)
  return(firstmoves)
}))
test <- dat %>% select(DogID, FirstMovement) %>% unnest(FirstMovement)

# visualize
ggplot(test, aes(x = hms::as_hms(FirstMovement))) + geom_histogram() + facet_wrap(~ DogID)







#### CONTINUE HERE!!!!





# Unnest again
dat <- unnest(dat, Data)

# Visualize status
ggplot(subset(dat, !is.na(Status)), aes(x = Time, y = Status, col = Date)) +
  geom_jitter(alpha = 0.5, pch = 20) +
  theme_minimal()

# Identify switchpoints and the direction of the switch
dat <- nest(dat, Data = -DogID)
dat <- mutate(dat, Data = map(Data, function(x) {
  x$Switch <- lag(x$Status) != x$Status
  x$SwitchTo <- NA
  indices <- which(x$Switch)
  for (i in indices) {
    if (x$Status[i-1] == "Resting" & x$Status[i] == "Moving") {
      x$SwitchTo[i] <- "ToMoving"
    } else if (x$Status[i-1] == "Moving" & x$Status[i] == "Resting") {
      x$SwitchTo[i] <- "ToResting"
    } else {
      x$SwitchTo[i] <- NA
    }
  }
  return(x)
}))

# Unnest again
dat <- unnest(dat, Data)

# Visualize the switches
ggplot(subset(dat, !is.na(SwitchTo)), aes(x = Time, y = SwitchTo)) +
  geom_jitter(pch = 20) +
  theme_minimal()

# Cut the data into distinct periods
dat <- nest(dat, Data = -DogID)
dat <- mutate(dat, Data = map(Data, function(x) {
  x$Period <- cumsum(ifelse(is.na(x$SwitchTo), 0, 1))
  x <- subset(x, Period > 0)
  return(x)
}))
dat <- unnest(dat, Data)

# Identify start and end of each period
per <- dat %>%
  group_by(DogID, Period) %>%
  summarize(
      PeriodStart = min(Timestamp)
    , PeriodEnd   = max(Timestamp)
    , Status      = unique(Status)
  )

dates <- unique(as_date(test$PeriodStart))
forplot <- lapply(dates, function(i) {
  sub <- subset(test, as_date(PeriodStart) == i | as_date(PeriodEnd) == i)
  if (as_date(sub$PeriodStart[1]) < i) {
    sub$PeriodStart[1] <- update(i, hours = 0, min = 0, sec = 1)
  }
  if (as_date(sub$PeriodEnd[nrow(sub)]) > i) {
    sub$PeriodEnd[nrow(sub)] <- update(i, hours = 23, min = 59, sec = 59)
  }
  return(sub)
}) %>% do.call(rbind, .)

forplot$StartDate <- as_date(forplot$PeriodStart)
forplot$StartTime <- hms::as_hms(forplot$PeriodStart)
forplot$EndDate <- as_date(forplot$PeriodEnd)
forplot$EndTime <- hms::as_hms(forplot$PeriodEnd)

ggplot(forplot, aes(x = StartTime, xend = EndTime, y = StartDate, yend = StartDate, col = Status)) +
  geom_segment() +
  facet_wrap(~ DogID, scales = "free")


ggplot(test, aes(x = StartTime, xend = EndTime, y = StartDate, yend = StartDate, col = Status)) +
  geom_segment()

ggplot(subset(dat, !is.na(Status)), aes(x = Time, y = Status, col = Period)) +
  geom_jitter(alpha = 0.5, pch = 20) +
  theme_minimal()

head(dat)

sub <- subset(dat, !is.na(SwitchTo))
sub$Period <- ifelse(sub$SwitchTo == "ToMoving", 1, 0)
sub$Period <- cumsum(sub$Period)

# CONTINUE HERE!!!
dat$Status <- pbmclapply(1:nrow(dat), function(x) {
  if (x < 6) {
    return(NA)
  } else {
    if (all(dat$Active[(x-5):x])) {
      return("Moving")
    } else if (all(!dat$Active[(x-5):x])) {
      return("Resting")
    } else {
      return(dat$Status[x-1])
    }
  }
}) %>% do.call(c, .)

ggplot(dat, aes(x = Hour, y = Moving)) + geom_jitter()

################################################################################
#### Exploration
################################################################################
# Get a feeling of the data
dim(dat)
nrow(dat) / 1e6
str(dat)
glimpse(dat)

# Histograms of activity data
dat %>%
  select(ActX, ActY, ActZ) %>%
  pivot_longer(ActX:ActZ, names_to = "Axis", values_to = "Activity") %>%
  ggplot(aes(x = Activity)) +
    geom_histogram(bins = 20) +
    facet_wrap(~ Axis, scales = "free")

# Check correlations of x and y axes
cor(dat$ActX, dat$ActY)

# It's clear that the Z axis is useless and that the Y and X axis are highly
# correlated. Thus, remove the z data and merge x and y data
dat$Act <- dat$ActX + dat$ActY
dat$ActX <- NULL
dat$ActY <- NULL
dat$ActZ <- NULL

# How about activity during the different times of the day?
dat %>%
  group_by(ToD) %>%
  summarize(Act = mean(Act)) %>%
  arrange(Act)

################################################################################
#### Average Moonlight Intensity per Night
################################################################################
# Let's compute average day and night activity for each cycle
avg <- dat %>%
  group_by(Cycle, ToD2) %>%
  summarize(
      DogID   = unique(DogID)
    , Date    = min(Date)
    , Act     = mean(Act)
    , .groups = "drop"
  ) %>%
  select(DogID, Date, Act, Cycle, ToD2) %>%
  pivot_wider(names_from = ToD2, values_from = Act)

# Let's also compute an "averaged" cycle location and join it to the data
avg <- dat %>%
  group_by(Cycle) %>%
  summarize(
      DogID   = unique(DogID)
    , Date    = min(Date)
    , x       = mean(x)
    , y       = mean(y)
    , .groups = "drop"
  ) %>%
  select(-c(DogID, Date)) %>%
  left_join(avg, ., by = c("Cycle"))

# We want to compute the lunar illumination for each night. We can use the
# moonlit package for this. Note, however, that the package computes nightly
# statistics for the night that is temporally closest. Hence, we need to amend
# the dates with a time in the afternoon to ensure that statistics are computed
# for the nights to come
avg$Timestamp <- as.POSIXct(avg$Date, tz = "UTC")
avg$Timestamp <- update(avg$Timestamp, hour = 20)

# Compute nightly statistics
moonstats <- pbmclapply(
    X                  = 1:nrow(avg)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {
  pdf(file = NULL)
  moon <- calculateMoonlightStatistics(
      date     = avg$Timestamp[x]
    , lon      = avg$x[x]
    , lat      = avg$y[x]
    , e        = 0.21
    , t        = "15 mins"
    , timezone = "UTC"
  )
  dev.off()
  return(moon)
}) %>% do.call(rbind, .)

# Only keep relevant data
moonstats <- select(moonstats, -c(sunset, sunrise, date))

# Bind the relevant data together
avg <- cbind(avg, moonstats)
avg$Timestamp <- NULL

# Store the data to file
write_csv(avg, "03_Data/02_CleanData/MoonlightData.csv")
