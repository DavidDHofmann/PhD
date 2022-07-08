################################################################################
#### Moolight Intensity
################################################################################
# Determination of moonlight intensity

# Clear R's brain
rm(list = ls())

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_7")

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(moonlit)     # For moonlight statistics
library(pbmcapply)   # To run stuff in parallel

# Reload cleaned activity data
dat <- read_csv("03_Data/02_CleanData/ActivityData.csv")

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
