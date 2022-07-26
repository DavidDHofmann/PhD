################################################################################
#### Testing
################################################################################
# Testing

# Clear R's brain
rm(list = ls())

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_7")
# setwd("C:/Users/david/Switchdrive/University/15. PhD/Chapter_7")

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(hms)         # To handle times
library(pracma)      # To identify peaks

# Reload cleaned activity data (note that the timestamps are all in UTC)
dat <- read_csv("03_Data/02_CleanData/ActivityDataMoonphase.csv")

################################################################################
#### NEED TO APPLY SOME FILTERS
################################################################################
# Look at residents only?
# What to do with individuals that were together?
dat <- subset(dat, State != "Disperser")
gc()

# Let's get an idea of the temporal resolution of the data
dat %>%
  count(Month, Year) %>%
  ggplot(aes(x = Year, y = Month, fill = n)) +
    geom_raster() +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal()

# Let's get an idea of the usual activity patterns by computing average activity
# for the different times.
means_hour <- dat %>%
  group_by(Time) %>%
  summarize(
      Mean    = mean(ActX)
    , SD      = sd(ActX)
    , .groups = "drop"
  ) %>%
  ungroup()

# Visualize it
ggplot(means_hour, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

# We can try to identify the peaks and valleys in this plot
peaks <- findpeaks(means_hour$Mean, minpeakdistance = 50, npeaks = 2)
peaks <- subset(means_hour, Mean %in% peaks[, 1])
valls <- findpeaks(-means_hour$Mean, minpeakdistance = 50, npeaks = 4)
valls <- subset(means_hour, -Mean %in% valls[, 1])
valls <- valls[-nrow(valls), ]

# Let's use them to distinguish morning and evening activity (i'll make the
# endtime for the morning burst a bit earlier to have a gap between morning and
# evening burst)
phase <- tibble(
    Burst     = factor(c("Morning", "Evening"), levels = c("Morning", "Evening"))
  , Starttime = c(valls$Time[1], valls$Time[2])
  , Endtime   = c(valls$Time[2] - duration(2, units = "hours"), valls$Time[3])
)

# Visualize again
ggplot(means_hour, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_rect(data = phase, inherit.aes = F, aes(xmin = Starttime, xmax = Endtime, ymin = -Inf, ymax = Inf, fill = Burst), alpha = 0.15) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_vline(data = peaks, aes(xintercept = Time), col = "gray20", lty = 2) +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = "bottom") +
  scale_fill_manual(values = c("cornflowerblue", "orange"))

# How does the activity pattern change depending on the month?
means_month <- dat %>%
  group_by(Month, Time) %>%
  summarize(
      Mean    = mean(ActX)
    , SD      = sd(ActX)
    , .groups = "drop"
  ) %>%
  ungroup()

# Find peaks each month
peaks <- means_month %>%
  nest(Data = -c(Month)) %>%
  mutate(Data = map(Data, function(x) {
    peaks <- findpeaks(x$Mean, minpeakdistance = 50, npeaks = 2)
    peaks <- subset(x, Mean %in% peaks[, 1])
    return(peaks)
  })) %>%
  unnest(Data)

# Visualize
ggplot(means_month, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_rect(data = phase, inherit.aes = F, aes(xmin = Starttime, xmax = Endtime, ymin = -Inf, ymax = Inf, fill = Burst), alpha = 0.15) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_vline(data = peaks, aes(xintercept = Time), col = "gray20", lty = 2) +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  facet_wrap(~ Month, ncol = 2, dir = "v") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = "bottom") +
  scale_fill_manual(values = c("cornflowerblue", "orange"))

# How does the activity pattern change depending on the year?
means_year <- dat %>%
  group_by(Year, Time) %>%
  summarize(
      Mean    = mean(ActX)
    , SD      = sd(ActX)
    , .groups = "drop"
  ) %>%
  ungroup()

# Find peaks each month
peaks <- means_year %>%
  nest(Data = -c(Year)) %>%
  mutate(Data = map(Data, function(x) {
    peaks <- findpeaks(x$Mean, minpeakdistance = 50, npeaks = 2)
    peaks <- subset(x, Mean %in% peaks[, 1])
    return(peaks)
  })) %>%
  unnest(Data)

# Visualize
ggplot(means_year, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_rect(data = phase, inherit.aes = F, aes(xmin = Starttime, xmax = Endtime, ymin = -Inf, ymax = Inf, fill = Burst), alpha = 0.15) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_vline(data = peaks, aes(xintercept = Time), col = "gray20", lty = 2) +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  facet_wrap(~ Year, ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = "bottom") +
  scale_fill_manual(values = c("cornflowerblue", "orange"))

# Use the cutoff times to assign to each activity fix "morning", "evening", or
# NA
dat <- dat %>% mutate(ActivityPhase = case_when(
    Time >= phase$Starttime[1] & Time <= phase$Endtime[1] ~ "Morning"
  , Time >= phase$Starttime[2] & Time <= phase$Endtime[2] ~ "Evening"
  , TRUE ~ "Ignore"
))

# Make sure it worked
sort(unique(dat[dat$ActivityPhase == "Morning", ]$Hour))
sort(unique(dat[dat$ActivityPhase == "Evening", ]$Hour))
sort(unique(dat[dat$ActivityPhase == "Ignore", ]$Hour))

# Let's compute average activity during those
sub <- subset(dat, ActivityPhase != "Ignore") %>%
  group_by(DogID, CollarID, Date, ActivityPhase) %>%
  summarize(
      meanActivity           = mean(ActX)
    , totalActivity          = sum(ActX)
    , minMoonlightIntensity  = min(minMoonlightIntensity)
    , meanMoonlightIntensity = mean(meanMoonlightIntensity)
    , maxMoonlightIntensity  = max(maxMoonlightIntensity)
    , .groups                = "drop"
  )

# Visualize
ggplot(sub, aes(x = maxMoonlightIntensity, y = meanActivity)) +
  geom_point(alpha = 0.2) +
  theme_minimal() +
  facet_wrap(~ ActivityPhase) +
  stat_smooth(method = "lm")
summary(lm(meanActivity ~ maxMoonlightIntensity, data = subset(sub, ActivityPhase == "Morning")))
summary(lm(meanActivity ~ maxMoonlightIntensity, data = subset(sub, ActivityPhase == "Evening")))

################################################################################
#### CONTINUE HERE!!!
################################################################################
# Let's compute the average activity during the morning and evening bursts for
# each day and individual
dat_summarized <- dat %>%
  group_by(DogID, ActivityPhase, Date) %>%
  summarize(
      Mean                  = mean(ActX)
    , SD                    = sd(ActX)
    , maxMoonlightIntensity = max(maxMoonlightIntensity)
    , .groups               = "drop"
  ) %>%
  ungroup()
hist(dat_summarized$Mean)
summary(lm(Mean ~ maxMoonlightIntensity, data = subset(test, ActivityPhase == "Evening")))
ggplot(test, aes(x = maxMoonlightIntensity, y = Mean)) +
  geom_point() +
  geom_smooth()

# Now let's try to figure out if this pattern changes depending on the moon
# illumination. For this, split the data into dark / illuminated nights based on
# the quantiles
dat %>%
  select(minMoonlightIntensity, meanMoonlightIntensity, maxMoonlightIntensity) %>%
  apply(2, summary) %>%
  round(6)
dat %>%
  sample_n(size = 1e6, replace = F) %>%
  select(minMoonlightIntensity, meanMoonlightIntensity, maxMoonlightIntensity) %>%
  pivot_longer(cols = 1:3, names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = Variable, y = Value)) +
    geom_boxplot()

# Let's split the data into "dark" and "illuminated" nights
quantiles <- quantile(dat$maxMoonlightIntensity, c(0.10, 0.90))
dat <- mutate(dat, NightIllumination = case_when(
    maxMoonlightIntensity < quantiles[1] ~ "Dark"
  , maxMoonlightIntensity > quantiles[2] ~ "Illuminated"
  , TRUE ~ "Ignore"
))

# Compute averages again
means <- dat %>%
  subset(NightIllumination != "Ignore") %>%
  group_by(Time, NightIllumination) %>%
  summarize(
      Mean    = mean(ActX)
    , SD      = sd(ActX)
    , .groups = "drop"
  ) %>%
  ungroup()

# Visualize them above each other
ggplot(means, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD, col = NightIllumination, fill = NightIllumination)) +
  geom_ribbon(alpha = 0.1, lwd = 0.2) +
  geom_line() +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  theme(axis.text.x = element_text(angle = 45)) +
  facet_wrap(~ NightIllumination, nrow = 2) +
  theme(legend.position = "bottom")

  # Visualize them on top of each other
ggplot(means, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD, col = NightIllumination, fill = NightIllumination)) +
  geom_ribbon(alpha = 0.1, lwd = 0.2) +
  geom_line() +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position = "bottom")

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
