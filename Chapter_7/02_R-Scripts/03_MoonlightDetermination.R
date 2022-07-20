################################################################################
#### Determining Moonlight Intensity During the Night
################################################################################
# Determination of moonlight intensity during the night closest to each activity
# burst

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
  , Date  = as_date(Timestamp)
  , Month = month(Timestamp)
  , Year  = year(Timestamp)
  , Time  = as_hms(Timestamp)
  , Hour  = hour(Timestamp)
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
