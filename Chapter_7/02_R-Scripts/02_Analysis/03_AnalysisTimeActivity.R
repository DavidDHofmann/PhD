################################################################################
#### Analysis
################################################################################
# Description: Analysis of the Activity data

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
# wd <- "C:/Users/david/switchdrive/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(hms)         # To handle times
# library(runner)      # To apply moving windows

# Load cleaned activity data (note that the timestamps are all in UTC)
dat <- read_csv("03_Data/02_CleanData/ActivityDataWithCovariates.csv")

############################################################
#### STOP HERE
############################################################
# For now, subset the data
dat_backup <- dat
dat <- subset(dat, DogID %in% c("Abel", "Oolong"))

# We now want to assess at what time of the day the dogs usually become active.
# Hence, we want to nest the data by dog/collar first
dat <- nest(dat, Data = -c("DogID", "CollarID", "BurstID"))
print(dat)

# Based on this, we now create a 30 min (i.e. 6 activity periods) moving window
# to find the time at which they start becoming active
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

# Visualize status
dat %>%
  unnest(Data) %>%
  subset(!is.na(Status)) %>%
  ggplot(aes(x = Time, y = Status, col = Date)) +
    geom_jitter(alpha = 0.5, pch = 20) +
    theme_minimal() +
    scale_color_viridis_c(option = "magma")

# Identify switchpoints at which individuals go from moving to resting or from
# resting to moving. Also determine the direction of the switch
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

# Visualize the switches
dat %>%
  unnest(Data) %>%
  subset(!is.na(SwitchTo)) %>%
  ggplot(aes(x = Time, y = SwitchTo)) +
    geom_jitter(pch = 20) +
    theme_minimal()

# Let's figure out at what time of the day the individuals start moving for the
# first time in the evening (afternoon)
first <- mutate(dat, FirstMovement = map(Data, function(x) {

  # Find all days through which we need to loop
  days <- unique(x$Date)
  firstmoves <- lapply(days, function(y) {
    # first <- x[x$Date == y & x$Time > phase$Endtime[phase$Burst == "Morning"] & x$Status == "Moving", ]
    first <- x[x$Date == y & x$Time > as_hms("12:00:00") & x$Status == "Moving", ]
    first <- first[1, ]
    return(first)
  }) %>% do.call(rbind, .)

  # Return the first moves
  firstmoves <- tibble(Day = days, FirstMovement = firstmoves$Timestamp)
  return(firstmoves)
}))

# Let's visualize them
first %>%
  select(DogID, FirstMovement) %>%
  unnest(FirstMovement) %>%
  subset(!is.na(FirstMovement)) %>%
  ggplot(aes(x = as_hms(FirstMovement))) +
    geom_histogram(bins = 30) +
    facet_wrap(~ DogID) +
    theme_minimal()

# Unnest the data again
dat <- unnest(dat, Data)

############ CONTINUE HERE!!!!

# Join this data with nightly statistics
stats <- dat %>%
  subset(Hour > 14) %>%
  select(DogID, Date, maxMoonlightIntensity) %>%
  mutate(maxMoonlightIntensity = round(maxMoonlightIntensity, 5)) %>%
  distinct()
test <- first %>%
  select(DogID, FirstMovement) %>%
  unnest(FirstMovement) %>%
  subset(!is.na(FirstMovement)) %>%
  subset(hour(FirstMovement) > 14) %>%
  left_join(., stats, by = c("DogID", "Day" = "Date"))
test$FirstMovement <- as.numeric(minutes(as_hms(test$FirstMovement)))
plot(as_hms(test$FirstMovement) ~ test$maxMoonlightIntensity)
abline(mod)
mod <- lm(FirstMovement ~ maxMoonlightIntensity, data = test)
summary(mod)
