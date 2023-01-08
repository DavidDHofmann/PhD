################################################################################
#### Analysis of Number of Activity Bursts
################################################################################
# Description: Analysis of the number of activity bursts depending on moonlight

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(hms)         # To handle times

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load data
dat <- read_rds("03_Data/02_CleanData/ActivityDataCovariatesStates.rds")
dat$ModelParams <- NULL
dat <- dat[1:2, ]

# Create tuples of night and day
dat <- mutate(dat, Data = map(Data, function(x) {
  x$Switch <- x$ToD2 != lag(x$ToD2)
  x$Switch[1] <- F
  x$Switch <- x$Switch & x$ToD2 == "Day"
  x$CycleDay  <- cumsum(x$Switch) + 1
  x$Switch <- NULL
  return(x)
}))
sub <- dat$Data[[1]]
sub$CycleDay

# Nest by ID and cycle
dat <- dat %>%
  unnest(Data) %>%
  nest(Data = -c(DogID, CycleDay))

# Per cycle, we are only interested in activity after 12:00 UTC (14:00 local
# time)
dat <- mutate(dat, Data = map(Data, function(x) {
  fir <- min(x$Date)
  fir <- ymd_hms(paste0(fir, " 12:00:00"))
  sub <- subset(x, Timestamp >= fir)
  return(sub)
}))

# Keep only entries where we have more than 0 rows
dat$Nrow <- sapply(dat$Data, nrow)
dat <- subset(dat, Nrow > 0)
dat$Nrow <- NULL

# Each tuple will be assigned the date of the start of the tuple
dat$Date <- lapply(dat$Data, function(x) {
  min(x$Date)
}) %>% do.call(c, .)

# Within each cycle, identify the number of activity bursts
dat <- mutate(dat, Data = map(Data, function(x) {
  x$Switch <- x$State != lag(x$State)
  x$Switch[1] <- F
  x$CycleState  <- cumsum(x$Switch) + 1
  x$Switch <- NULL
  return(x)
}))

# Now we can assess how often each animal was active in each cycle
dat <- mutate(dat, Nact = map_int(Data, function(x) {
  nactive <- x %>%
    subset(State == 3) %>%
    select(CycleState) %>%
    distinct() %>%
    nrow()
  return(nactive)
}))

# Add nightly stats
moon <- read_csv("03_Data/02_CleanData/Moonlight.csv")
dat <- left_join(dat, moon, by = c("Date" = "date"))
print(names(dat))

subset(dat, maxMoonlightIntensity > 0.2) %>%
  slice(1) %>%
  pull(Date)

plot(Nact ~ maxMoonlightIntensity, data = dat)
par(mfrow = c(2, 1))
summary(dat$N)
hist(dat$Nact[dat$maxMoonlightIntensity < 0.2], freq = F)
hist(dat$Nact[dat$maxMoonlightIntensity >= 0.2], freq = F)
boxplot(dat$Nact[dat$maxMoonlightIntensity < 0.2])
boxplot(dat$Nact[dat$maxMoonlightIntensity >= 0.2])

# Can do some plots using the function below (note that a darker moon means that
# it is more illuminated)
plotActivityMoon <- function(index) {
  subdat <- dat_nested %>%
    slice(index) %>%
    unnest(Data)
  p1 <- ggplot(subdat, aes(x = Timestamp, y = MoonAngle)) +
      geom_hline(yintercept = 0, col = "gray", lty = 2) +
      geom_line() +
      geom_moon(data = subdat[10, ], aes(x = Timestamp, y = 0, ratio = 1), col = "black", fill = "white") +
      geom_moon(data = subdat[10, ], aes(x = Timestamp, y = 0, ratio = MoonPhase), fill = "black") +
      theme_minimal() +
      theme(panel.grid.minor = element_blank()) +
      ylim(c(-90, 90))
  p2 <- ggplot(subdat, aes(x = Timestamp, y = Act, col = as.factor(State))) +
      geom_point() +
      geom_vline(aes(xintercept = max(Sunset))) +
      geom_vline(aes(xintercept = max(Sunset) - hours(2)), lty = 2) +
      geom_vline(aes(xintercept = max(Sunset) + hours(2)), lty = 2) +
      geom_vline(aes(xintercept = max(Sunrise))) +
      geom_vline(aes(xintercept = max(Sunrise) - hours(2)), lty = 2) +
      geom_vline(aes(xintercept = max(Sunrise) + hours(2)), lty = 2) +
      theme_minimal() +
      scale_color_viridis_d() +
      theme(legend.position = c(0.1, 0.8)) +
      xlab("") +
      ylim(c(0, 510))
  ggarrange(p2, p1, nrow = 2, align = "hv", heights = c(0.8, 0.2))
}

# Try it
dat_nested <- select(dat, DogID, Data)
library(gggibbous)
library(ggpubr)
plotActivityMoon(index = 90)
