################################################################################
#### Analysis
################################################################################
# How does lunar intensity correlate with wild dog activity?

# Clear R's brain
rm(list = ls())

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_7")

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates

# Reload cleaned activity data
dat <- read_csv("03_Data/02_CleanData/MoonlightData.csv")

plot(dat$Night ~ dat$meanMoonlightIntensity)
ggplot(dat, aes(x = meanMoonlightIntensity, y = Night)) +
  geom_point() +
  scale_y_log10() +
  scale_x_sqrt()

ggplot(dat, aes(x = meanMoonlightIntensity, y = Day)) +
  geom_point() +
  scale_y_log10() +
  scale_x_sqrt()
