################################################################################
#### Compare Residents to Dispersers
################################################################################
# Description:

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(tidyverse)
library(amt)

# Load all wild dog data
dat <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_2/03_Data/01_RawData/POPECOL/Cleaned_GPSData.csv")

# Remove NAs
dat <- subset(dat, !is.na(x) & !is.na(y))

# Coerce data to tracks
dat <- make_track(dat, .x = x, .y = y, .t = Timestamp, .id = DogName, all_cols = T)
