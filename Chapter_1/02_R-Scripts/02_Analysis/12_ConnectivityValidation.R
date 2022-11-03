################################################################################
#### Connectivity Validation
################################################################################
# Description: Here, we use independent dispersal data to validate our
# predictions of landscape connectivity. Specifically, we will use
# step-selection functions to determine whether dispersing individuals have a
# tendency to stick to areas of high connectivity.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(lubridate)      # To handle timestamps
library(wilddogr)       # To download dispersal data
library(amt)            # Some movement metrics
library(davidoff)       # Access to custom functions
library(terra)          # To handle spatial data

################################################################################
#### Download Additional Dispersal Data
################################################################################
# Let's download data from additional dispersers from dropbox
files <- dog_files(rvc = F)

# We are only interested in some of the dispersers
keep <- c("Appalachia", "Aspen", "Carson", "Chiounard", "Dell", "Earth"
  , "Encinitas", "Pula", "Rattan", "Ripley", "Saturday", "Sishen")
dispersers <- subset(files, DogName %in% keep)

# Download their data
downloaded <- dog_download(dispersers
  , outdir    = tempdir()
  , printpath = T
  , clean     = T
)

################################################################################
#### CONTINUE HERE!!!
################################################################################
# Something appears to be off with the step durations

################################################################################
#### Observed Steps
################################################################################
# Load the data and subset to their dispersal phases only. We will also resample
# the data to 4 hours and only keep those fixes that were collected at certain
# times. Finally, compute the duration of each step and determine bursts.
table(round(dat$dt_))
index <- round(dat$dt_) == 3
dat <- downloaded %>%
  read_csv() %>%
  subset(State == "Disperser") %>%
  resampleFixes(hours = 4, start = 1, tol = 0.5) %>%
  mutate(Timestamp = Timestamp - hours(2)) %>%
  mutate(TimestampRounded = round_date(Timestamp, "1 hour")) %>%
  subset(hour(TimestampRounded) %in% c(3, 7, 15, 19, 23)) %>%
  group_by(DogName) %>%
  mutate(dt_ = difftime(Timestamp, lag(Timestamp), units = "hours")) %>%
  ungroup() %>%
  mutate(NewBurst = ifelse(
    dt_ > 8.25 |
    (dt_ > 4.25 & hour(TimestampRounded) != 15) |
    is.na(dt_), yes = 1, no = 0)
  ) %>%
  mutate(BurstID = cumsum(NewBurst)) %>%
  dplyr::select(-c(NewBurst)) %>%
  subset(!is.na(x) & !is.na(y))

# We can only work with bursts that contain at least three fixes
dat <- dat %>%
  group_by(DogName, BurstID) %>%
  nest() %>%
  mutate(Nrow = map(data, nrow) %>% do.call(rbind, .)) %>%
  subset(Nrow >= 3) %>%
  unnest(data) %>%
  ungroup() %>%
  dplyr::select(-c(Nrow, dt_)) %>%
  group_by(DogName, BurstID) %>%
  mutate(BurstID = cur_group_id()) %>%
  ungroup()

# Use amt to run step selection analysis
ssf <- dat %>%
  make_track(.
    , .x      = x
    , .y      = y
    , .t      = Timestamp
    , id      = BurstID
    , crs     = CRS("+init=epsg:4326")
    , State   = State
    , DogName = DogName
  ) %>%
  transform_coords(CRS("+init=epsg:32734")) %>%
  nest(data = -"id") %>%
  mutate(data = map(data, function(x){
    x %>%
      steps(keep_cols = "start") %>%
      mutate(., dt_ = as.numeric(dt_, units = "hours"))
  })) %>%
  unnest(data) %>%
  random_steps(n_control = 24)
