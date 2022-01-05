################################################################################
#### Step Selection Function - Generation of Random Steps
################################################################################
# Description: In this script we coerce our gps data to steps and generate
# random steps for the (time varying) integrated step selection analysis.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load packages
library(tidyverse)    # For data wrangling
library(lubridate)    # To handle dates
library(pbmcapply)    # For multicore abilities
library(davidoff)     # Custom functions
library(terra)        # For handling spatial data
library(raster)       # For handling spatial data
library(sp)           # For handling spatial data
library(amt)          # To fit distributions

# Set seed for reproducability (need to change random sampler for parallel)
set.seed(1234)

################################################################################
#### Data Cleaning
################################################################################
# Load the gps data of dispersers and create rounded timestamps, then nest data
# by each dog
data <- "03_Data/02_CleanData/00_General_Dispersers.csv" %>%
  read_csv() %>%
  subset(State == "Disperser") %>%
  mutate(TimestampRounded = round_date(Timestamp, "1 hour")) %>%
  nest(data = -DogName)

# Resample GPS data to two (minimally) 2 hours
cat("Resampling GPS data to 2 hours...\n")
data$data <- suppressMessages(
  pbmclapply(1:nrow(data)
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x) {
      resFix(data$data[[x]], hours = 2, start = 1, tol = 0.5)
    }
  )
)

# Can unnest the data again
data <- data %>% unnest(data)

# Keep only fixes on our regular scheme
data <- subset(data, hour(TimestampRounded) %in% c(3, 7, 15, 19, 23))

# # Check for duplicates (there should be none)
# table(duplicated(data[, c("DogName", "TimestampRounded")]))

# We're going to convert the fixes to steps, yet we want to indicate a new
# burst if a step takes longer than 4.25 hours (or longer than 8.25 hours
# between 07:00 and 15:00). Let's thus calculate the timelags.
data <- data %>%
  group_by(DogName) %>%
  mutate(dt_ = Timestamp - lag(Timestamp)) %>%
  mutate(dt_ = as.numeric(dt_, units = "hours")) %>%
  ungroup()

# Indicate when a new burst starts (a new burst also starts whenever dt_ = NA)
# In case we're using iSSF
data <- data %>%
  mutate(NewBurst = ifelse(
    dt_ > 8.25 |
    (dt_ > 4.25 & hour(TimestampRounded) != 15) |
    is.na(dt_), yes = 1, no = 0)
  ) %>%
  mutate(BurstID = cumsum(NewBurst)) %>%
  dplyr::select(-c(NewBurst))

# We can only work with bursts that contain at least three fixes
data <- data %>%
  group_by(DogName, BurstID) %>%
  nest() %>%
  mutate(Nrow = map(data, nrow) %>% do.call(rbind, .)) %>%
  subset(Nrow >= 3) %>%
  unnest(data) %>%
  ungroup() %>%
  dplyr::select(-c(Nrow, dt_))

# Create burst IDs that are unique
data <- data %>%
  group_by(DogName, BurstID) %>%
  mutate(BurstID = cur_group_id()) %>%
  ungroup()

# Nest data by burst and compute step metrics
cat("Computing step metrics...\n")
tracks <- data %>%
  nest(data = -BurstID) %>%
  mutate(data = map(data, function(x) {
    coords <- cbind(x$x, x$y)
    coords <- vect(coords, crs = "+init=epsg:4326")
    coords <- project(coords, "+init=epsg:32734")
    coords <- as.data.frame(crds(coords))
    metrics <- step_met(coords$x, coords$y, x$Timestamp, degree = F)
    return(cbind(x, metrics))
  })) %>% unnest(data)

# Remove steps with missing turning angles
tracks <- subset(tracks, !is.na(relta_))

# Generate unique step id and indicate that steps are observed (case_) steps
tracks <- tracks %>% mutate(
    step_id_ = 1:nrow(.)
  , case_    = T
)

# Indicate if a fix was taken during time of activity or inactivity (we define
# inactivity as anything between 07:00 and 15:00)
tracks$inactive <- hour(tracks$TimestampRounded) == 7

# Remove unnecessary columns
tracks <- dplyr::select(tracks, -c(x, y, Timestamp, DOP, State, Source))

# # Visualize the steps
# ggplot(tracks, aes(x = x1_, y = y1_, col = BurstID, group = BurstID)) +
#   geom_path() +
#   geom_point() +
#   coord_sf() +
#   theme_minimal() +
#   theme(legend.position = "bottom")
#
# # Visualize step length
# ggplot(tracks, aes(x = sl_, col = inactive)) +
#   geom_density()

################################################################################
#### Generation of Random Steps
################################################################################
# We can't work with 0 step lengths or step speeds
tracks$sl_[tracks$sl_ == 0] <- 1

# Fit a Gamma distribution to the step speed (again, fitting a gamma to step
# speeds or fitting a gamma to (normalized) step lengths yields the same under
# iSSF) -> almost at least
sl <- fit_distr(tracks$sl_, "gamma")

# # Let's visualize the fit
# x <- seq(0, max(tracks$sl_), 1)
# y <- dgamma(x, scale = sl$params$scale, shape = sl$params$shape)
# hist(tracks$sl_, freq = F, breaks = 100)
# lines(y ~ x, col = "red")

# Write the distribution to file
dir.create("03_Data/03_Results", showWarnings = F)
write_rds(sl, "03_Data/03_Results/99_GammaDistribution.rds")

# Define the number of random steps
nsteps <- 24

# Now we use the fitted gamma distribution to sample step lengths. For the
# turning angles we follow Avgar et al. 2016 and use a uniform distribution
cat("Generating random steps...\n")
randomSteps <- pbmclapply(1:nrow(tracks), mc.cores = 1, ignore.interactive = T, function(i) {

  # Draw random turning angles
  relta_new <- runif(nsteps, min = -pi, max = pi)

  # Draw random step lengths with the fitted parameters
  sl_new <- rgamma(nsteps
    , shape = sl$params$shape
    , scale = sl$params$scale
  )

  # Put the step lengths and turning angles into a new dataframe and calculate the
  # new endpoints of each random step. We also indicate that the steps are control
  # steps (i.e. 0)
  steps <- tracks[rep(i, nsteps), ] %>%
    mutate(.
      , absta_  = absta_ + (relta_new - relta_)
      , relta_  = relta_new
      , sl_     = sl_new
      , case_   = 0
      , x2_     = x1_ + sin(absta_) * sl_
      , y2_     = y1_ + cos(absta_) * sl_
    )

  # Return the steps
  return(steps)
})

# Collapse the list of dataframes into a single dataframe
randomSteps <- do.call(rbind, randomSteps)

# We need to make sure that the absolute turning angle ranges from 0 to 2 * pi
randomSteps$absta_[randomSteps$absta_ > 2 * pi] <-
  randomSteps$absta_[randomSteps$absta_ > 2 * pi] - 2 * pi
randomSteps$absta_[randomSteps$absta_ < 0] <-
  randomSteps$absta_[randomSteps$absta_ < 0] + 2 * pi

# Merge the dataframes of the observed and random steps
ssf <- rbind(tracks, randomSteps) %>%
  arrange(step_id_, desc(case_)) %>%
  transform(case_ = as.logical(case_))

# Reorder the columns a bit
ssf <- select(ssf, DogName, Sex, BurstID, step_id_, everything(), -TimestampRounded)

# Store to file
write_csv(ssf, "03_Data/02_CleanData/00_General_Dispersers_SSF.csv")
ssf <- read_csv("03_Data/02_CleanData/00_General_Dispersers_SSF.csv")

# Coerce the steps to spatial lines
cat("Creating spatial lines from observed and random steps...\n")
lines <- lineTrack(ssf, CRS("+init=epsg:32734"))

# Transform the lines to WGS84
cat("Reprojecting steps to WGS84...\n")
lines <- spTransform(lines, CRS("+init=epsg:4326"))

# Coerce the duration column to a format that can be stored (numeric rather
# than difftime)
lines$dt_ <- as.numeric(lines$dt_)

# Store the lines
cat("Storing steps...\n")
writeOGR(lines
  , dsn       = "03_Data/02_CleanData"
  , layer     = filename
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

################################################################################
#### Session Information
################################################################################
# Create folder
dir.create("02_Rscripts/99_SessionInformation/02_Analysis", showWarnings = F)

# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/02_Analysis/00_SSF.rds")
cat("Done :)\n")
