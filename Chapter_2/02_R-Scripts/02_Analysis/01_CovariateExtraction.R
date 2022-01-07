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
library(fitdistrplus) # To fit distributions
library(ggpubr)       # To put multiple plots together

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

# Check out the resampled data. Note that for Ripley there are only three fixes
print(data, n = 30)

# Unnest the data again
data <- data %>% unnest(data)

# Keep only fixes on our regular scheme
data <- subset(data, hour(TimestampRounded) %in% c(3, 7, 15, 19, 23))

# Check for duplicates (there should be none)
table(duplicated(data[, c("DogName", "TimestampRounded")]))

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
tracks <- dplyr::select(tracks, -c(Timestamp, DOP, State, Source))

# Visualize the steps
ggplot(tracks, aes(x = x, y = y, group = BurstID)) +
  geom_path(size = 0.1) +
  geom_point(size = 0.5) +
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "none")

# Visualize step lengths and turning angles during activity and inactivity
p1 <- ggplot(tracks, aes(x = sl_, col = inactive, fill = inactive)) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  ggtitle("Step Lengths") +
  xlab("Step Length (m)") +
  ylab("Density")
p2 <- ggplot(tracks, aes(x = relta_, col = inactive, fill = inactive)) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  ggtitle("Turning Angles") +
  xlab("Turning Angle (rad)") +
  ylab("Density")

# Visualize step lengths and turning angles for different sexes
p3 <- ggplot(tracks, aes(x = sl_, col = Sex, fill = Sex)) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  ggtitle("Step Lengths") +
  xlab("Step Length (m)") +
  ylab("Density")
p4 <- ggplot(tracks, aes(x = relta_, col = Sex, fill = Sex)) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  ggtitle("Turning Angles") +
  xlab("Turning Angle (rad)") +
  ylab("Density")

# Put plots together
ggarrange(p1, p2, p3, p4)

################################################################################
#### Generation of Random Steps
################################################################################
# We can't work with 0 step lengths so we'll force each step to at least 1m
tracks$sl_[tracks$sl_ == 0] <- 1

# Fit a Gamma distribution to the step speed (again, fitting a gamma to step
# speeds or fitting a gamma to (normalized) step lengths yields the same under
# iSSF) -> almost at least
# sl <- fit_distr(tracks$sl_, "gamma")
sl <- fitdist(tracks$sl_, "gamma", method = "mle", lower = 0)
sl <- list(
    scale = 1 / sl$estimate[["rate"]]
  , shape = sl$estimate[["shape"]]
)

# Let's visualize the fit
x <- seq(0, max(tracks$sl_), 1)
y <- dgamma(x, scale = sl$scale, shape = sl$shape)
hist(tracks$sl_, freq = F, breaks = 100)
lines(y ~ x, col = "red")

# Write the distribution to file
dir.create("03_Data/03_Results", showWarnings = F)
write_rds(sl, "03_Data/03_Results/99_GammaDistribution.rds")

# Define the number of random steps
nsteps <- 100

# Now we use the fitted gamma distribution to sample step lengths. For the
# turning angles we follow Avgar et al. 2016 and use a uniform distribution
cat("Generating random steps...\n")
pb <- txtProgressBar(min = 0, max = nrow(tracks), style = 3)
randomSteps <- lapply(1:nrow(tracks), function(i) {

  # Draw random turning angles
  relta_new <- runif(nsteps, min = -pi, max = pi)

  # Draw random step lengths with the fitted parameters
  sl_new <- rgamma(nsteps
    , shape = sl$shape
    , scale = sl$scale
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
  setTxtProgressBar(pb, i)
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

# I'd like to work exclusively with lon-lat data instead of utm, hence, let's
# reproject start and endcoordinates to lon-lat
xy_from <- reprojCoords(
    xy   = cbind(ssf$x1_, ssf$y1_)
  , from = "+init=epsg:32734"
  , to   = "+init=epsg:4326"
)
xy_to <- reprojCoords(
    xy   = cbind(ssf$x2_, ssf$y2_)
  , from = "+init=epsg:32734"
  , to   = "+init=epsg:4326"
)
ssf[, c("x1_", "y1_")] <- xy_from
ssf[, c("x2_", "y2_")] <- xy_to
ssf[, c("x", "y")] <- NULL

# Reorder the columns a bit
ssf <- dplyr::select(ssf, DogName, Sex, BurstID, step_id_, everything(), -TimestampRounded)

# Store to file
write_csv(ssf, "03_Data/02_CleanData/00_General_Dispersers_SSF.csv")

################################################################################
#### Session Information
################################################################################
# Create folder into which the session info goes to
dir.create("02_R-Scripts/99_SessionInformation/02_Analysis", showWarnings = F)

# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/02_Analysis/00_SSF.rds")
cat("Done :)\n")
