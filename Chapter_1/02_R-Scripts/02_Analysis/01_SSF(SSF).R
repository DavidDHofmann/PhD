################################################################################
#### Step Selection Function - Generation of Random Steps
################################################################################
# Description: In this script we use our gps data to create steps and propose
# potential alternative steps

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load packages
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(lubridate)    # To handle dates
library(amt)          # To coerce gps fixes to steps
library(raster)       # To handle spatial data
library(rgdal)        # To handle spatial data
library(pbmcapply)    # For multicore abilities

# Specify if you're going for iSSF or TiSSF
method <- "iSSF"
# method <- "TiSSF"

# Set seed for reproducability
set.seed(1234)

################################################################################
#### Data Cleaning
################################################################################
# Load the gps data
data <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL.csv")

# We only need dispersers' data
data <- subset(data, State == "Disperser")

# Let's create a timestamp that is rounded to the nearest hour
data$TimestampRounded <- round_date(data$Timestamp, "1 hour")

# Make sure there are no duplicates!
table(duplicated(data[, c("DogName", "TimestampRounded")]))
table(hour(data$TimestampRounded))

# We either use iSSF or TiSSF. Hence, we may need to make the data regular. For
# this, we will only be using fixes that were collected during the following
# hours: 3, 7, 15, 19, 23
if (method == "iSSF"){
  data <- subset(data, hour(TimestampRounded) %in% c(3, 7, 15, 19, 23))
}

# We want to convert the fixes to steps. However, we don't want to consider any
# steps that take longer than 8.25 hours. Let's therefore identify consecutive
# bursts with fixes spread < 8.25 for each individual
data <- data %>%
  group_by(DogName) %>%
  mutate(dt_ = Timestamp - lag(Timestamp)) %>%
  mutate(dt_ = as.numeric(dt_, units = "hours")) %>%
  ungroup()

# A new burst starts whenever a step takes longer than 8.25 hours
data <- data %>%
  mutate(NewBurst = ifelse(dt_ > 8.25 | is.na(dt_), yes = 1, no = 0)) %>%
  mutate(BurstID = cumsum(NewBurst)) %>%
  dplyr::select(-c(NewBurst))

# We can only work with bursts that contain at least three fixes (since we will
# need to calculate turning angles)
data <- data %>%
  group_by(DogName, BurstID) %>%
  nest() %>%
  mutate(Nrow = map(data, nrow) %>% do.call(rbind, .)) %>%
  subset(Nrow >= 3) %>%
  unnest(data) %>%
  ungroup() %>%
  dplyr::select(-c(Nrow, dt_))

# Create bursts that are unique
data$BurstID <- data %>% group_indices(DogName, BurstID)

# Now we can coerce the data to proper steps. Note that we pretend that the
# burst is the animal. This allows us to only calculate turning angles and step
# lengths between consecutive steps.
tracks <- data %>%
  make_track(.
    , .x      = x
    , .y      = y
    , .t      = Timestamp
    , id      = BurstID
    , crs     = CRS("+init=epsg:4326")
    , State   = State
    , DogName = DogName
  ) %>%

  # Transform the tracks to utm
  transform_coords(CRS("+init=epsg:32734")) %>%

  # Nest the tracks
  nest(-"id") %>%

  # Turn to a step representation (look up the amt vignette for details)
  mutate(data = map(data, function(x){
    x %>%

      # The function creates steps from the fixes in each tibble row. The option
      # "keep_cols" allows to keep the time of day (tod) column that we added
      steps(keep_cols = "start") %>%

      # Transform the difftime column to a numeric column to avoid that there
      # are heterogeneous units
      mutate(., dt_ = as.numeric(dt_, units = "hours"))
  })) %>%

  # Unnest the tibble
  unnest(data) %>%

  # Multiply turning angles with negative one (for some reason the package
  # calculates the turning angles conuterclockwise but we want them clockwise)
  mutate(ta_ = ta_ * (-1)) %>%

  # Add a column indicating the absolute turning angle (important for the
  # ellipses). We can use the function we created above.
  mutate(absta_ = absAngle(.)) %>%

  # We can only work with steps for which we have a turning angle. Let's get rid
  # of any steps where the turning angle is NA
  filter(!is.na(ta_))

# Generate step id and indicate if case or control step
tracks <- tracks %>% mutate(
    step_id_ = 1:nrow(.)
  , case_    = 1
)

################################################################################
#### Generation of Random Steps
################################################################################
# We want to normalize the step length to 4 hours (for TiSSF). However, remember
# that we assume steps between 7 and 15 oclock to be half as quick as the other
# steps. Let's thus identify the steps that were realized between 7 and 15
# o'clock
index <- hour(round_date(tracks$t1_, "1 hour")) >= 7 &
  hour(round_date(tracks$t2_, "1 hour")) <= 15 &
  hour(round_date(tracks$t1_, "1 hour")) < hour(round_date(tracks$t2_, "1 hour"))

# Create a corrected step duration (we assume that time runs half as fast from 7
# to 15 oclock)
tracks$dtcorr_ <- tracks$dt_
tracks$dtcorr_[index] <- tracks$dt_[index] / 2

# Normalize step lengths to 4 hours (only important for TiSSF)
tracks$normsl_ <- tracks$sl_ / tracks$dtcorr_ * 4

# To fit a gamma distribution, we need to make all steps larger than zero. So
# let's add one meter to all steps that are currently zero (only four steps)
tracks$sl_[tracks$sl_ == 0] <- 1
tracks$normsl_[tracks$normsl_ == 0] <- 1

# Fit a Gamma distribution to the observed step lengths
if (method == "TiSSF"){
  sl <- fit_distr(tracks$normsl_, "gamma")
} else {
  sl <- fit_distr(tracks$sl_, "gamma")
}

# Also calculate the rate parameter
sl$params$rate <- 1 / sl$params$scale

# Let's visualize the fit
x <- seq(0, max(tracks$sl_), 1)
y <- dgamma(x, rate = sl$params$rate, shape = sl$params$shape)
# hist(tracks$normsl_, freq = F, breaks = 100)
hist(tracks$sl_, freq = F, breaks = 100)
lines(y ~ x, col = "red")

# Write the distribution to file
if (method == "TiSSF"){
  write_rds(sl, "03_Data/03_Results/99_GammaDistribution(TiSSF).rds")
} else {
  write_rds(sl, "03_Data/03_Results/99_GammaDistribution(iSSF).rds")
}

# Define the number of random steps
nsteps <- 24

# Now we use the fitted gamma distribution to sample step lengths. For the
# turning angles we follow Avgar et al. 2016 and use a uniform distribution
randomSteps <- pbmclapply(
    X                  = 1:nrow(tracks)
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(x){

  # Draw random turning angles
  ta_new <- runif(nsteps, min = -pi, max = pi)

  # Draw random step lengths with the fitted parameters. In case we want to
  # apply TiSSF, we need to adjust the rate parameter of the distribution
  # accordingly
  sl_new <- rgamma(nsteps
    , shape = sl$params$shape
    , rate  = sl$params$rate
  )

  # In case we go for the TiSSF, we need to adjust the step length according to
  # the (corrected) step time
  if (method == "TiSSF"){
    sl_new <- sl_new * tracks$dtcorr_[x] / 4
  }

  # Put the step lengths and turning angles into a new dataframe and calculate
  # the new endpoints of each random step. We also indicate that the steps are
  # control steps (i.e. 0)
  randomSteps <- tracks[rep(x, nsteps), ] %>%
    mutate(.
      , absta_  = absta_ + (ta_new - ta_)
      , ta_     = ta_new
      , sl_     = sl_new
      , case_   = 0
      , x2_     = x1_ + sin(absta_) * sl_
      , y2_     = y1_ + cos(absta_) * sl_
    )

  # Return the sampled steps
  return(randomSteps)

}) %>% do.call(rbind, .)

# We need to make sure that the absolute turning angle ranges from 0 to 2 * pi
randomSteps$absta_[randomSteps$absta_ > 2 * pi] <-
  randomSteps$absta_[randomSteps$absta_ > 2 * pi] - 2 * pi
randomSteps$absta_[randomSteps$absta_ < 0] <-
  randomSteps$absta_[randomSteps$absta_ < 0] + 2 * pi

# Merge the dataframes of the observed and random steps
ssf <- rbind(tracks, randomSteps) %>%
  arrange(step_id_, desc(case_)) %>%
  transform(case_ = as.logical(case_))

# Coerce the steps to spatial lines. We can use the function we defined earlier
# for this
lines <- lineTrack(ssf, CRS("+init=epsg:32734"))

# Transform the lines to WGS84
lines <- spTransform(lines, CRS("+init=epsg:4326"))

# Coerce the duration column to a format that can be stored (numeric rather than
# difftime)
lines$dt_ <- as.numeric(lines$dt_)

# Store the object to a shapefile
if (method == "iSSF"){
  filename <- "00_General_Dispersers_Popecol(iSSF)"
} else {
  filename <- "00_General_Dispersers_Popecol(TiSSF)"
}

# Store the stuff into both output directories
writeOGR(lines
  , dsn       = "03_Data/02_CleanData"
  , layer     = filename
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
