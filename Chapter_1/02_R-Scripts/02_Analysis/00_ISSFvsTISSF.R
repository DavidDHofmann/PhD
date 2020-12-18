################################################################################
#### Step Selection Function - Generation of Random Steps
################################################################################
# Description: In this script we coerce our gps data to steps and generate
# random steps for the integratest step selection analysis

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
library(cowplot)      # For nice plots
library(ggpubr)       # To add summary stats to ggplot

# Set seed for reproducability
set.seed(1234)

# Load the gps data
data <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL.csv")

# Let's check the total number of dispersers
sub <- subset(data, State == "Disperser")
subset(sub, DogName == "Appalachia")
subset(sub, DogName == "Ripley")
dispersers <- unique(sub$DogName)
length(dispersers)

################################################################################
#### Data Cleaning
################################################################################
# Loop over the two methods and prepare the data
methods <- c("iSSF", "TiSSF")
tracks <- lapply(methods, function(x){

  # Select method
  method <- x

  # Load the gps data
  data <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL.csv")

  # We only need dispersers' data
  data <- subset(data, State == "Disperser")

  # Let's create a timestamp that is rounded to the nearest hour
  data$TimestampRounded <- round_date(data$Timestamp, "1 hour")

  # Make sure there are no duplicates!
  table(duplicated(data[, c("DogName", "TimestampRounded")]))
  table(hour(data$TimestampRounded))

  # Count number of dogs
  length(unique(data$DogName))

  # Resample fixes to 2 hours
  data <- data %>% group_by(DogName) %>% nest()

  # Resample
  data$data <- suppressMessages(
    pbmclapply(1:nrow(data)
      , ignore.interactive = T
      , mc.cores           = detectCores() - 1
      , FUN                = function(x){
        resFix2(data$data[[x]], hours = 2, start = 1, tol = 0.5)
      }
    )
  )
  data <- unnest(data)

  # We either use iSSF or TiSSF. Hence, we may need to make the data regular.
  # For this, we will only be using fixes that were collected during the
  # following hours: 3, 7, 15, 19, 23
  if (method == "iSSF"){
    data <- subset(data, hour(TimestampRounded) %in% c(3, 7, 15, 19, 23))
  }
  else {
    data <- subset(data, hour(TimestampRounded) %in% seq(1, 23, 2))
  }

  # We want to convert the fixes to steps. However, we don't want to consider
  # any steps that take longer than 8.25 hours. Let's therefore identify
  # consecutive bursts with fixes spread < 8.25 for each individual
  data <- data %>%
    group_by(DogName) %>%
    mutate(dt_ = Timestamp - lag(Timestamp)) %>%
    mutate(dt_ = as.numeric(dt_, units = "hours")) %>%
    ungroup()

  # A new burst starts whenever a step takes longer than 8.25 hours or if a step
  # takes longer than 4.25 hours unless its realized at 15:00 (actually lagged)
  if (method == "iSSF"){
    data <- data %>%
      mutate(NewBurst = ifelse(dt_ > 8.25 | (dt_ > 4.25 & hour(TimestampRounded) != 15) | is.na(dt_), yes = 1, no = 0)) %>%
      mutate(BurstID = cumsum(NewBurst)) %>%
      dplyr::select(-c(NewBurst))
  } else {
    data <- data %>%
      mutate(NewBurst = ifelse(dt_ > 8.25 | is.na(dt_), yes = 1, no = 0)) %>%
      mutate(BurstID = cumsum(NewBurst)) %>%
      dplyr::select(-c(NewBurst))
  }

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

  # Add a row indicating the method
  tracks$method <- method

  # Return the tracks
  return(tracks)
}) %>% do.call(rbind, .)

# Add some further columns and remove undesired columns
tracks <- tracks %>%
  mutate(
      speed_      = sl_ / dt_
    , dt_rounded  = round(dt_)
    , inactive    = hour(
        round_date(t1_, "1 hour")) >= 7 &
        hour(round_date(t2_, "1 hour")) <= 15 &
        hour(round_date(tracks$t1_, "1 hour")) <
        hour(round_date(tracks$t2_, "1 hour"))
  ) %>%
  dplyr::select(-c(id, x1_, x2_, y1_, y2_, State))

# We can also calculate an adjusted speed for the time between 7 and 15
tracks$dt_adjusted <- ifelse(tracks$inactive, yes = tracks$dt_ / 2, tracks$dt_)
tracks$speed_adjusted <- tracks$sl_ / tracks$dt_adjusted

# Check out steptimes
table(tracks$dt_rounded)
subset(tracks, dt_rounded == 0)

################################################################################
#### Descriptive Stats
################################################################################
# Compare number of steps
tracks %>%
  group_by(method) %>%
  summarize(NoFixes = n()) %>%
  mutate(difference = NoFixes - lag(NoFixes))

# Compare number of individuals
tracks %>%
  group_by(method) %>%
  summarize(NoDogs = length(unique(DogName)))

################################################################################
#### Comparison of step speeds over all animals
################################################################################
# Check out the distribution of speeds (true speeds)
p1 <- tracks %>%
  ggplot(aes(x = speed_, col = method, fill = method)) +
  geom_histogram(alpha = 0.5) +
  ggtitle("True Speeds")

# Check out the distribution of speeds (true adjusted speeds)
p2 <- tracks %>%
  ggplot(aes(x = speed_adjusted, col = method, fill = method)) +
  geom_histogram(alpha = 0.5) +
  ggtitle("Adjusted Speeds")

# Compare distribution of true and adjusted speeds
ggarrange(p1, p2, nrow = 2)

################################################################################
#### Comparison of step speeds per animal
################################################################################
# Check out the distribution of true speeds for each animal
p1 <- tracks %>%
  ggplot(aes(x = speed_, col = method, fill = method)) +
  geom_histogram(alpha = 0.5) +
  facet_wrap("DogName") +
  ggtitle("True Speeds")

# Check out the distribution of adjusted speeds for each animal
p2 <- tracks %>%
  ggplot(aes(x = speed_adjusted, col = method, fill = method)) +
  geom_histogram(alpha = 0.5) +
  facet_wrap("DogName") +
  ggtitle("Adjusted Speeds")

# Compare distribution of true and adjusted speeds
ggarrange(p1, p2, nrow = 2)

################################################################################
#### Comparison of step speeds depending on the step duration
################################################################################
# Compare how the true step speeds change with the step duration
p1 <- tracks %>%
  group_by(method, dt_rounded, inactive) %>%
  summarize(speed_ = mean(speed_)) %>%
  ggplot(aes(x = dt_rounded, y = speed_, col = inactive)) +
  geom_jitter(data = tracks, size = 0.5, position = position_dodge(0.5)) +
  geom_line() +
  facet_wrap("method") +
  ggtitle("True Speeds")

# Compare how the true step speeds change with the step duration
p2 <- tracks %>%
  group_by(method, dt_rounded, inactive) %>%
  summarize(speed_adjusted = mean(speed_adjusted)) %>%
  ggplot(aes(x = dt_rounded, y = speed_adjusted, col = inactive)) +
  geom_jitter(data = tracks, size = 0.5, position = position_dodge(0.5)) +
  geom_line() +
  facet_wrap("method") +
  ggtitle("Adjusted Speeds")

# Compare distribution of true and adjusted speeds
ggarrange(p1, p2, nrow = 2)

################################################################################
#### Fit model to step speed
################################################################################
# Try to fit a model to the step speed (only for TiSSF)
library(lme4)
library(lmerTest)
sub <- subset(tracks, method == "TiSSF")
mod1 <- lm(speed_ ~ dt_, data = sub)
mod2 <- lmer(speed_ ~ dt_ + (1|DogName), data = sub)
mod3 <- lmer(speed_ ~ dt_ + inactive + (1|DogName), data = sub)
mod4 <- lmer(speed_ ~ dt_ + I(dt_**2) + inactive + (1|DogName), data = sub)
mod4 <- lmer(speed_ ~ dt_ + I(dt_**2) + inactive + (1|DogName), data = sub)
summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)

# Predict on the data gain
newdat <- expand_grid(
    DogName = unique(tracks$DogName)
  , dt_ = seq(0, 8, 0.5)
  , inactive = c(T, F)
)
newdat$speed_ <- predict(mod4, newdat)

# Visualize
newdat %>%
  ggplot(aes(x = dt_, y = speed_, col = inactive)) +
  geom_line() +
  geom_jitter(data = sub, width = 0.3, size = 0.4) +
  facet_wrap("DogName")

################################################################################
#### Gamma Distribution Comparison
################################################################################
# Let's fit a gamma distribution for each of the scenarios
sub <- subset(tracks, method == "TiSSF")
comparison <- expand_grid(
    steptime = unique(sub$dt_rounded)
  , inactive = c(T, F)
  , shape    = NA
  , scale    = NA
)

# Loop through them
for (i in 1:nrow(comparison)){
  sl_ <- sub$sl_[sub$dt_rounded == comparison$steptime[i] & sub$inactive == comparison$inactive[i]]
  if (length(sl_) > 1){
    fit <- fit_distr(sl_, dist_name = "gamma")
    comparison$shape[i] <- fit$params$shape
    comparison$scale[i] <- fit$params$scale
  } else {
    comparison$shape[i] <- NA
    comparison$scale[i] <- NA
  }
}

# Remove NA
comparison <- na.omit(comparison)

# Visualize shape
comparison %>%
  ggplot(aes(x = steptime, y = shape, col = inactive)) +
  geom_point()

# Visualize scale
comparison %>%
  ggplot(aes(x = steptime, y = scale, col = inactive)) +
  geom_point()
