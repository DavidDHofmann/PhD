################################################################################
#### Integrated Step Selection Analysis of African Wild Dog Data
################################################################################
# Description: Use ISSF analysis estimate preferences on the simulated data

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(survival)      # To run conditional logistic regression
library(sf)            # For plotting
library(davidoff)      # Custom functions
library(lubridate)     # To handle dates

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/00_Functions.r")

# Suppress scientific notation
options(scipen = 999)

# Load observed movement data of dispersing wild dogs
obs <- "03_Data/01_RawData/WildDogs.csv" %>%
  read_csv() %>%
  subset(State == "Disperser") %>%
  mutate(TimestampRounded = round_date(Timestamp, "1 hour")) %>%
  group_by(DogName) %>%
  nest()

# Resample data to 1 hour
obs$data <- suppressMessages(
  pbmclapply(1:nrow(obs)
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x){
      resFix(obs$data[[x]], hours = 1, start = 1, tol = 0.5)
    }
  )
)

# Unnest and keep only fixes on the attempted schedule
obs <- obs %>%
  unnest(data) %>%
  subset(hour(TimestampRounded) %in% c(3, 7, 15, 19, 23))

# Determine "step_number" and "step_ids"
obs$step_id <- 1:nrow(obs)
obs <- obs %>% mutate(step_number = row_number())

# Reproject data to utm
coords <- obs[, c("x", "y")]
coordinates(coords) <- c("x", "y")
crs(coords) <- "+init=epsg:4326"
coords <- spTransform(coords, "+init=epsg:32734")
obs$x <- coordinates(coords)[, c("x")]
obs$y <- coordinates(coords)[, c("y")]

# # Check for missing fixes
# test <- data %>%
#   group_by(DogName) %>%
#   nest() %>%
#   mutate(Schedule = map(data, function(x) {
#     range <- range(x$TimestampRounded)
#     attempted <- tibble(Attempted = seq(range[1], range[2], by = "4 hour"))
#     attempted <- left_join(attempted, x, by = c("Attempted" = "TimestampRounded"))
#     attempted <- subset(attempted, hour(Attempted) != 11)  # At 11 we neve expect a fix
#     return(attempted)
#   }))
# testi <- dplyr::select(test, Schedule) %>% unnest(Schedule)

# In some individuals the missingness is not unplanned but due to either a
# different schedule or because there was actually no collar on the animal
# during that time.


# Specify the different design combinations through which we want to run. Note
# that a forgiveness of 1 refers to a regular step selection function
dat <- expand_grid(
    Missingness = seq(0, 0.8, by = 0.1)  # Fraction of the fixes that is removed
  , Forgiveness = 1:5                    # Allowed lag of steps (in steps)
  , Replicate   = 1:100                  # Number of replicates for each combination
  , AdjustDists = c(T, F)                # Whether to use a dynamic distribution for step lengths
)

# Adjust column names slightly
obs <- rename(obs, ID = DogName)

################################################################################
#### Functions
################################################################################
# Function to create a dataset with rarified observations
rarifyData <- function(data, missingness) {
  rarified <- data[sort(sample(1:nrow(data), size = nrow(data) * (1 - missingness))), ]
  return(rarified)
}

# Function to compute bursts (per ID, depending on the fogriveness)
computeBursts <- function(data, forgiveness) {

  # Nest data by id
  data_bursted <- data %>%
    group_by(ID) %>%
    nest() %>%

    # Can only work with individuals where there is more than 1 fix
    mutate(NFixes = map(data, nrow)) %>%
    subset(NFixes > 1) %>%
    select(-NFixes) %>%

    # Compute bursts by id. A new burst is defined if the duration is >
    # forgiveness
    mutate(data = map(data, function(x) {
      x$duration <- difftime(lead(x$TimestampRounded), x$TimestampRounded, units = "hours")
      x$duration[hour(x$TimestampRounded) == 7] <- x$duration[hour(x$TimestampRounded) == 7] - 4
      x$irregular <- x$duration > forgiveness * 4 # To account for the fact that we have four-hourly steps
      x$burst <- NA
      x$burst[1] <- 1
      x$burst[2:nrow(x)] <- lag(cumsum(x$irregular) + 1)[-1]
    return(x)

  # Unnest the data again
  })) %>% unnest(data)

  # Return the "bursted" data
  return(data_bursted)
}

# Function to compute step metrics per burst
computeMetrics <- function(data) {

  # Nest by ID and burst
  data_metrics <- data %>%
    group_by(ID, burst) %>%
    nest() %>%

    # Can only work with bursts longer than 1 fix
    mutate(NFixes = map(data, nrow)) %>%
    subset(NFixes > 1) %>%
    select(-NFixes) %>%

    # Compute step metrics
    mutate(data = map(data, function(z) {
      metrics <- stepMet(x = z$x, y = z$y)
      metrics <- cbind(z, metrics)
      return(metrics)
    })) %>%

    # Unnest again and tidy up
    unnest(data) %>%
    dplyr::select(burst, step_number, step_id, everything()) %>%
    ungroup()

  # Return the data containing step metrics
  return(data_metrics)
}

# Function to generate random steps
computeSSF <- function(data, n_rsteps, adjust_dists) {

  # Indicate case steps
  data$case <- 1

  # Cannot work with steps that have no turning angle
  data <- subset(data, !is.na(relta))

  # Create a new dataframe for alternative steps
  rand <- data[rep(1:nrow(data), each = n_rsteps), ]

  # Indicate that they are control steps (case = 0)
  rand$case <- 0

  # Sample new step lengths and turning angles
  if (adjust_dists) {
    rand$sl <- sapply(1:nrow(rand), function(z) {
      rgamma(n = 1
        , scale = dists_means$mean[dists_means$duration == rand$duration[z] & dists_means$Parameter == "Scale"]
        , shape = dists_means$mean[dists_means$duration == rand$duration[z] & dists_means$Parameter == "Shape"]
      )
    })
    rand$relta_new <- sapply(1:nrow(rand), function(z) {
      rvonmises(n = 1
        , kappa = dists_means$mean[dists_means$duration == rand$duration[z] & dists_means$Parameter == "Kappa"]
        , mu    = dists_means$mean[dists_means$duration == rand$duration[z] & dists_means$Parameter == "Mu"]
        , by    = 0.01
      )
    })
  } else {
    rand$sl <- rgamma(n = nrow(rand)
      , scale = sl_dist["scale"] * rand$duration
      , shape = sl_dist["shape"]
    )
    rand$relta_new <- rvonmises(n = nrow(rand)
      , kappa = ta_dist["kappa"]
      , mu    = ta_dist["mu"]
      , by    = 0.01
    )
  }

  # Calculate new "absolute" turning angle
  rand$relta_diff <- rand$relta - rand$relta_new
  rand$absta <- rand$absta - rand$relta_diff
  rand$absta[rand$absta > 2 * pi] <-
    rand$absta[rand$absta > 2 * pi] - 2 * pi
  rand$absta[rand$absta < 0] <-
    rand$absta[rand$absta < 0] + 2 * pi
  rand$relta <- rand$relta_new

  # Remove undesired stuff
  rand$relta_new <- NULL
  rand$relta_diff <- NULL

  # Put steps together
  all <- rbind(data, rand)
  all <- arrange(all, ID, step_id, desc(case))

  # Calculate new endpoints
  all$x_to <- all$x + sin(all$absta) * all$sl
  all$y_to <- all$y + cos(all$absta) * all$sl

  # Sort and ungroup
  all <- dplyr::select(
      all
    , ID
    , step_number
    , step_id
    , x
    , y
    , x_to
    , y_to
    , everything()
  )
  all <- ungroup(all)
  return(all)
}

# Function to extract covariates along steps and compute step covariates
computeCovars <- function(data, covariates, multicore = F) {

  # Run the covariate extraction
  if (multicore) {
    extracted <- pbmclapply(1:nrow(data), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
      ints <- interpolatePoints(
          x1 = data$x[x]
        , x2 = data$x_to[x]
        , y1 = data$y[x]
        , y2 = data$y_to[x]
        , by = 1
      )
      extr <- raster::extract(covariates, ints)
      extr <- colMeans(extr)
      return(extr)
    })

  } else {
    extracted <- lapply(1:nrow(data), function(x) {
      ints <- interpolatePoints(
          x1 = data$x[x]
        , x2 = data$x_to[x]
        , y1 = data$y[x]
        , y2 = data$y_to[x]
        , by = 1
      )
      extr <- raster::extract(covariates, ints)
      extr <- colMeans(extr)
      return(extr)
    })
  }
  extracted <- as.data.frame(do.call(rbind, extracted))

  # Bind with other data
  data <- cbind(data, extracted)

  # Calculate movement metrics
  data$log_sl <- log(data$sl)
  data$cos_ta <- cos(data$relta)

  # Return the results
  return(data)
}

# Function to run the step selection model
runModel <- function(data) {

  # Run the step selection model
  mod <- clogit(case ~
    + sl
    + log_sl
    + cos_ta
    + water
    + elev
    + dist
    + strata(step_id)
    , data = data
  )

  # Extract model coefficients
  ci <- confint(mod, level = 0.95)
  coefs <- summary(mod)$coefficients
  coefs <- data.frame(
      Coefficient = rownames(coefs)
    , Estimate    = coefs[, "coef"]
    , SE          = coefs[, "se(coef)"]
    , Z           = coefs[, "z"]
    , LCI         = ci[, 1]
    , UCI         = ci[, 2]
  )
  rownames(coefs) <- NULL

  # Return them
  return(coefs)
}

# Try out the functions to see how they work
testing <- rarifyData(obs, missingness = 0.5)
testing <- computeBursts(testing, forgiveness = 2)
testing <- computeMetrics(testing)
testing <- computeSSF(testing, n_rsteps = 10, adjust_dists = F)
testing <- computeCovars(testing, covars, multicore = T)
testing <- runModel(testing)
testing

################################################################################
#### Fit Step Length Distributions to Step Durations
################################################################################