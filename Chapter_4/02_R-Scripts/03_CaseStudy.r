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
library(amt)           # To fit distributions

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
    Missingness    = seq(0, 0.5, by = 0.1)                  # Fraction of the fixes that is removed
  , Forgiveness    = 1:5                                    # Allowed lag of steps (in steps)
  , Replicate      = 1:100                                  # Number of replicates for each combination
  , Distributions  = c("uncorrected", "naive", "dynamic")   # Which distributions to use
)

# Adjust column names slightly
obs <- rename(obs, ID = DogName)

# Load all covariate layers
water <- stack("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd")
trees <- stack("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/01_LandCover_TreeCover_MODIS.grd")
shrub <- stack("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/01_LandCover_NonTreeVegetation_MODIS.grd")
human <- stack("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluenceBuff_FACEBOOK.grd")
human <- human[["Buffer_5000"]]

# Extract dates from layernames
dates_water <- ymd(substr(names(water), start = 2, stop = 12))
dates_trees <- ymd(substr(names(trees), start = 2, stop = 12))
dates_shrub <- ymd(substr(names(shrub), start = 2, stop = 12))

# All dates should be the same, so we can safely keep only one of the objects
# to reduce redundancy
if (all.equal(dates_water, dates_trees) & all.equal(dates_trees, dates_shrub)) {
    dates <- dates_water
    rm(dates_water, dates_trees, dates_shrub)
  } else {
    stop("Dates don't match")
}

# Prepare a list that stores a ppp layer for water (Code 1) for each floodmap
water_ppp <- suppressMessages(
  pbmclapply(
      X                   = 1:nlayers(water)
    , mc.cores            = detectCores() / 2
    , ignore.interactive  = T
    , FUN                 = function(x) {
      points <- rasterToPoints(water[[x]], fun = function(z) {z == 1}, spatial = T)
      points <- spTransform(points, CRS("+init=epsg:32734"))
      points <- as(points, "ppp")
      return(points)
      gc()
  })
)

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
computeSSF <- function(data, n_rsteps, distributions) {

  # Indicate case steps
  data$case <- 1

  # Cannot work with steps that have no turning angle
  data <- subset(data, !is.na(relta))

  # Create a new dataframe for alternative steps
  rand <- data[rep(1:nrow(data), each = n_rsteps), ]

  # Indicate that they are control steps (case = 0)
  rand$case <- 0

  # Sample new step lengths and turning angles according to the specified
  # distributions
  if (distributions == "uncorrected") {
    rand$sl <- rgamma(n = nrow(rand)
      , scale = dists$uncorrected$sl$scale
      , shape = dists$uncorrected$sl$shape
    )
    rand$relta_new <- rvonmises(n = nrow(rand)
      , kappa = dists$uncorrected$ta$kappa
      , mu    = dists$uncorrected$ta$mu
      , by    = 0.01
    )
  } else if (distributions == "naive") {
    rand$sl <- rgamma(n = nrow(rand)
      , scale = dists$uncorrected$sl$scale * rand$duration
      , shape = dists$uncorrected$sl$shape
    )
    rand$relta_new <- rvonmises(n = nrow(rand)
      , kappa = dists$uncorrected$ta$kappa
      , mu    = dists$uncorrected$ta$mu
      , by    = 0.01
    )
  } else if (distributions == "dynamic") {
    rand$sl <- rgamma(n = nrow(rand)
      , scale = dists$dynamic$sl$scale[match(rand$duration, dists$dynamic$sl$duration)]
      , shape = dists$dynamic$sl$shape[match(rand$duration, dists$dynamic$sl$duration)]
    )
    rand$relta_new <- sapply(1:nrow(rand), function(z) {
      rvonmises(n = 1
        , kappa = dists$dynamic$ta$kappa[dists$dynamic$ta$duration == rand$duration[z]]
        , mu    = dists$dynamic$ta$mu[dists$dynamic$ta$duration == rand$duration[z]]
        , by    = 0.01
      )
    })
  } else {
    stop("Provide valid input for the desired distributions")
  }

  # Calculate new "absolute" turning angle
  rand$relta_diff <- rand$relta - rand$relta_new
  rand$absta <- rand$absta - rand$relta_diff
  rand$absta[rand$absta > 2 * pi] <- rand$absta[rand$absta > 2 * pi] - 2 * pi
  rand$absta[rand$absta < 0 * pi] <- rand$absta[rand$absta < 0 * pi] + 2 * pi
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

      # Interpolate between points
      ints <- interpolatePoints(
          x1 = data$x[x]
        , x2 = data$x_to[x]
        , y1 = data$y[x]
        , y2 = data$y_to[x]
        , by = 250
      )

      # Generate spatial points
      ints <- SpatialPoints(ints)
      crs(ints) <- "+init=epsg:32734"

      # Prepare interpolated points as ppp object for later
      ppp <- as(ints, "ppp")

      # Reproject interpolated points to lonlat
      ints <- spTransform(ints, "+init=epsg:4326")
      ints <- coordinates(ints)

      # Determine the index of the layer closest in date
      index <- which.min(abs(as.Date(data$Timestamp[x]) - dates))[1]

      # Run the extraction on the respective covariate layers
      extr <- data.frame(
          Water           = mean(raster::extract(water[[index]], ints))
        , Trees           = mean(raster::extract(trees[[index]], ints))
        , Shrubs          = mean(raster::extract(shrub[[index]], ints))
        , HumanInfluence  = mean(raster::extract(human, ints))
        , DistanceToWater = mean(nncross(ppp, water_ppp[[index]])$dist)
      )
      return(extr)
    })

  } else {
    extracted <- lapply(1:nrow(data), function(x) {

      # Interpolate between points
      ints <- interpolatePoints(
          x1 = data$x[x]
        , x2 = data$x_to[x]
        , y1 = data$y[x]
        , y2 = data$y_to[x]
        , by = 250
      )

      # Generate spatial points
      ints <- SpatialPoints(ints)
      crs(ints) <- "+init=epsg:32734"

      # Prepare interpolated points as ppp object for later
      ppp <- as(ints, "ppp")

      # Reproject interpolated points to lonlat
      ints <- spTransform(ints, "+init=epsg:4326")
      ints <- coordinates(ints)

      # Determine the index of the layer closest in date
      index <- which.min(abs(as.Date(data$Timestamp[x]) - dates))[1]

      # Run the extraction on the respective covariate layers
      extr <- data.frame(
          Water           = mean(raster::extract(water[[index]], ints))
        , Trees           = mean(raster::extract(trees[[index]], ints))
        , Shrubs          = mean(raster::extract(shrub[[index]], ints))
        , HumanInfluence  = mean(raster::extract(human, ints))
        , DistanceToWater = mean(nncross(ppp, water_ppp[[index]])$dist)
      )
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

# Fit a "base" step length distribution that we can use for the uncorrected case
fitted <- obs %>%
  rarifyData(missingness = 0) %>%
  computeBursts(forgiveness = 1) %>%
  computeMetrics() %>%
  mutate(sl = ifelse(sl == 0, 1, sl)) %>%
  pull(sl) %>%
  fit_distr(dist_name = "gamma")

# Put the parameters into a specialized list
sl_dist <- list(shape = fitted$params$shape, scale = fitted$params$scale)   # According to Avgar et al. 2016
ta_dist <- list(kappa = 0, mu = 0)                                          # According to Avgar et al. 2016
dists <- list(uncorrected = list(sl = sl_dist, ta = ta_dist))

# Try out the functions to see how they work
testing <- rarifyData(obs[1:1000, ], missingness = 0.5)
testing <- computeBursts(testing, forgiveness = 2)
testing <- computeMetrics(testing)
testing <- computeSSF(testing, n_rsteps = 10, distributions = "uncorrected")
testing <- computeCovars(testing, covars, multicore = T)
testing <- runModel(testing)
testing

################################################################################
#### Fit Step Length Distributions to Different Step Durations
################################################################################
# Parametrize separate step length distributions and turning angle distributions
# for the different accepted step durations. Replicate 1000 times.
dists_dynamic <- pbmclapply(
    X                  = 1:1000
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {
    params <- obs %>%
      rarifyData(missingness = 0.5) %>%
      computeBursts(forgiveness = max(dat$Forgiveness)) %>%
      computeMetrics() %>%
      mutate(sl = ifelse(sl == 0, 1, sl)) %>%
      dplyr::select(duration, sl, relta) %>%
      subset(duration <= max(dat$Forgiveness * 4)) %>%
      group_by(duration) %>%
      nest() %>%
      mutate(DistParams = map(data, function(y) {
        fit_sl <- fit_distr(y$sl, "gamma")
        fit_ta <- fit_distr(y$relta, "vonmises")
        fit <- data.frame(
            shape = fit_sl$params$shape
          , scale = fit_sl$params$scale
          , kappa = fit_ta$params$kappa
          , mu    = fit_ta$params$mu
        )
        return(fit)
      })) %>%
      dplyr::select(duration, DistParams) %>%
      unnest(DistParams) %>%
      arrange(duration)
  return(params)
}) %>% do.call(rbind, .)

# Summarize values
dists_means <- dists_dynamic %>%
  pivot_longer(shape:mu, names_to = "Parameter", values_to = "Value") %>%
  group_by(duration, Parameter) %>%
  summarize(mean = mean(Value), sd = sd(Value), .groups = "drop")

# Put values together with the uncorrected distributions
dists$dynamic <- dists_means %>%
  dplyr::select(duration, Parameter, mean) %>%
  pivot_wider(names_from = Parameter, values_from = mean) %>%
  list(
      sl = select(., duration, shape, scale)
    , ta = select(., duration, kappa, mu)
  )

# Visualize everything
dists_dynamic %>%
  pivot_longer(shape:mu, names_to = "Parameter", values_to = "Value") %>%
  ggplot(aes(x = duration, y = Value)) +
    geom_jitter(width = 0.1, alpha = 0.2, size = 0.5) +
    geom_point(data = dists_means, aes(x = duration, y = mean), col = "orange", size = 5) +
    geom_line(data = dists_means, aes(x = duration, y = mean), col = "orange") +
    facet_wrap(~ Parameter, nrow = 2, scales = "free") +
    theme_minimal() +
    xlab("Step Duration") +
    ylab("Parameter Estimate")
