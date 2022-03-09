################################################################################
#### Integrated Step Selection Analysis of African Wild Dog Data
################################################################################
# Description: Use ISSF analysis estimate preferences on the simulated data

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(terra)         # To handle spatial data
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(survival)      # To run conditional logistic regression
library(lubridate)     # To handle dates
library(amt)           # To fit distributions
library(maptools)      # Access to spatstat functions
library(spatstat)      # Access to spatstat functions
library(glmmTMB)       # To fit mixed effects clogit

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/00_Functions.r")

# Suppress scientific notation
options(scipen = 999)

################################################################################
#### GPS Data Cleaning
################################################################################
# Load observed movement data of dispersing wild dogs and nest data by ID
obs <- "03_Data/01_RawData/WildDogs.csv" %>%
  read_csv() %>%
  subset(State == "Disperser") %>%
  mutate(TimestampRounded = round_date(Timestamp, "1 hour")) %>%
  group_by(DogName) %>%
  rename(ID = DogName) %>%
  nest()

# Resample the GPS data to 1 hour
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
  subset(hour(TimestampRounded) %in% c(3, 7, 15, 19, 23)) %>%
  ungroup()

# Determine "step_number" and "step_ids"
obs <- obs %>%
  mutate(step_id = row_number()) %>%
  group_by(ID) %>%
  mutate(step_number = row_number()) %>%
  ungroup()

# Let's define an extent from the GPS locations
ext <- extent(c(min(obs$x), max(obs$x), min(obs$y), max(obs$y))) + 1
ext <- as(ext, "SpatialPolygons")
crs(ext) <- "+init=epsg:4326"

# Reproject data to utm
obs[, c("x", "y")] <- reprojCoords(xy = obs[, c("x", "y")]
  , from = "+init=epsg:4326"
  , to   = "+init=epsg:32734"
)

# Specify the different design combinations through which we want to run. Note
# that a forgiveness of 1 refers to a regular step selection function
dat <- expand_grid(
    Missingness    = seq(0, 0.5, by = 0.1)                                        # Fraction of the fixes that is removed
  , Forgiveness    = 1:5                                                          # Allowed lag of steps (in steps)
  , Replicate      = 1:1                                                          # Number of replicates for each combination
  , Distributions  = c("uncorrected", "naive", "dynamic", "multistep", "model")   # Which distributions to use
)

################################################################################
#### Covariate Layers
################################################################################
# Load all covariate layers
water <- stack("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd")
trees <- stack("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/01_LandCover_TreeCover_MODIS.grd")
shrub <- stack("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/01_LandCover_NonTreeVegetation_MODIS.grd")
human <- stack("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluenceBuff_FACEBOOK.grd")
human <- human[["Buffer_5000"]]

# # Subset for now
# water <- water[[1:2]]
# trees <- trees[[1:2]]
# shrub <- shrub[[1:2]]

# We can crop them to a smaller extent
water <- crop(water, ext)
trees <- crop(trees, ext)
shrub <- crop(shrub, ext)
human <- crop(human, ext)

# Normalize the human influence layer to values between 0 and 1
human <- human / maxValue(human)

# Read all data into memory if possible
water <- readAll(water)
trees <- readAll(trees)
shrub <- readAll(shrub)
human <- readAll(human)

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
      x$duration[hour(x$TimestampRounded) == 7] <- x$duration[hour(x$TimestampRounded) == 7] - 4    # Steps that start at 7 are allowed to be 4 hours longer
      x$irregular <- x$duration > forgiveness * 4                                                   # To account for the fact that we have four-hourly steps
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

# Function to generate random steps. The approach parameter is used to specify
# the approach that should be used to generate random steps
computeSSF <- function(data, n_rsteps, approach) {

  # Generate a new column that indicates that the steps are "observed" steps
  data$case <- 1

  # Cannot work with steps that have no turning angle, so remove them
  data <- subset(data, !is.na(relta))

  # Create a new dataframe into which we can put alternative/random steps
  rand <- data[rep(1:nrow(data), each = n_rsteps), ]

  # Indicate that these steps are random steps (case = 0)
  rand$case <- 0

  ##############################################################################
  #### Approach 1 - Uncorrected
  ##############################################################################
  # Step lengths sampled from "minimal" step-duration distributions
  if (approach == "uncorrected") {
    rand$sl <- rgamma(n = nrow(rand)
      , scale = dists$uncorrected$sl$scale
      , shape = dists$uncorrected$sl$shape
    )
    rand$relta_new <- rvonmises(n = nrow(rand)
      , kappa = dists$uncorrected$ta$kappa
      , mu    = dists$uncorrected$ta$mu
      , by    = 0.01
    )

  ##############################################################################
  #### Approach 2 - Naive
  ##############################################################################
  # Step lengths adjusted merely by multiplying with the step duration
  } else if (approach == "naive") {
    rand$sl <- rgamma(n = nrow(rand)
      , scale = dists$uncorrected$sl$scale * rand$duration
      , shape = dists$uncorrected$sl$shape
    )
    rand$relta_new <- rvonmises(n = nrow(rand)
      , kappa = dists$uncorrected$ta$kappa
      , mu    = dists$uncorrected$ta$mu
      , by    = 0.01
    )

  ##############################################################################
  #### Approach 3 - Dynamic
  ##############################################################################
  # Step lengths and turning angles sampled from distributions fitted to
  # different step-durations
  } else if (approach == "dynamic" | approach == "model") {
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

  ##############################################################################
  #### Approach 4 - Multistep
  ##############################################################################
  # Step lengths and turning angles resulting from multiple random steps, where
  # the number of steps matches the step duration of the observed step
  } else if (approach == "multistep") {
    rand <- rand %>% group_by(duration) %>% nest()
    rand$data <- lapply(1:nrow(rand), function(z) {

      # Extract important info of the current iteration (i.e. the "z")
      duration <- rand$duration[z]
      x        <- rand$data[[z]]$x
      y        <- rand$data[[z]]$y
      relta    <- rand$data[[z]]$relta
      absta    <- rand$data[[z]]$absta

      # Generate sequence of random steps
      for (i in 1:duration) {
        sl <- rgamma(n = nrow(rand$data[[z]])
          , scale = dists$uncorrected$sl$scale
          , shape = dists$uncorrected$sl$shape
        )
        relta_new <- rvonmises(n = nrow(rand$data[[z]])
          , kappa = dists$uncorrected$ta$kappa
          , mu    = dists$uncorrected$ta$mu
          , by    = 0.01
        )

        # In case we are looking at the first step, we can only calculate the
        # new absolute turning angle by first deriving the difference in
        # relative turning angles
        if (i == 1) {
          relta_diff <- relta_new - relta
          absta <- absta + relta_diff
        } else {
          absta <- absta + relta_new
        }

        # Compute new endpoints
        absta <- ifelse(absta < 0, 2 * pi + absta, absta)
        absta <- ifelse(absta > 2 * pi, absta - 2 * pi, absta)
        x <- x + sin(absta) * sl
        y <- y + cos(absta) * sl
        relta <- relta_new
      }

      # Compute step length and relative turning angle from the origins to the
      # end coordinate of the last multi-random-steps
      dx <- x - rand$data[[z]]$x
      dy <- y - rand$data[[z]]$y
      sl <- sqrt(dx ** 2 + dy ** 2)
      absta <- atan2(dy, dx)
      absta <- (absta - pi / 2) * (-1)
      absta <- ifelse(absta < 0, 2 * pi + absta, absta)
      absta <- ifelse(absta > 2 * pi, absta - 2 * pi, absta)
      relta <- absta - rand$data[[z]]$absta
      relta <- ifelse(relta > +pi, relta - 2 * pi, relta)
      relta <- ifelse(relta < -pi, 2 * pi + relta, relta)
      rand$data[[z]]$sl <- sl
      rand$data[[z]]$relta_new <- relta
      return(rand$data[[z]])
    })
    rand <- unnest(rand, data)
    rand <- ungroup(rand)
    rand <- arrange(rand, ID, step_id, desc(case))
  } else {
    stop("Provide valid input for the desired approach")
  }

  # Calculate new "absolute" turning angle
  rand$relta_diff <- rand$relta_new - rand$relta
  rand$absta <- rand$absta + rand$relta_diff
  rand$absta <- ifelse(rand$absta < 0, 2 * pi + rand$absta, rand$absta)
  rand$absta <- ifelse(rand$absta > 2 * pi, rand$absta - 2 * pi, rand$absta)
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
computeCovars <- function(data, multicore = F) {

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

  # Put extracted covariates into a dataframe and bind them to the original data
  extracted <- as.data.frame(do.call(rbind, extracted))
  data <- cbind(data, extracted)

  # Ensure that step lengths cover a minimal distance
  data$sl[data$sl == 0] <- min(data$sl[data$sl != 0])

  # Specify phases of activity
  data$Activity <- ifelse(hour(data$TimestampRounded) == 7, "Inactive", "Active")

  # Calculate movement metrics
  data$log_sl <- log(data$sl)
  data$cos_ta <- cos(data$relta)

  # Compute step length in kilometers
  data$sl_km <- data$sl / 1000
  data$log_sl_km <- log(data$sl_km)

  # Return the results
  return(data)
}

# Function to run the step selection model
runModel <- function(data) {

  # Calculate distance to water in kilometers
  data$DistanceToWater <- data$DistanceToWater / 1000

  # Compute square rooted distance
  data$SqrtDistanceToWater <- sqrt(data$DistanceToWater)

  # Prepare model formula
  form <- (case ~
    + cos_ta
    + sl_km
    + log_sl_km
    + Water
    + SqrtDistanceToWater
    # + HumanInfluence
    # + Trees
    # + Shrubs
    # + Activity:sl_km
    # + Activity:log_sl_km
    + Water:sl_km
    # + Water:log_sl_km
    + cos_ta:sl_km
    # + cos_ta:log_sl_km
    + (1|step_id)
    + (0 + cos_ta|ID)
    + (0 + sl_km|ID)
    + (0 + log_sl_km|ID)
    + (0 + Water|ID)
    + (0 + SqrtDistanceToWater|ID)
    # + (0 + Trees|ID)
    # + (0 + HumanInfluence|ID)
    # + (0 + Shrubs|ID)
  )

  # The model might not converge well, thus, we need to write an error handling
  # function
  result <- tryCatch(
      expr    = {
        mod <- glmm_clogit(form, data)
        res <- list("Success", mod)
        return(res)
      }
    , error   = function(e) {return(list(e, mod))}
    , warning = function(w) {return(list(w, mod))}
  )

  # Run model
  return(result)

}

# Before we can test the functions above, we need some distributions from which
# we can sample step lengths and turning angles to generate random steps. Thus,
# let's fit step-length and turning angle distributions to the minimum step
# duration
metrics <- obs %>%
  rarifyData(missingness = 0) %>%
  computeBursts(forgiveness = 1) %>%
  computeMetrics() %>%
  subset(duration == 1 * 4) %>%
  select(sl, relta) %>%
  tibble()
sl_dist <- fit_distr(metrics$sl, dist_name = "gamma")$params
ta_dist <- fit_distr(metrics$relta, dist_name = "vonmises")$params
dists <- list(uncorrected = list(sl = sl_dist, ta = ta_dist))

# Try out the functions to see how they work
testing <- rarifyData(obs, missingness = 0.5)
testing <- computeBursts(testing, forgiveness = 1)
testing <- computeMetrics(testing)
testing <- computeSSF(testing, n_rsteps = 10, approach = "naive")
testing <- computeCovars(testing, multicore = T)
testing <- runModel(testing)
summary(testing)

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
