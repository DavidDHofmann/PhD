################################################################################
#### Integrated Step Selection Analysis of Simulated Data
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
library(amt)           # To fit distributions

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/00_Functions.r")

# Suppress scientific notation
options(scipen = 999)

# Load observed movement data
obs <- read_csv("03_Data/01_RawData/ObservedMovements.csv")

# Define the original distributions
sl_dist <- list(shape = 3, scale = 1)
ta_dist <- list(kappa = 0.5, mu = 0)
dists <- list(uncorrected = list(sl = sl_dist, ta = ta_dist))

# Load covariate layers
covars <- stack("03_Data/01_RawData/CovariateLayers.grd")
covars <- readAll(covars)

# Specify the different design combinations through which we want to run. Note
# that a forgiveness of 1 refers to a regular step selection function
dat <- expand_grid(
    Missingness    = seq(0, 0.8, by = 0.1)                  # Fraction of the fixes that is removed
  , Forgiveness    = 1:5                                    # Allowed lag of steps (in steps)
  , Replicate      = 1:100                                  # Number of replicates for each combination
  , Distributions  = c("uncorrected", "naive", "dynamic")   # Which distributions to use
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

    # Compute bursts by id. A new burst is defined if the duration is >
    # forgiveness
    mutate(data = map(data, function(x) {
      x$duration <- lead(x$step_number) - x$step_number
      x$irregular <- x$duration > forgiveness
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
    + forest
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
      dplyr::select(duration, sl, relta) %>%
      subset(duration <= max(dat$Forgiveness)) %>%
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

################################################################################
#### Fit Models
################################################################################
# Define the number of random steps
n_rsteps <- 10

# Go through the design and run step selection analysis with the specified
# parameters
dat$Coefs <- pbmclapply(
  # dat$Coefs <- lapply(
    X                  = 1:nrow(dat)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {

  # Extract important information from dataframe
  miss <- dat$Missingness[[x]]
  forg <- dat$Forgiveness[[x]]
  dist <- dat$Distributions[[x]]

  # Prepare data for ssf
  data <- obs %>%
    rarifyData(missingness = miss) %>%
    computeBursts(forgiveness = forg) %>%
    computeMetrics() %>%
    computeSSF(n_rsteps = n_rsteps, distributions = dist) %>%
    computeCovars(covars, multicore = F)

  # Run conditional logistic regression
  coefs <- runModel(data)

  # Return the coefficients
  return(coefs)
})

# Store results to file
write_rds(dat, "03_Data/03_Results/Models.rds")
dat <- read_rds("03_Data/03_Results/Models.rds")

# Visualize Results
unnest(dat, Coefs) %>%
  group_by(Missingness, Forgiveness, AdjustDists, Coefficient) %>%
  summarize(
    , SD           = sd(Estimate)
    , Estimate     = mean(Estimate)
    , LCI = Estimate - 2 * SD
    , UCI = Estimate + 2 * SD
    , .groups      = "drop"
  ) %>%
  ggplot(aes(x = Missingness, y = Estimate, col = as.factor(Forgiveness), fill = as.factor(Forgiveness))) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5, lwd = 0.5) +
    # geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.05, alpha = 0.5) +
    # geom_point() +
    # geom_line() +
    facet_wrap(~ AdjustDists + Coefficient, scales = "free", nrow = 2) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_viridis_d() +
    scale_color_viridis_d()
    # scale_fill_manual(values = c("orange", "cornflowerblue")) +
    # scale_color_manual(values = c("orange", "cornflowerblue"))

################################################################################
#### CAN WE USE SIMEX TO BACKTRANFORM?
################################################################################
################################################################################
#### What if we fit a separate step length distribution?
################################################################################
