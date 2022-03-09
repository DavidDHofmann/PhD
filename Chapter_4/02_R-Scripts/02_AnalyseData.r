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
glimpse(obs)

# Load covariate layers
covars <- stack("03_Data/01_RawData/CovariateLayers.grd")
covars <- readAll(covars)
glimpse(covars)

# Extend covariate layers slightly
covars <- extendRaster(covars, extent(c(-50, 350, -50, 350)))
plot(covars)

# Specify the different design combinations through which we want to run. Note
# that a forgiveness of 1 refers to a regular step selection function
dat <- expand_grid(
    Missingness    = seq(0, 0.5, by = 0.1)                                         # Fraction of the fixes that is removed
  , Forgiveness    = 1:5                                                           # Allowed lag of steps (in steps)
  , Replicate      = 1:100                                                         # Number of replicates for each combination
  , Approach       = c("uncorrected", "naive", "dynamic", "multistep", "model")    # Which approach to use to generate random steps
)

################################################################################
#### Functions
################################################################################
# Function to create a dataset with rarified observations (the same function is
# also used to subsample x individuals from all individuals)
rarifyData <- function(data, missingness, n_id = NULL) {
  if (is.null(n_id)) {
    keep <- unique(data$ID)
  } else {
    keep <- sample(unique(data$ID), n_id, replace = F)
  }
  data_sub <- subset(data, ID %in% keep)
  rarified <- data_sub[sort(sample(1:nrow(data_sub), size = nrow(data_sub) * (1 - missingness))), ]
  return(rarified)
}

# Function to compute bursts (per ID, depending on the fogriveness). A new burst
# always starts if the step-duration exceeds the forgiveness.
computeBursts <- function(data, forgiveness) {

  # Nest data by id (a burst cannot expand across multiple ids)
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

# Function to compute step metrics (step length, relative turning angle,
# absolute turning angle). Step metrics are calculated on bursts.
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
computeCovars <- function(data, covariates, multicore = F) {

  # Run the covariate extraction on multiple cores
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

  # Run the covariate extraction on a single core
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

  # Put extracted covariates into a dataframe and bind them to the original data
  extracted <- as.data.frame(do.call(rbind, extracted))
  data <- cbind(data, extracted)

  # Ensure that step lengths cover a minimal distance
  data$sl[data$sl == 0] <- min(data$sl[data$sl != 0])

  # Calculate derived movement metrics
  data$log_sl <- log(data$sl)
  data$cos_ta <- cos(data$relta)

  # Return the results
  return(data)
}

# Function to run the step selection model using two different approaches
runModel <- function(data, approach) {

  # Run the step selection model
  if (approach == "model" & length(unique(data$duration)) > 1) {
    data$duration <- as.factor(data$duration)
    mod <- clogit(case ~
      + sl:duration
      + log_sl:duration
      + cos_ta:duration
      + forest
      + elev
      + dist
      + strata(step_id)
      , data = data
    )
  } else {
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
  }

  # Extract model coefficients and put them into a dataframe. Also, compute
  # confidence intervals.
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

# Before we can test the functions above, we need some distributions from which
# we can sample step lengths and turning angles to generate random steps. Thus,
# let's fit step-length and turning angle distributions to the minimum step
# duration
metrics <- obs %>%
  rarifyData(missingness = 0, n_id = 100) %>%
  computeBursts(forgiveness = 1) %>%
  computeMetrics() %>%
  subset(duration == 1) %>%
  select(sl, relta) %>%
  tibble()
sl_dist <- fit_distr(metrics$sl, dist_name = "gamma")$params
ta_dist <- fit_distr(metrics$relta, dist_name = "vonmises")$params
dists <- list(uncorrected = list(sl = sl_dist, ta = ta_dist))

# Try out the functions to see how they work
testing <- rarifyData(obs, missingness = 0.5, n_id = 50)
testing <- computeBursts(testing, forgiveness = 2)
testing <- computeMetrics(testing)
testing <- computeSSF(testing, n_rsteps = 10, approach = "uncorrected")
testing <- computeCovars(testing, covars, multicore = T)
testing <- runModel(testing, approach = "uncorrected")
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
      rarifyData(missingness = 0.5, n_id = 100) %>%
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

# So that we can independently store the results from different iterations,
# let's prepare a filename for each row of the design table
dat$Filename <- paste0(
    "03_Data/03_Results/ModelResults/Simulation/"
  , "M", sprintf("%02d", as.integer(dat$Missingness * 100)), "_"
  , "F", sprintf("%02d", dat$Forgiveness), "_"
  , "R", sprintf("%03d", dat$Replicate), "_"
  , dat$Approach
  , ".rds"
)

# Create the respective folder
dir.create("03_Data/03_Results/ModelResults", showWarnings = F)
dir.create("03_Data/03_Results/ModelResults/Simulation", showWarnings = F)

# Let's randomize the design matrix
dat <- dat[sample(nrow(dat)), ]

# Subset to rows that haven't been run yet
dat <- subset(dat, !file.exists(Filename))

# Go through the design and run step selection analysis with the specified
# parameters
pbmclapply(
  # dat$Coefs <- lapply(
    X                  = 1:nrow(dat)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {

  # Extract important information from dataframe
  miss <- dat$Missingness[[x]]
  forg <- dat$Forgiveness[[x]]
  appr <- dat$Approach[[x]]
  file <- dat$Filename[[x]]

  # Prepare data for ssf
  data <- obs %>%
    rarifyData(missingness = miss, n_id = 100) %>%
    computeBursts(forgiveness = forg) %>%
    computeMetrics() %>%
    computeSSF(n_rsteps = n_rsteps, approach = appr) %>%
    computeCovars(covars, multicore = F)

  # Run conditional logistic regression
  coefs <- runModel(data, approach = appr)

  # Store results to file
  write_rds(coefs, file)
  return(NULL)
})

################################################################################
#### CONTINUE HERE
################################################################################

# Store results to file
write_rds(dat, "03_Data/03_Results/Models.rds")
dat <- read_rds("03_Data/03_Results/Models.rds")

# Let's also prepare a dataframe containing the "truth"
truth <- data.frame(
    Coefficient = c("sl", "log_sl", "cos_ta", "forest", "elev", "dist")
  , Estimate    = c(0, 0, 0, -1, 0.5, -20)
)

# Visualize Results
unnest(dat, Coefs) %>%
  mutate(Approach = Distributions) %>%
  group_by(Missingness, Forgiveness, Approach, Coefficient) %>%
  summarize(
    , SD       = sd(Estimate)
    , Estimate = mean(Estimate)
    , LCI      = Estimate - 2 * SD
    , UCI      = Estimate + 2 * SD
    , .groups  = "drop"
  ) %>%
  ggplot(aes(x = Missingness, y = Estimate, col = as.factor(Forgiveness), fill = as.factor(Forgiveness))) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5, lwd = 0.5) +
    geom_hline(data = truth, aes(yintercept = Estimate), lty = 2, lwd = 0.5) +
    # geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.05, alpha = 0.5) +
    # geom_point() +
    # geom_line() +
    facet_wrap(~ Approach + Coefficient, scales = "free", nrow = 3) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_viridis_d() +
    scale_color_viridis_d()
    # scale_fill_manual(values = c("orange", "cornflowerblue")) +
    # scale_color_manual(values = c("orange", "cornflowerblue"))
