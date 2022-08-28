################################################################################
#### Simulating Movement Data with Known Preferences
################################################################################
# Description: Simulating movement on virtual landscape and estimate habitat and
# movement preferences using different flavors of step selection.

# Clear R's brain
rm(list = ls())

# Load required packages
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(NLMR)          # To simulate covariates
library(raster)        # To handle spatial data
library(sf)            # For plotting spatial features
library(amt)           # To fit distributions
library(survival)      # To run conditional logistic regression
library(crawl)         # To impute missing fixes

# Need to ensure the following versions are installed!
if (packageVersion("RandomFieldsUtils") != "0.5") {
  stop("Please install RandomFieldsUtils version 0.5")
}
if (packageVersion("RandomFields") != "3.3.1") {
  stop("Please install RandomFields version 3.3.1")
}
# devtools::install_version("RandomFieldsUtils", version = "0.5", repos = "http://cran.us.r-project.org")
# devtools::install_version("RandomFields", version = "3.3.1", repos = "http://cran.us.r-project.org", dependencies = F)

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")
# setwd("C:/Users/david/switchdrive/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Simulation Parameters
################################################################################
n            <- 300                                                          # Resolution of the covariate layers
n_rsteps     <- 10                                                           # Number of random steps to be generated
n_steps      <- 100                                                          # Number of consecutive steps to be simulated
n_inds       <- 100                                                          # How many individuals to simulate for each treatment
formula      <- ~ forest + elev + dist                                       # Formula used to predict step-selection scores
prefs        <- c(-1, 0.5, -20)                                              # Preferences of the simulated individuals
sl_dist      <- list(name = "gamma", params = list(shape = 3, scale = 1))    # Step-length parameters
ta_dist      <- list(name = "vonmises", params = list(kappa = 0.5, mu = 0))  # Turning-angle parameters
stop         <- F                                                            # Should the simulation terminate at boundaries?
maxattempts  <- 5                                                            # Maximum number of times to repeat an analysis when it fails to converge

# Let's specify the extent on which animals are allowed to move
ext <- extent(0, n, 0, n)
ext <- as(ext, "SpatialPolygons")

# Let's also specify an extent within which individuals will be released
ext2 <- extent(ext) - 100
ext2 <- as(ext2, "SpatialPolygons")

# And an extent to which we will expand covariates before extracting covariate
# values
ext3 <- extent(c(-50, 350, -50, 350))
ext3 <- as(ext3, "SpatialPolygons")

# Generate a point of attraction at the center of the study area and compute the
# distance to it
r <- raster(nrows = n, ncols = n, xmn = 0, xmx = n, ymn = 0, ymx = n)
cent <- SpatialPoints(t(c(n / 2, n / 2)))
dist <- distanceFromPoints(r, cent)
dist <- (dist - cellStats(dist, min)) / (cellStats(dist, max) - cellStats(dist, min))

################################################################################
#### Main Simulation and Analysis Functions
################################################################################
# Function to simulate the different covariates and return them as a stack
simCovars <- function(autocorr_range = 10, proportion_forest = 0.5) {

  # Simulate forest layer
  forest <- nlm_gaussianfield(
      ncol           = n
    , nrow           = n
    , autocorr_range = autocorr_range
    , mag_var        = 1
    , nug            = 0
  )
  cutoff <- quantile(forest, 1 - proportion_forest)
  forest <- forest > cutoff
  forest <- focal(forest
    , w        = matrix(rep(1, 9), nrow = 3)
    , fun      = function(x) {sum(x) > 4.5}
    , pad      = T
    , padValue = 0
  )

  # Simulate elevation layer
  elev <- nlm_gaussianfield(
      ncol           = n
    , nrow           = n
    , autocorr_range = autocorr_range
    , mag_var        = 1
    , nug            = 0
  )

  # Put them into a stack
  covars <- stack(forest, elev, dist)
  names(covars) <- c("forest", "elev", "dist")
  return(covars)
}

# Function to simulate movement on a simulated landscape. Returns "observed" GPS
# data
simMove <- function(covars, n_id, multicore = F) {

  # Simulate movement for every individual
  if (multicore) {
    sims <- pbmclapply(1:n_id, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {

      # Simulate movement
      sim <- move(
          xy       = matrix(runif(2, xmin(ext2), xmax(ext2)), ncol = 2)
        , covars   = covars
        , formula  = formula
        , prefs    = prefs
        , sl_dist  = sl_dist
        , ta_dist  = ta_dist
        , ext      = ext
        , n_steps  = n_steps
        , n_rsteps = n_rsteps
        , stop     = stop
      )

      # Add an artifial timestamp and the ID to the data
      sim$timestamp <- lubridate::ymd_hms("2000-01-01 00:00:00") + hours(1:nrow(sim))
      sim$ID <- x
      return(sim)
    })
  } else {
    sims <- lapply(1:n_id, function(x) {

      # Simulate movement
      sim <- move(
          xy       = matrix(runif(2, xmin(ext2), xmax(ext2)), ncol = 2)
        , covars   = covars
        , formula  = formula
        , prefs    = prefs
        , sl_dist  = sl_dist
        , ta_dist  = ta_dist
        , ext      = ext
        , n_steps  = n_steps
        , n_rsteps = n_rsteps
        , stop     = stop
      )

      # Add an artifial timestamp, the ID, as well as a unique step_id to the
      # data
      sim$timestamp <- lubridate::ymd_hms("2000-01-01 00:00:00") + hours(1:nrow(sim))
      sim$ID <- x
      return(sim)
    })
  }

  # Bind simulations into a single dataframe and assign a unique step id
  sims <- do.call(rbind, sims)
  sims$step_id <- 1:nrow(sims)

  # In reality, we don't observe step lengths etc.
  sims$absta <- NULL
  sims$ta    <- NULL
  sims$sl    <- NULL
  return(sims)
}

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

# Function to impute missing fixes to get a regularized dataframe. For
# simplicity, only conduct a single imputation.
imputeFixes <- function(data) {
  inds <- unique(data$ID)
  pred <- lapply(inds, function(x) {

    # Fit model but suppress the unneccessary messages
    sink(tempfile())
    mod <- suppressMessages(suppressWarnings(crwMLE(
        mov.model   = ~ 1
      , data        = data[data$ID == x, ]
      , Time.name   = "timestamp"
      , coord       = c("x", "y")
      , timeStep    = "hour"
      # , method      = "L-BFGS-B"
      , control     = list(maxit = 50, trace = 0, REPORT = 1)
      , initialSANN = list(maxit = 100, trace = 1, REPORT = 1, temp = 0.5, tmax = 10)
      , attemps     = 200
    )))
    sink()

    # Make prediction
    pre <- suppressMessages(suppressWarnings(crwPredict(
        object.crwFit = mod
      , predTime      = "1 hour"
      , flat          = TRUE
    )))

    # Clean the output
    pre <- data.frame(
        ID          = pre$ID
      , x           = pre$mu.x
      , y           = pre$mu.y
      , step_number = 1:nrow(pre)
      , step_id     = NA
      , timestamp   = pre$timestamp
    )

    # Indicate if a fix was imputed or not
    pre$imputed <- !(pre$timestamp %in% data$timestamp[data$ID == x])

    # Return the imputed dataframe
    return(pre)
  })

  # Bind all individuals together, assign unique step id, return all
  pred <- do.call(rbind, pred)
  pred$step_id <- 1:nrow(pred)
  return(pred)
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

# Function to fit step length and turning angle distributions. Statically, as
# well as dynamically to different step durations
fitDists <- function(data, replicate = 100, max_forgiveness = 5, dynamic = T, multicore = F) {

  # Fit distributions to a step duration of one
  metrics <- data %>%
    rarifyData(missingness = 0) %>%
    computeBursts(forgiveness = 1) %>%
    computeMetrics() %>%
    subset(duration == 1) %>%
    select(sl, relta) %>%
    tibble()
  sl_dist <- fit_distr(metrics$sl, dist_name = "gamma")$params
  ta_dist <- fit_distr(metrics$relta, dist_name = "vonmises")$params
  dists <- list(uncorrected = list(sl = sl_dist, ta = ta_dist))

  # Fit distributions to various step durations
  if (dynamic) {
    if (multicore) {
      dists_dynamic <- pbmclapply(1:replicate, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
        params <- data %>%
          rarifyData(missingness = 0.5) %>%
          computeBursts(forgiveness = max_forgiveness) %>%
          computeMetrics() %>%
          dplyr::select(duration, sl, relta) %>%
          subset(duration <= max_forgiveness) %>%
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
    } else {
      dists_dynamic <- lapply(1:replicate, function(x) {
        params <- data %>%
          rarifyData(missingness = 0.5) %>%
          computeBursts(forgiveness = max_forgiveness) %>%
          computeMetrics() %>%
          dplyr::select(duration, sl, relta) %>%
          subset(duration <= max_forgiveness) %>%
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
    }

    # Summarize values
    dists_means <- dists_dynamic %>%
      pivot_longer(shape:mu, names_to = "Parameter", values_to = "Value") %>%
      group_by(duration, Parameter) %>%
      summarize(mean = mean(Value), sd = sd(Value), .groups = "drop") %>%
      dplyr::select(duration, Parameter, mean) %>%
      pivot_wider(names_from = Parameter, values_from = mean)
    dists_sds <- dists_dynamic %>%
      pivot_longer(shape:mu, names_to = "Parameter", values_to = "Value") %>%
      group_by(duration, Parameter) %>%
      summarize(mean = mean(Value), sd = sd(Value), .groups = "drop") %>%
      dplyr::select(duration, Parameter, sd) %>%
      pivot_wider(names_from = Parameter, values_from = sd) %>%
      setNames(paste0(names(.), "_sd"))
    dists_means <- cbind(dists_means, dists_sds[, -1])

    # Put values together with the uncorrected distributions
    dists$dynamic <- list(
          sl = select(dists_means, duration, shape, scale, shape_sd, scale_sd)
        , ta = select(dists_means, duration, kappa, mu, kappa_sd, mu_sd)
      )
  }

  # Return the distributions
  return(dists)
}

# Function to generate random steps. The approach parameter is used to specify
# the approach that should be used to generate random steps
computeSSF <- function(data, n_rsteps, approach, dists) {

  # # TESTING
  # data <- testing[2, ]
  # data$x <- 0
  # data$y <- 0
  # # # data$duration <- c(1, 3)
  # data$duration <- 3
  # data$sl <- 1
  # data$absta <- 0
  # data$relta <- 0
  # n_rsteps <- 1
  # approach <- "multistep"

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
  if (approach == "uncorrected" | approach == "imputed") {
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
  #### Approach 3 - Dynamic or Modelled
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
  # the number of generated steps matches the step duration of the observed step
  } else if (approach == "multistep") {
    rand <- rand %>% group_by(duration) %>% nest()
    rand$data <- lapply(1:nrow(rand), function(z) {

      # Extract important info of the current iteration (i.e. the "z")
      duration <- rand$duration[z]
      x        <- rand$data[[z]]$x
      y        <- rand$data[[z]]$y
      relta    <- rand$data[[z]]$relta
      absta    <- rand$data[[z]]$absta

      # Generate sequence of random steps (if duration = 1, generate only one
      # step, if duration = 2, generate two, etc.)
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

        # # TESTING
        # sl <- c(1, 1, 1)[i]
        # relta_new <- c(pi/2, 0, -pi/2)[i]
        # sl <- 1
        # relta_new <- 0

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
      # relta <- absta - rand$data[[z]]$absta           # THIS IS CAUSING ISSUES!!!!
      relta <- rand$data[[z]]$relta - (rand$data[[z]]$absta - absta)
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
computeCovars <- function(data, covars, multicore = F) {

  # Extend covariates slightly to ensure that random steps are not outside our
  # borders
  covars <- extendRaster(covars, ext3)

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
      extr <- raster::extract(covars, ints)
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
      extr <- raster::extract(covars, ints)
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

# ################################################################################
# #### Ensure All Functions Work as Intended
# ################################################################################
# # Simulate covariates and across the layers
# cov <- simCovars(autocorr_range = 50, proportion_forest = 0.5)
# obs <- simMove(cov, n_id = 10, multicore = T)
# as.data.frame(cov, xy = T) %>%
#   gather(key = covariate, value = value, 3:5) %>%
#   ggplot(aes(x = x, y = y, fill = value)) +
#     geom_raster() +
#     geom_sf(data = st_as_sf(ext2), inherit.aes = F, fill = NA, col = "white", lty = 2) +
#     geom_path(data = obs, aes(x = x, y = y, col = as.factor(ID)), inherit.aes = F) +
#     scale_fill_viridis_c(option = "viridis") +
#     coord_sf() +
#     theme_minimal() +
#     facet_wrap("covariate") +
#     theme(
#         axis.title.y    = element_text(angle = 0, vjust = 0.5)
#       , legend.position = "none"
#     )
#
# # Rarify data
# obs_missing <- rarifyData(obs, missingness = 0.2, n_id = 9)
#
# # Try to impute missing fixes
# obs_imputed <- imputeFixes(obs_missing)
# ggplot(obs_imputed, aes(x = x, y = y, group = ID, col = imputed)) +
#   geom_path(lwd = 0.2) +
#   geom_point(size = 0.2) +
#   coord_sf() +
#   theme_minimal() +
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
#
# # Compute bursts from the (non-imputed) data -> Try with different forgiveness
# # values!
# obs_bursted <- computeBursts(obs_missing, forgiveness = 2)
# ggplot(obs_bursted, aes(x = x, y = y, group = ID, col = as.factor(burst))) +
#   geom_path(lwd = 0.2) +
#   geom_point(size = 0.2) +
#   coord_sf() +
#   theme_minimal() +
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
#
# # Fit step- and turning-angle distributions to original and rarified data
# fit_distris <- fitDists(obs_bursted, max_forgiveness = 2, dynamic = F)
# fit_distris <- fitDists(obs_bursted, replicate = 10, max_forgiveness = 2, dynamic = T)
#
# # Compute step metrics
# obs_metrics <- computeMetrics(obs_bursted)
#
# # Generate random steps
# obs_stepsel <- computeSSF(obs_metrics, n_rsteps = 10, approach = "uncorrected", dists = fit_distris)
#
# # Extract covariates
# obs_covaris <- computeCovars(obs_stepsel, cov, multicore = T)
#
# # Run the model
# mod <- runModel(obs_covaris, approach = "uncorrected")
# mod
#
# # Clear up space and remov unused objects
# rm(
#     cov
#   , obs
#   , mod
#   , obs_missing
#   , obs_imputed
#   , obs_bursted
#   , obs_metrics
#   , obs_stepsel
#   , obs_covaris
#   , fit_distris
# )
# gc()

################################################################################
#### Proper Analysis
################################################################################
# Specify the different design combinations through which we want to run
dat <- expand_grid(
    Missingness    = seq(0, 0.5, by = 0.1)                                                 # Fraction of the fixes that is removed
  , Forgiveness    = 1:5                                                                   # Allowed step-duration
  , AutocorrRange  = c(10, 50, 100)                                                        # Autocorrelation range in the covariates
  , Replicate      = 1:10                                                                 # Number of replicates for each combination
  , Approach       = c("uncorrected", "naive", "dynamic", "model", "multistep", "imputed") # Which approach to use to analyse the data
)

# Nest by replicate and autocorrelation range
dat <- nest(dat, Data = -c(Replicate, AutocorrRange))

# Let's prepare a filename for each row of the (nested) design table
dat$Filename <- paste0(
    "03_Data/ModelResults/"
  , "A", sprintf("%02d", dat$AutocorrRange), "_"
  , "R", sprintf("%03d", dat$Replicate)
  , ".rds"
)

# Create the respective folder
dir.create("03_Data/ModelResults", showWarnings = F)

# Write the design to file
write_rds(dat, "03_Data/SimulationDesign.rds")

# Let's randomize the design matrix
dat <- dat[sample(nrow(dat)), ]

# Subset to rows that haven't been run yet
dat_sub <- subset(dat, !file.exists(Filename))

# Loop through the design and run the simulations
dat_sub$Results <- lapply(1:nrow(dat_sub), function(x) {

  # Extract the authocorrelation range and the filename
  auto <- dat_sub$AutocorrRange[x]
  file <- dat_sub$Filename[x]
  desi <- dat_sub$Data[[x]]

  # Simulate covariates and movement
  cat("Simulating virtual landscape and dpsersal movements...\n")
  cov <- simCovars(autocorr_range = auto, proportion_forest = 0.5)
  obs <- simMove(covars = cov, n_id = n_inds, multicore = T)

  # Fitting step length and turning angle distributions (statically and
  # dynamically)
  cat("Fitting step length and turning angle distributions...\n")
  dis <- fitDists(obs
    , replicate       = 10
    , max_forgiveness = max(desi$Forgiveness)
    , dynamic         = T
    , multicore       = T
  )

  # Analyse the data using the various approaches
  desi$Coefs <- pbmclapply(
      X                  = 1:nrow(desi)
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(y) {

    # Extract relevant design parameters
    miss <- desi$Missingness[[y]]
    forg <- desi$Forgiveness[[y]]
    appr <- desi$Approach[[y]]

    # The model may not converge. Repeat the analysis if this is the case
    attempts <- 1
    coefs <- NA
    while (!is.data.frame(coefs) & attempts <= maxattempts) {

      # Rareify the data
      obs_sub <- rarifyData(obs, missingness = miss, n_id = NULL)

      # If required, impute fixes
      if (appr == "imputed") {
        obs_sub <- imputeFixes(obs_sub)
      }

      # Compute bursts and on those the relevant step metrics
      obs_sub <- computeBursts(obs_sub, forgiveness = forg)
      obs_sub <- computeMetrics(obs_sub)

      # Generate random steps
      obs_ssf <- computeSSF(obs_sub
        , n_rsteps = n_rsteps
        , approach = appr
        , dists    = dis
      )

      # Extract covariates
      obs_ssf <- computeCovars(obs_ssf
        , covars    = cov
        , multicore = F
      )

      # Run the model
      coefs <- tryCatch(runModel(obs_ssf, approach = appr)
        , warning = function(w) {return (NA)}
        , error   = function(e) {return (NA)}
      )

      # Update number of attempts
      attempts <- attempts + 1
    }

    # Return the results
    return(coefs)
  })

  # Unnest the design
  desi <- unnest(desi, Coefs)

  # Write results to file
  write_rds(desi, file)

  # Store results to file
  return(NULL)
})
