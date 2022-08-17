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
library(crawl)         # To impute missing fixes

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")
# setwd("C:/Users/david/switchdrive/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Suppress scientific notation
options(scipen = 999)

# Load simulated data
obs <- read_rds("03_Data/01_RawData/SimulatedMovement.rds")
cov <- read_rds("03_Data/01_RawData/SimulatedCovariates.rds")

# Generate unique step_ids
obs <- obs %>%
  unnest(Simulations) %>%
  mutate(step_id = 1:n()) %>%
  nest(Simulations = -c(AutocorrRange, ID))

# Add an artifical timestamp to the simulated GPS data
obs <- obs %>%
  mutate(Simulations = pbmclapply(Simulations
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x) {
      x$timestamp <- lubridate::ymd_hms("2000-01-01 00:00:00") + hours(1:nrow(x))
    return(x)
  }))

# Extend covariate layers slightly
cov <- cov %>%
  mutate(Covariates = pbmclapply(Covariates
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x) {
      covars <- extendRaster(x, extent(c(-50, 350, -50, 350)))
      covars <- writeRaster(covars, tempfile())
    return(covars)
  }))

# Specify the different design combinations through which we want to run
dat <- expand_grid(
    Missingness    = seq(0, 0.5, by = 0.1)                                                 # Fraction of the fixes that is removed
  , Forgiveness    = 1:5                                                                   # Allowed step-duration
  , AutocorrRange  = unique(obs$AutocorrRange)                                             # Autocorrelation range in the covariates
  , Replicate      = 1:100                                                                 # Number of replicates for each combination
  , Approach       = c("uncorrected", "naive", "dynamic", "model", "multistep", "imputed") # Which approach to use to analyse the data
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
computeCovars <- function(data, covariates, multicore = F) {

  # Pair simulated movement data with the appropriate covariate layer
  data <- data %>%
    nest(Data = -c(AutocorrRange, ID)) %>%
    left_join(., cov, by = c("ID", "AutocorrRange"))

  # Run the covariate extraction on multiple cores
  if (multicore) {
    data$Data <- pbmclapply(1:nrow(data), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
      gps_data <- data$Data[[x]]
      cov_data <- data$Covariates[[x]]
      extracted <- lapply(1:nrow(gps_data), function(y) {
        ints <- interpolatePoints(
            x1 = gps_data$x[y]
          , x2 = gps_data$x_to[y]
          , y1 = gps_data$y[y]
          , y2 = gps_data$y_to[y]
          , by = 1
        )
        extr <- raster::extract(cov_data, ints)
        extr <- colMeans(extr)
        return(extr)
      })
      extracted <- as.data.frame(do.call(rbind, extracted))
      gps_data <- cbind(gps_data, extracted)
      return(gps_data)
    })

  # Run the covariate extraction on a single core
  } else {
    data$Data <- lapply(1:nrow(data), function(x) {
      gps_data <- data$Data[[x]]
      cov_data <- data$Covariates[[x]]
      extracted <- lapply(1:nrow(gps_data), function(y) {
        ints <- interpolatePoints(
            x1 = gps_data$x[y]
          , x2 = gps_data$x_to[y]
          , y1 = gps_data$y[y]
          , y2 = gps_data$y_to[y]
          , by = 1
        )
        extr <- raster::extract(cov_data, ints)
        extr <- colMeans(extr)
        return(extr)
      })
      extracted <- as.data.frame(do.call(rbind, extracted))
      gps_data <- cbind(gps_data, extracted)
      return(gps_data)
    })
  }

  # Unnest again
  data <- data %>% select(-Covariates) %>% unnest(Data)

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

# Function to fit step length and turning angle distributions. Statically, as
# well as dynamically to different step durations
fitDists <- function(data, replicate = 100) {

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
  dists_dynamic <- lapply(1:replicate, function(x) {
    params <- data %>%
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

  # Return the distributions
  return(dists)
}

# Let's test if the above functions work as they should
test <- unnest(obs, Simulations)
test <- subset(test, AutocorrRange == 1)
test <- rarifyData(test, missingness = 0.5, n_id = 50)
testi <- imputeFixes(test)
dists <- fitDists(test, replicate = 10)
test <- computeBursts(test, forgiveness = 1)
test <- computeMetrics(test)
testing1 <- computeSSF(test, n_rsteps = 10, approach = "uncorrected", dists = dists)
testing2 <- computeSSF(test, n_rsteps = 10, approach = "multistep", dists = dists)
par(mfrow = c(2, 1))
hist(testing1$relta[testing1$case == 0])
hist(testing2$relta[testing2$case == 0])
fit_distr(testing1$relta[testing1$case == 0], dist_name = "vonmises")$params
fit_distr(testing2$relta[testing2$case == 0], dist_name = "vonmises")$params
testing <- computeCovars(testing1, covars, multicore = T)
testing <- runModel(testing, approach = "uncorrected")
testing

# Try to impute missing fixes
testing <- imputeFixes(rarifyData(obs, missingness = 0.5, n_id = 10))
head(testing)

# Function to fit step length and turn angle distributions uncorrected and
# dynamically
data <- obs %>% subset(AutocorrRange == 1) %>% unnest(Simulations)

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
  , "A", sprintf("%02d", dat$AutocorrRange), "_"
  , "R", sprintf("%03d", dat$Replicate), "_"
  , dat$Approach
  , ".rds"
)

# Create the respective folder
dir.create("03_Data/03_Results/ModelResults", showWarnings = F)
dir.create("03_Data/03_Results/ModelResults/Simulation", showWarnings = F)

# Write the design to file
write_rds(dat, "03_Data/03_Results/SimulationDesign.rds")

# Let's randomize the design matrix
dat <- dat[sample(nrow(dat)), ]

# Subset to rows that haven't been run yet
dat_sub <- subset(dat, !file.exists(Filename))

# Go through the design and run step selection analysis with the specified
# parameters
pbmclapply(
  # dat$Coefs <- lapply(
    X                  = 1:nrow(dat_sub)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {

  # Extract important information from dataframe
  miss <- dat_sub$Missingness[[x]]
  forg <- dat_sub$Forgiveness[[x]]
  auto <- dat_sub$AutocorrRange[[x]]
  appr <- dat_sub$Approach[[x]]
  file <- dat_sub$Filename[[x]]

  # Prepare data for ssf
  data <- obs %>%
    subset(AutocorrRange == auto) %>%
    unnest(Simulations) %>%
    rarifyData(missingness = miss, n_id = 100)

  # Compute turning angle and step length distributions
  dists <- fitDists(data, replicate = 100)

  # Impute if necessary
  if (appr == "imputed") {
    data <- imputeFixes(data)
  }

  # Continue with analysis
  data <- data %>% computeBursts(forgiveness = forg) %>%
    computeMetrics() %>%
    computeSSF(n_rsteps = n_rsteps, approach = appr, dists = dists) %>%
    computeCovars(covars, multicore = F)

  # Run conditional logistic regression
  coefs <- runModel(data, approach = appr)

  # Store results to file
  write_rds(coefs, file)
  return(NULL)
})
