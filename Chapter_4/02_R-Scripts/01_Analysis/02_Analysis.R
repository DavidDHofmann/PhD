################################################################################
#### Analysis of Simulated Data using Different Approaches
################################################################################
# Description: Analysis of the generated data using different approaches

# Clear R's brain
rm(list = ls())

# Load required packages
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(raster)        # To handle spatial data
library(sf)            # For plotting spatial features
library(amt)           # To fit distributions
library(survival)      # To run conditional logistic regression
library(crawl)         # To impute missing fixes

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_4"
setwd(wd)

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Reload simulated data
dat <- read_rds("03_Data/Simulation.rds")

# Load all simulated data into memory
dat$Covariates <- lapply(dat$Covariates, readAll)
dat$Movement   <- lapply(dat$MovementFilename, read_rds)

# Remove any columns that are not needed anymore
dat <- select(dat, -c(CovariateFilename, MovementFilename))

################################################################################
#### Global Model Parameters
################################################################################
n_rsteps     <- 10                     # Number of random steps to be generated
formula      <- ~ forest + elev + dist # Formula used to predict step-selection scores
maxattempts  <- 5                      # Maximum number of times to repeat an analysis when it fails to converge

################################################################################
#### Main Analysis Functions
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

  # # Extend covariates slightly to ensure that random steps are not outside our
  # # borders
  # covars <- extendRaster(covars, ext3)

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
# # Extract some example data
# cov <- dat$Covariates[[1]]
# obs <- dat$Movement[[1]]
#
# # Rarify data
# obs_missing <- rarifyData(obs, missingness = 0.3)
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
# # Compute step metrics
# obs_metrics <- computeMetrics(obs_bursted)
#
# # Fit step- and turning-angle distributions to original and rarified data
# fit_distris <- fitDists(obs_bursted, max_forgiveness = 2, dynamic = F)
# fit_distris <- fitDists(obs_bursted, replicate = 10, max_forgiveness = 2, dynamic = T)
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
# Let's indicate the filenames into which we will store the analysis results
dat$ResultsFilename <- paste0(
    "03_Data/ModelResults/"
  , "A", sprintf("%02d", dat$AutocorrRange), "_"
  , "R", sprintf("%03d", dat$Replicate)
  , ".rds"
)

# Create the respective folder
dir.create("03_Data/ModelResults", showWarnings = F)

# For each set of covariates and simulations we will need to loop through
# different treatments. Let's create a dataframe to keep track of those
design <- expand_grid(
    Missingness    = seq(0, 0.5, by = 0.1)                                                 # Fraction of the fixes that is removed
  , Forgiveness    = 1:5                                                                   # Allowed step-duration
  , Approach       = c("uncorrected", "naive", "dynamic", "model", "multistep", "imputed") # Which approach to use to analyse the data
)

# Write the design to file (don't save the covariates and simualated movements
# again though)
write_rds(select(dat, -c(Covariates, Movement)), "03_Data/Analysis.rds")

# Let's randomize the order in which we will go through the different
# simulations
dat <- dat[sample(nrow(dat)), ]

# Subset to rows that haven't been run yet
dat_sub <- subset(dat, !file.exists(ResultsFilename))

# Loop through the simulations and apply the different analyses
lapply(1:nrow(dat_sub), function(x) {

  # Print progress
  cat("Iteration", x, "out of", nrow(dat_sub), "\n")

  # Extract relevant information
  cov <- dat_sub$Covariates[[x]]
  obs <- dat_sub$Movement[[x]]
  fil <- dat_sub$ResultsFilename[x]

  # Analyse the data using the various approaches
  desi <- design
  desi$Results <- pbmclapply(
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

      # Rarify the data
      obs_sub <- rarifyData(obs, missingness = miss)

      # Fit step length distributions
      if (appr %in% c("dynamic", "model")) {
          dis <- fitDists(obs_sub, dynamic = T, replicate = 10, max_forgiveness = forg)
        } else {
          dis <- fitDists(obs_sub, dynamic = F, replicate = 10, max_forgiveness = forg)
      }

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

    # Also keep track of the distribution used to draw random steps
    results <- list(Dists = dis, Coefs = coefs)

    # Return the results
    return(results)
  })

  # Write results to file
  write_rds(desi, fil)

  # Store results to file
  return(NULL)
})
