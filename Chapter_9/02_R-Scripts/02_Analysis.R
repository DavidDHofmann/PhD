################################################################################
#### Simulation and Analysis of Hidden Markov Data with Step Selection Functions
################################################################################
# Clear R's brain
rm(list = ls())

# Suppress scientific notation
options(scipen = 999)

# Load required packages
library(tidyverse)     # For data wrangling
library(raster)        # To handle spatial data
library(amt)           # To fit distributions
library(pbmcapply)     # To run stuff in parallel

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_9"
setwd(wd)

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Analysis Parameters
################################################################################
# Simulation setup
n_rsteps  <- 25
formula   <- ~ elev + dist

# Let's also reload the true simulation parameters
params_true <- read_rds("03_Data/SimulationParameters.rds")

################################################################################
#### Functions
################################################################################
# Function to generate random steps
computeSSF <- function(data, n_rsteps, dists) {

  # Generate a new column that indicates that the steps are "observed" steps
  data$case <- 1

  # Cannot work with steps that have no turning angle, so remove them
  data <- subset(data, !is.na(ta))

  # Create a new dataframe into which we can put alternative/random steps
  rand <- data[rep(1:nrow(data), each = n_rsteps), ]

  # Indicate that these steps are random steps (case = 0)
  rand$case <- 0

  # Step lengths sampled from "minimal" step-duration distributions
  rand$sl <- rgamma(n = nrow(rand)
    , scale = dists$sl$scale
    , shape = dists$sl$shape
  )
  rand$ta_new <- rvonmises(n = nrow(rand)
    , kappa = dists$ta$kappa
    , mu    = dists$ta$mu
    , by    = 0.01
  )

  # Calculate new "absolute" turning angle
  rand$ta_diff <- rand$ta_new - rand$ta
  rand$absta <- rand$absta + rand$ta_diff
  rand$absta <- ifelse(rand$absta < 0, 2 * pi + rand$absta, rand$absta)
  rand$absta <- ifelse(rand$absta > 2 * pi, rand$absta - 2 * pi, rand$absta)
  rand$ta <- rand$ta_new

  # Remove undesired stuff
  rand$ta_new <- NULL
  rand$ta_diff <- NULL

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
  rownames(all) <- NULL
  return(all)
}

# Function to extract covariates
computeCovars <- function(data, covars, multicore = F) {

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
  data$cos_ta <- cos(data$ta)

  # Return the results
  return(data)
}

# Function to sample parameter vector
sampleParams <- function(perfect = F) {
  if (perfect) {
    params_star <- c(
          qlogis(params_true[1:2])
        , log(params_true[3:8])
        , params_true[9:12]
    )
  } else {
    params_star <- c(
          qlogis(runif(n = 2, min = 0, max = 1))  # transition
        , log(runif(n = 4, min = 1, max = 5))     # step length
        , log(runif(n = 2, min = 0, max = 1))     # tunring angle
        , runif(n = 2, min = -2, max = 2)         # elevation
        , runif(n = 2, min = -20, max = 20)       # distance
    )
  }
  names(params_star) <- names(params_true)
  return(params_star)
}

# Function to backtransform parameters
backtransformParams <- function(theta_est) {
  c(plogis(theta_est[1:2]), exp(theta_est[3:8]), theta_est[9:12])
}

# Function to calculate the negative loglikelihood
negLoglik <- function(form, params_star, ssf) {

  # Extract parameters
  trans <- plogis(params_star[1:2])
  shape <- exp(params_star[3:4])
  scale <- exp(params_star[5:6])
  kappa <- exp(params_star[7:8])
  beta  <- matrix(params_star[9:12], nrow = 2, byrow = F)

  # Transition matrix and steady state
  t <- matrix(c(trans[1], 1 - trans[2], 1 - trans[1], trans[2]), nrow = 2)
  d <- solve(t(diag(2) - t + 1), c(1, 1))

  # Model matrix
  mat <- model.matrix(form, ssf)
  mat <- mat[, -1]

  # Step selection scores if in state 1
  score <- as.vector(exp(mat %*% beta[1, ]))
  probs_st1 <- aggregate(score, by = list(ssf$step_id), FUN = function(x) {
    x[1] / sum(x)
  })$x

  # Step selection scores if in state 2
  score <- as.vector(exp(mat %*% beta[2, ]))
  probs_st2 <- aggregate(score, by = list(ssf$step_id), FUN = function(x) {
    x[1] / sum(x)
  })$x

  # Put them together
  probs_st <- cbind(probs_st1, probs_st2)

  # Step length probabilities
  probs_sl <- cbind(
      dgamma(ssf$sl[ssf$case == 1], shape = shape[1], scale = scale[1])
    , dgamma(ssf$sl[ssf$case == 1], shape = shape[2], scale = scale[2])
  )

  # Turning angle probabilities
  probs_ta <- cbind(
      dvonmises(ssf$ta[ssf$case == 1], mu = 0, kappa = kappa[1])
    , dvonmises(ssf$ta[ssf$case == 1], mu = 0, kappa = kappa[2])
  )

  # Put all probabilities together
  probs <- probs_st * probs_sl * probs_ta

  # Forward algorithm
  foo <- d %*% diag(probs[1, ])
  l <- log(sum(foo))       # to avoid numerical issues
  phi <- foo / sum(foo)    # to avoid numerical issues
  for (i in 2:nrow(probs)) {
    foo <- phi %*% t %*% diag(probs[i, ])
    l <- l + log(sum(foo)) # to avoid numerical issues
    phi <- foo / sum(foo)  # to avoid numerical issues
  }

  # Return the negative log-likelihood
  return(-l)

}

# Function to decode the states
viterbi <- function(form, theta, ssf) {

  # Transition matrix and steady state
  t <- matrix(c(theta[1], 1 - theta[2], 1 - theta[1], theta[2]), nrow = 2)
  d <- solve(t(diag(2) - t + 1), c(1, 1))

  # Distributional parameters for the two states
  shape <- theta[3:4]
  scale <- theta[5:6]
  kappa <- theta[7:8]
  beta  <- matrix(theta[9:12], nrow = 2, byrow = F)

  # Model matrix
  mat <- model.matrix(form, ssf)
  mat <- mat[, -1]

  # Step selection scores if in state 1
  score <- as.vector(exp(mat %*% beta[1, ]))
  probs_st1 <- aggregate(score, by = list(ssf$step_id), FUN = function(x) {
    x[1] / sum(x)
  })$x

  # Step selection scores if in state 2
  score <- as.vector(exp(mat %*% beta[2, ]))
  probs_st2 <- aggregate(score, by = list(ssf$step_id), FUN = function(x) {
    x[1] / sum(x)
  })$x

  # Put them together
  probs_st <- cbind(probs_st1, probs_st2)

  # Step length probabilities
  probs_sl <- cbind(
      dgamma(ssf$sl[ssf$case == 1], shape = shape[1], scale = scale[1])
    , dgamma(ssf$sl[ssf$case == 1], shape = shape[2], scale = scale[2])
  )

  # Turning angle probabilities
  probs_ta <- cbind(
      dvonmises(ssf$ta[ssf$case == 1], mu = 0, kappa = kappa[1])
    , dvonmises(ssf$ta[ssf$case == 1], mu = 0, kappa = kappa[2])
  )

  # Put all probabilities together
  probs <- probs_st * probs_sl * probs_ta

  # Backwards algorithm
  n <- nrow(probs)
  xi <- matrix(0, n, 2)
  foo <- d * probs[1, ]
  xi[1, ] <- foo / sum(foo)
  for (i in 2:n) {
    foo <- apply(xi[i - 1, ] * t, 2, max) * probs[i, ]
    xi[i, ] <- foo / sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for (i in (n - 1):1) {
    iv[i] <- which.max(t[, iv[i + 1]] * xi[i, ])
  }
  return(iv)
}

# Function to run through the different models
runAnalysis <- function(sim, covars, multicore = F) {

  # Prepare dataframe to fill
  results <- tibble(
      Parameter     = names(params_true)
    , Value         = params_true
    , StatesKnown   = rep(NA, length(Parameter))
    , StatesIgnored = rep(NA, length(Parameter))
    , Nlm           = rep(NA, length(Parameter))
    , Optim         = rep(NA, length(Parameter))
  )

  # If we know the states, we can analyse the data separately
  for (i in 1:2) {

    # Fit state dependent distributions
    cat("Running sub-analysis:", i, "/ 3 \n")
    dist_sl <- fit_distr(sim$sl[sim$state == i], dist_name = "gamma")$params
    dist_ta <- fit_distr(sim$ta[sim$state == i], dist_name = "vonmises")$params
    dists <- list(sl = dist_sl, ta = dist_ta)

    # Prepare random steps and extract covariates
    ssf <- computeSSF(sim[sim$state == i, ], n_rsteps, dists)
    ssf <- computeCovars(ssf, covars, multicore = multicore)

    # Run model
    mod <- fit_clogit(case ~ elev + dist + strata(step_id), data = ssf)
    coefs <- coef(mod)
    if (i == 1) {
        results$StatesKnown[results$Parameter == "shape1"] <- dists$sl$shape
        results$StatesKnown[results$Parameter == "scale1"] <- dists$sl$scale
        results$StatesKnown[results$Parameter == "kappa1"] <- dists$ta$kappa
        results$StatesKnown[results$Parameter == "beta_elev1"] <- coefs["elev"]
        results$StatesKnown[results$Parameter == "beta_dist1"] <- coefs["dist"]
      } else {
        results$StatesKnown[results$Parameter == "shape2"] <- dists$sl$shape
        results$StatesKnown[results$Parameter == "scale2"] <- dists$sl$scale
        results$StatesKnown[results$Parameter == "kappa2"] <- dists$ta$kappa
        results$StatesKnown[results$Parameter == "beta_elev2"] <- coefs["elev"]
        results$StatesKnown[results$Parameter == "beta_dist2"] <- coefs["dist"]
    }
  }

  # If we dont know the states, we can't fit separate distributions
  cat("Running sub-analysis: 3 / 3 \n")
  dist_sl <- fit_distr(sim$sl, dist_name = "gamma")$params
  dist_ta <- fit_distr(sim$ta, dist_name = "vonmises")$params
  dists <- list(sl = dist_sl, ta = dist_ta)

  # Prepare random steps and extract covariates
  ssf <- computeSSF(sim, n_rsteps, dists)
  ssf <- computeCovars(ssf, covars, multicore = multicore)

  # Run model
  mod <- fit_clogit(case ~ elev + dist + strata(step_id), data = ssf)
  coefs <- coef(mod)

  # Add to results table
  results$StatesIgnored[results$Parameter %in% c("shape1", "shape2")] <- dists$sl$shape
  results$StatesIgnored[results$Parameter %in% c("scale1", "scale2")] <- dists$sl$scale
  results$StatesIgnored[results$Parameter %in% c("kappa1", "kappa2")] <- dists$ta$kappa
  results$StatesIgnored[results$Parameter %in% c("beta_elev1", "beta_elev2")] <- coefs["elev"]
  results$StatesIgnored[results$Parameter %in% c("beta_dist1", "beta_dist2")] <- coefs["dist"]

  # Sample initial parameter vector
  cat("Optimizing ML. This may take a moment...\n")
  params_star <- sampleParams(perfect = F)

  # Optimize using nlm
  ml_nlm <- suppressWarnings(nlm(negLoglik, params_star, form = formula, ssf = ssf))
  ml_nlm <- backtransformParams(ml_nlm$est)

  # Apply the viterbi algorithm to decode states and compute how often we were
  # right.
  states    <- viterbi(formula, ml_nlm, ssf[ssf$case == 1, ])
  n_correct <- sum(states == ssf$state[ssf$case == 1])

  # If we have more than half correct, we can leave the labels, otherwise we
  # should switch them
  if (n_correct < nrow(ssf[ssf$case == 1, ]) / 2) {
    ml_nlm <- ml_nlm[c(2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11)]
  }
  results$Nlm <- ml_nlm

  # Optimize using optim
  ml_optim <- suppressWarnings(optim(par = params_star, fn = negLoglik, form = formula, ssf = ssf))
  ml_optim <- backtransformParams(ml_optim$par)

  # Apply the viterbi algorithm to decode states and compute how often we were
  # right.
  states    <- viterbi(formula, ml_optim, ssf[ssf$case == 1, ])
  n_correct <- sum(states == ssf$state[ssf$case == 1])

  # If we have more than half correct, we can leave the labels, otherwise we
  # should switch them
  if (n_correct < nrow(ssf[ssf$case == 1, ]) / 2) {
    ml_optim <- ml_optim[c(2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11)]
  }
  results$Optim <- ml_optim

  # Return
  return(results)
}

################################################################################
#### Analysis
################################################################################
# Reload the design
design <- read_rds("03_Data/Design.rds")

# I want to replicate the analysis for each individual
design <- expand_grid(design, Replicate = 1:10)

# Prepare filenames for model results
design$FilenameModelResults <- paste0("03_Data/ModelResults/Results_ID", sprintf("%03d", design$ID), "_R", sprintf("%03d", design$Replicate), ".rds")

# Create the respective directory
dir.create("03_Data/ModelResults", showWarnings = F)

# Write the design to file
write_rds(select(design, ID, FilenameModelResults), "03_Data/Analysis.rds")

# Loop through the design and run analysis on the data
lapply(
    X   = 1:nrow(design)
  , FUN = function(x) {

  # Load covariates and movements
  covars          <- stack(design$FilenameCovariates[x])
  covars          <- readAll(covars)
  sim             <- read_rds(design$FilenameMovement[x])
  filename_results <- design$FilenameModelResults[x]

  # Run models (rerun if it fails, but not more than 5 times in total)
  if (!file.exists(filename_results)) {
    success <- F
    trials <- 1
    cat("Running analyses: Iteration", x, "/", nrow(design), "\n")
    while (!success & trials <= 5) {
      results <- tryCatch(runAnalysis(sim, covars, multicore = T)
        , error = function(e) {return(e)}
      )
      success <- is.data.frame(results)
      trials <- trials + 1
    }
    write_rds(results, filename_results)
  }
  return(NULL)
})

# Load all results
design$Results <- lapply(design$FilenameModelResults, read_rds)
design$Results[[14]]
hi <- design %>% select(ID, Replicate, Results) %>% unnest(Results)
subset(hi, Parameter == "beta_dist1") %>%
  pivot_longer(StatesKnown:Optim, names_to = "Method", values_to = "Estimate") %>%
  group_by(ID, Method) %>%
  summarize(Estimate = mean(Estimate)) %>%
  group_by(Method) %>%
  summarize(Estimate = mean(Estimate))

test <- do.call(rbind, design$Results)
test <- subset(test, Parameter %in% c("beta_elev1", "beta_elev2", "beta_dist1", "beta_dist2"))
test %>%
  subset(Parameter == "beta_elev2") %>%
  pull(StatesKnown) %>%
  mean()
