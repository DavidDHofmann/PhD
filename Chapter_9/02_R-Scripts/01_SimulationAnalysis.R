################################################################################
#### Simulation and Analysis of Hidden Markov Data with Step Selection Functions
################################################################################
# Clear R's brain
rm(list = ls())

# Suppress scientific notation
options(scipen = 999)

# Load required packages
library(tidyverse)     # For data wrangling
library(NLMR)          # To simulate covariates
library(raster)        # To handle spatial data
library(amt)           # To fit distributions
library(pbmcapply)     # To run stuff in parallel

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
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_9"
setwd(wd)

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Simulation Parameters
################################################################################
# Simulation setup
n         <- 300   # Resolution of the covariate layers
n_rsteps  <- 25    # Number of random steps to be generated
n_steps   <- 250   # Number of consecutive steps to be simulated
stop      <- F     # Should the simulation terminate at boundaries?

# State-dependent parameters
formula   <- ~ elev + dist
beta_elev <- c(-1, 2)
beta_dist <- c(-20, -20)
beta      <- cbind(beta_elev, beta_dist)
shape     <- c(1, 2)
scale     <- c(1, 4)
kappa     <- c(0.2, 0.9)
trans     <- matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2)

# Visualize the parametric distributions
sl <- tibble(
    x  = seq(0, 20, 0.01)
  , y1 = dgamma(x, shape = shape[1], scale = scale[1])
  , y2 = dgamma(x, shape = shape[2], scale = scale[2])
) %>% pivot_longer(y1:y2)
ta <- tibble(
    x  = seq(-pi, pi, 0.01)
  , y1 = dvonmises(x, mu = 0, kappa = kappa[1])
  , y2 = dvonmises(x, mu = 0, kappa = kappa[2])
) %>% pivot_longer(y1:y2)
ggplot(sl, aes(x = x, y = value, col = name)) + geom_line() + theme_minimal()
ggplot(ta, aes(x = x, y = value, col = name)) + geom_line() + theme_minimal()

# Put all relevant parameters into a single vector
params_true <- c(trans[1, 1], trans[2, 2], shape, scale, kappa, beta)
names(params_true) <- c("t11", "t22", "shape1", "shape2", "scale1", "scale2"
  , "kappa1", "kappa2", "beta_elev1", "beta_elev2", "beta_dist1", "beta_dist2")
print(params_true)

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

################################################################################
#### Functions
################################################################################
# Function to simulate covariates
simCovars <- function() {

  # Generate a point of attraction at the center of the study area and compute the
  # distance to it
  r <- raster(nrows = n, ncols = n, xmn = 0, xmx = n, ymn = 0, ymx = n)
  cent <- SpatialPoints(t(c(n / 2, n / 2)))
  dist <- distanceFromPoints(r, cent)
  dist <- (dist - cellStats(dist, min)) / (cellStats(dist, max) - cellStats(dist, min))

  # And a continuous covariate layer
  elev <- nlm_gaussianfield(ncol = n, nrow = n)

  # Put the covariates together into a stack
  covs <- stack(dist, elev)
  names(covs) <- c("dist", "elev")

  # Return
  return(covs)

}

# Function to simulate a single trajectory
simMove <- function(
      xy        = NULL    # Source point (in matrix form -> n * 2)
    , initstate = NULL    # Initial state of the animal
    , covars    = NULL    # Stack of covariate layers (names need to match formula!)
    , formula   = NULL    # Model formula used to predict selection score
    , prefs     = NULL    # Preferences used to predict selection score
    , shape     = NULL    # Shape parameters
    , scale     = NULL    # Scale parameters
    , kappa     = NULL    # Concentration parameters
    , trans     = NULL    # Matrix of transition probabilities
    , ext       = NULL    # Extent on which animals are allowed to move
    , n_steps   = 10      # Number of steps simulated
    , n_rsteps  = 25      # Number of random steps proposed at each step
    , stop      = TRUE    # Should the simulation stop at boundaries?
    , progress  = F       # Should a progress bar be shown?
  ) {

  # # For testing only
  # xy        <- cbind(150, 150)
  # covars    <- covars
  # formula   <- ~ elev + dist
  # prefs     <- beta
  # shape     <- shape
  # scale     <- scale
  # kappa     <- kappa
  # trans     <- trans
  # initstate <- sample(1:2, 1)
  # n_rsteps  <- 25
  # n_steps   <- 200
  # ext       <- ext
  # stop      <- F

  # Create a new dataframe based on the source point. Note that we draw random
  # turning angles to start off
  track <- data.frame(
      x     = c(NA, xy[, 1])
    , y     = c(NA, xy[, 2])
    , absta = c(runif(1, min = 0, max = 2 * pi), NA)
    , ta    = NA
    , sl    = NA
    , state = initstate
  )

  # Simulate random steps
  if (progress) {
    pb <- txtProgressBar(min = 0, max = n_steps, style = 3)
  }
  for (i in 2:(n_steps + 1)) {

    # For testing only
    # i <- 2

    # Sample a new state
    track$state[i] <- sample(c(1, 2), size = 1, prob = trans[track$state[i - 1], ])

    # Draw random turning angles
    ta_new <- rvonmises(n_rsteps
      , mu    = 0
      , kappa = kappa[track$state[i]]
    )

    # Draw random step lengths
    sl_new <- rgamma(n_rsteps
      , shape = shape[track$state[i]]
      , scale = scale[track$state[i]]
    )

    # Make sure that the steps cover at least a minimal distance (this is
    # relevant if we need to compute the log of it)
    sl_new[sl_new < 0.0001] <- 0.0001

    # Put the step lengths and turning angles into a new dataframe. These are
    # our proposed random steps.
    rand <- data.frame(
        absta  = track$absta[i - 1] + ta_new
      , ta     = ta_new
      , sl     = sl_new
    )

    # We need to make sure that the absolute turning angle ranges from 0 to 2 *
    # pi
    rand$absta[rand$absta > 2 * pi] <-
      rand$absta[rand$absta > 2 * pi] - 2 * pi
    rand$absta[rand$absta < 0] <-
      rand$absta[rand$absta < 0] + 2 * pi

    # Calculate new endpoints
    rand$x <- track$x[i] + sin(rand$absta) * rand$sl
    rand$y <- track$y[i] + cos(rand$absta) * rand$sl

    # Create spatial points from endpoints
    coordinates(rand) <- c("x", "y")

    # Depending on the answer in the beginning, the loop breaks if one of the
    # new coordinates is outside the map boundaries
    if (stop) {
        if (nrow(rand[ext, ]) != n_rsteps) {
          break
        }
      } else {
        rand <- rand[ext, ]
    }

    # Coerce back to regular dataframe
    rand <- as.data.frame(rand)

    # Prepare a "line" for each random step. We first need the coordinates of
    # the steps for this
    begincoords <- track[i, c("x", "y")]
    endcoords   <- rand[, c("x", "y")]

    # Interpolate coordinates and extract covariates
    extracted <- sapply(1:nrow(endcoords), function(x) {
      line <- interpolatePoints(
          x1 = begincoords[1, 1]
        , x2 = endcoords[x, 1]
        , y1 = begincoords[1, 2]
        , y2 = endcoords[x, 2]
        , by = 0.1
      )
      extr <- raster::extract(covars, line)
      extr <- colMeans(extr)
      return(extr)
    })

    # Bind with data on random steps
    rand <- cbind(rand, t(extracted))

    # Calculate cos_ta and log_sl
    rand$cos_ta <- cos(rand$ta)
    rand$log_sl <- log(rand$sl)

    # Prepare model matrix (and remove intercept)
    mat <- model.matrix(formula, rand)
    mat <- mat[ , -1]

    # Calculate selection scores
    score <- exp(mat %*% prefs[track$state[i], ])

    # Convert scores to probabilities
    probs <- score / sum(score)

    # Sample one of the steps based on predicted probabilities
    rand <- rand[sample(nrow(rand), 1, prob = probs), ]

    # Add the step to our track
    track$absta[i] <- rand$absta
    track$ta[i] <- rand$ta
    track$sl[i] <- rand$sl
    track[i + 1, "x"] <- rand$x
    track[i + 1, "y"] <- rand$y
    if (progress) {
      setTxtProgressBar(pb, i)
    }
  }

  # Add "to" coordinates
  track$x_to <- lead(track$x)
  track$y_to <- lead(track$y)

  # Return track, yet remove initial pseudo-fix
  track <- na.omit(track)
  track$step_number <- 1:nrow(track)
  return(track)
}

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
  params_star <- sampleParams(perfect = F)

  # Optimize using nlm
  ml_nlm <- suppressWarnings(nlm(negLoglik, params_star, form = formula, ssf = ssf))
  ml_nlm <- backtransformParams(ml_nlm$est)
  results$Nlm <- ml_nlm

  # Optimize using optim
  ml_optim <- suppressWarnings(optim(par = params_star, fn = negLoglik, form = formula, ssf = ssf))
  ml_optim <- backtransformParams(ml_optim$par)
  results$Optim <- ml_optim

  # Return
  return(results)
}

# ################################################################################
# #### Verifying the Functions Work
# ################################################################################
# # Simulate covariates
# covars <- simCovars()
#
# # Simulate movement
# sim <- simMove(
#     xy        = matrix(runif(2, xmin(ext2), xmax(ext2)), ncol = 2)
#   , covars    = covars
#   , formula   = formula
#   , prefs     = beta
#   , shape     = shape
#   , scale     = scale
#   , kappa     = kappa
#   , trans     = trans
#   , initstate = sample(1:2, 1)
#   , n_rsteps  = n_rsteps
#   , n_steps   = n_steps
#   , stop      = stop
#   , ext       = ext
#   , progress  = T
# )
# sim$ID <- 1
# sim$step_id <- 1:nrow(sim)
#
# # Run models
# runAnalysis(sim, multicore = T)

################################################################################
#### Replicated Analysis
################################################################################
# Design
design <- expand_grid(
    ID        = 1:10
  , Replicate = 1:10
)

# Loop
design$Results <- pbmclapply(
    X                  = 1:nrow(design)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {

  # Simulate covariates
  covars <- simCovars()

  # Simulate movement
  sim <- simMove(
      xy        = matrix(runif(2, xmin(ext2), xmax(ext2)), ncol = 2)
    , covars    = covars
    , formula   = formula
    , prefs     = beta
    , shape     = shape
    , scale     = scale
    , kappa     = kappa
    , trans     = trans
    , initstate = sample(1:2, 1)
    , n_rsteps  = n_rsteps
    , n_steps   = n_steps
    , stop      = stop
    , ext       = ext
    , progress  = F
  )
  sim$ID <- 1
  sim$step_id <- 1:nrow(sim)

  # Run models (rerun if it fails, but not more than 5 times in total)
  success <- F
  trials <- 1
  while (!success & trials <= 5) {
    results <- tryCatch(runAnalysis(sim, covars, multicore = F)
      , error = function(e) {return(e)}
    )
    success <- is.data.frame(results)
    trials <- trials + 1
  }
  return(results)
})

design$Results[[1]]
print(design, n = 100)
unnest(design, Results)
