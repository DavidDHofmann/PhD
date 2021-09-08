################################################################################
#### Simulating Movement Data with Known Preferences
################################################################################
# Description: Use ISSF analysis to simulate movement trajectories with known
# preferences. We will simulate movement through the different landscapes

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(Rcpp)          # For faster point interpolation
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(sf)            # For plotting spatial features
library(rgeos)         # For manipulating spatial objects
library(ggpubr)        # To put plots together

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_6")

# Function to interpolate between spatial points
sourceCpp("02_R-Scripts/interpolatepoints.cpp")

# Load covariates and other necessary data from previous session
load("03_Data/Landscape.Rdata")

# Make sure the covariates are stored in memory
inMemory(covars$Layers[[1]])

################################################################################
#### Useful Functions
################################################################################
# Function to determine the pdf of a mixed von mises distribution
dvonmises <- function(x, kappa, mu){
  exp(kappa * cos(x - mu)) / (2 * pi * besselI(kappa, nu = 0))
}

# Function to randomly sample from a mixed von mises distribution
rvonmises <- function(n, kappa, mu, by = 0.01){
  x <- seq(-pi, +pi, by = by)
  probs <- dvonmises(x, kappa = kappa, mu = mu)
  random <- sample(x, size = n, prob = probs, replace = T)
  return(random)
}

################################################################################
#### Function To Simulate Movement
################################################################################
# Function to simulate movement
move <- function(
      xy       = NULL    # Source point (in matrix form -> n * 2)
    , covars   = NULL    # Stack of covariate layers (names need to match!)
    , formula  = NULL    # Model formula used to predict selection score
    , prefs    = NULL    # Preferences used to predict selection score
    , sl_dist  = NULL    # Parameters describing the step length distribution
    , ta_dist  = NULL    # Parameters describing the turning angle distribution
    , ext      = NULL    # Extent on which animals are allowed to move
    , n_steps  = 10      # Number of steps simulated
    , n_rsteps = 25      # Number of random steps proposed at each iteration
    , stop     = TRUE    # Should the simulation stop at boundaries?
  ){

  # If no extent is provided, use the extent of the covariates
  if (is.null(ext)){
    ext <- extent(covars)
  }
  ext <- as(ext, "SpatialPolygons")

  # Create a new dataframe based on the source point. Note that we draw random
  # turning angles to start off
  track <- data.frame(
      x     = c(NA, xy[, 1])
    , y     = c(NA, xy[, 2])
    , absta = c(runif(1, min = 0, max = 2 * pi), NA)
    , ta    = NA
    , sl    = NA
  )

  # Simulate random steps
  for (i in 2:n_steps){

    # Draw random turning angles
    ta_new <- rvonmises(n_rsteps
      , mu    = ta_dist$params$mu
      , kappa = ta_dist$params$kappa
    )

    # Draw random step lengths
    sl_new <- rgamma(n_rsteps
      , shape = sl_dist$params$shape
      , scale = sl_dist$params$scale
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
    if (stop){
      if (nrow(rand[ext, ]) != n_rsteps){
        break
      }
    } else {
      rand <- rand[ext, ]
    }

    # Coerce back to regular dataframe
    rand <- as.data.frame(rand, xy = T)
    rand$xy <- NULL

    # Prepare a "line" for each random step. We first need the coordinates of
    # the steps for this
    begincoords <- track[i, c("x", "y")]
    endcoords   <- rand[, c("x", "y")]

    # Interpolate coordinates and extract covariates
    extracted <- sapply(1:nrow(endcoords), function(x){
      line <- interpolatePointsC(
          x1 = begincoords[1, 1]
        , x2 = endcoords[x, 1]
        , y1 = begincoords[1, 2]
        , y2 = endcoords[x, 2]
        , by = 1
      )
      extr <- raster::extract(covars, line)
      extr <- colMeans(extr)
      return(extr)
    })

    # Bind with other data
    rand <- cbind(rand, t(extracted))

    # Calculate cos_ta and log_sl
    rand$cos_ta <- cos(rand$ta)
    rand$log_sl <- log(rand$sl)

    # Prepare model matrix (and remove intercept)
    mat <- model.matrix(formula, rand)
    mat <- mat[ , 2:ncol(mat)]

    # Calculate selection scores
    score <- exp(mat %*% prefs)

    # Convert scores to probabilities
    probs <- score / sum(score)

    # Keep only the step with the highest score
    rand <- rand[sample(nrow(rand), 1, prob = probs), ]

    # Add the step to our track
    track$absta[i] <- rand$absta
    track$ta[i] <- rand$ta
    track$sl[i] <- rand$sl
    track[i + 1, "x"] <- rand$x
    track[i + 1, "y"] <- rand$y
  }

  # Assign step numbers
  track$step_number <- 0:(nrow(track) - 1)

  # Return track, yet remove initial pseudo-fix
  return(track[-1, ])
}

################################################################################
#### Simulate Single Trajectory
################################################################################
# Simulation Parameters
formula <- ~ water + elev + dist + cos_ta
prefs     <- c(-2, 0.5, -5, 1)
sl_dist   <- list(name = "gamma", params = list(shape = 3, scale = 1))
ta_dist   <- list(name = "vonmises", params = list(kappa = 0, mu = 0))
n_rsteps  <- 25
n_steps   <- 100
stop      <- T

# Plot the selected step length and turning angle distributions
ta <- rvonmises(n = 1000, kappa = ta_dist$params$kappa, mu = ta_dist$params$mu)
sl <- rgamma(n = 1000, shape = sl_dist$params$shape, scale = sl_dist$params$scale)
hist(ta, col = "cornflowerblue", border = "white", main = "Turning Angle Distribution")
hist(sl, col = "cornflowerblue", border = "white", main = "Step Length Distribution")

# Simulate a test trajectory
sim <- move(
    xy       = coordinates(spsample(nps, type = "random", n = 1))
  , covars   = covars$Layers[[1]]
  , formula  = formula
  , prefs    = prefs
  , sl_dist  = sl_dist
  , ta_dist  = ta_dist
  , ext      = ext
  , n_steps  = n_steps
  , n_rsteps = n_rsteps
  , stop     = stop
)

# Check distribution of step lengths and turning angles again
hist(sim$ta, col = "cornflowerblue", border = "white", main = "Turning Angle Distribution")
hist(sim$sl, col = "cornflowerblue", border = "white", main = "Step Length Distribution")

# Visualize the simulation
plot(covars$Layers[[1]][[1]], asp = 1.03)
points(sim$y ~ sim$x, type = "o", pch = 16, cex = 0.2)
points(sim$y[1] ~ sim$x[1], type = "o", pch = 16, cex = 2, col = "green")

################################################################################
#### Multiple Trajectories
################################################################################
# Let's specify the design matrix
design <- expand_grid(
    Covariates = unique(covars$Type)
  , Replicate  = unique(covars$Replicate)
  , Repell     = c(T, F)
)

# Simulation Parameters
formula <- ~ water + elev + dist + cos_ta
prefs     <- c(-2, 0.5, -5, 1)
sl_dist   <- list(name = "gamma", params = list(shape = 3, scale = 1))
ta_dist   <- list(name = "vonmises", params = list(kappa = 0, mu = 0))
n_rsteps  <- 25
n_steps   <- 500
stop      <- T
ndisp     <- 100       # Per source area

# Sample starting points
design$SourcePoints <- lapply(1:nrow(design), function(y){
  pts <- lapply(1:length(nps), function(x){
    spsample(nps[x, ], type = "random", n = ndisp)
  }) %>% do.call(rbind, .)
})

# Let's also give each row a unique filename so that we can store the outputs to
# result occasionally
dir.create("03_Data/Simulation", showWarnings = F)
design$Filename <- paste0("03_Data/Simulation/Simulation_No", 1:nrow(design), ".Rdata")

# Loop through the design matrix and simulate individuals according to the
# design parameters
lapply(1:nrow(design), function(i){

  # Identify the needed covariate layer
  cov <- subset(covars, Type == design$Covariates[i] & Replicate == design$Replicate[i])
  cov <- cov$Layers[[1]]

  # Prepare its extent
  ext <- as(extent(cov), "SpatialPolygons")

  # Extract the source points
  pts <- design$SourcePoints[[i]]
  pts <- coordinates(pts)

  # Simulate track for each individual
  sims <- pbmclapply(
      X                  = 1:nrow(pts)
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x){

    # Simulate trajectory
    sim <- move(
        xy       = matrix(pts[x, ], ncol = 2)
      , covars   = cov
      , formula  = formula
      , prefs    = prefs
      , sl_dist  = sl_dist
      , ta_dist  = ta_dist
      , ext      = ext
      , n_steps  = n_steps
      , n_rsteps = n_rsteps
      , stop     = design$Repell[i]
    )

    # Assign unique simulation ID
    sim$ID <- x

    # Return the simulation
    return(sim)
  })

  # Bind simulations
  sims <- do.call(rbind, sims)
  rownames(sims) <- NULL

  # Assign unique ID to each step
  sims$step_id <- 1:nrow(sims)

  # Store the simulation to file
  save(sims, file = design$Filename[i])

})

# Remove the source points from the design matrix
design <- design %>% dplyr::select(-SourcePoints)

# Store the design
save(design, file = "03_Data/Design.Rdata")
