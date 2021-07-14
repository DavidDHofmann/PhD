################################################################################
#### Simulating Movement Data with Known Preferences
################################################################################
# Description: Use ISSF analysis to simulate movement trajectories with known
# preferences

# Clear R's brain
rm(list = ls())

# Load required packages
library(RandomFields)  # To simulate covariates
library(raster)        # To handle spatial data
library(Rcpp)          # For faster point interpolation
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling

# Function to interpolate between spatial points
sourceCpp("/home/david/ownCloud/Dokumente/Bibliothek/Wissen/R-Scripts/interpolatepoints.cpp")

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

################################################################################
#### Simulate Covariates
################################################################################
# Generate random covariate (elevation)
elev <- RMexp(var = 5, scale = 10) +
  RMnugget(var = 1) +
  RMtrend(mean = 0)
elev <- RFsimulate(elev, x = 1:300, y = 1:300)
elev <- raster(elev)

# Add a point of attraction
center <- coordinates(elev)
center <- colMeans(center)
center <- SpatialPoints(t(center))

# Calculate distance to center
dist <- distanceFromPoints(elev, center)

# Normalize both covariates
elev <- (elev - cellStats(elev, mean)) / cellStats(elev, sd)
dist <- (dist - cellStats(dist, mean)) / cellStats(dist, sd)

# Put covariate layers into a stack
covars <- stack(elev, dist)
names(covars) <- c("elev", "dist")

# Let's also specify an extent on which the animals are allowed to move
ext <- extent(c(50, 250, 50, 250))
ext <- as(ext, "SpatialPolygons")

# Visualize all covariates + extent
rasterVis::levelplot(covars, at = seq(-5, 5, length = 50))

################################################################################
#### Function To Simulate Movement
################################################################################
# Function to simulate movement
move <- function(
      xy       = NULL
    , covars   = NULL
    , ext      = NULL
    , prefs    = NULL
    , sl_dist  = NULL
    , n_steps  = 10
    , n_rsteps = 25
    , stop     = TRUE
  ){

  # If no extent is provided, use the extent of the covariates
  if (is.null(ext)){
    ext <- extent(covars)
  }
  ext <- as(ext, "SpatialPolygons")

  # create a new dataframe based on the source point. Note that we draw random
  # turning angles to start off
  track <- data.frame(
      x     = c(NA, xy[, 1])
    , y     = c(NA, xy[, 2])
    , absta = c(runif(1, min = 0, max = 2 * pi), NA)
    , ta    = c(runif(1, min = 0, max = 2 * pi), NA)
    , sl    = NA
  )

  # Simulate random steps
  for (i in 2:n_steps){

    # Prepare an empty list in which we can store the random steps
    rand <- list()

    # Draw random turning angles
    ta_new <- runif(n_rsteps
      , min = -pi
      , max = +pi
    )

    # Draw random step lengths
    sl_new <- rgamma(n_rsteps
      , shape = sl_dist["shape"]
      , scale = sl_dist["scale"]
    )

    # Make sure that the steps cover at least a minimal distance
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
    mat <- model.matrix(~ elev + dist + sl + log_sl + cos_ta, rand)
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
#### Single Trajectory
################################################################################
# Simulation Parameters
stop      <- F
n_rsteps  <- 25
n_steps   <- 200
sl_dist   <- c(shape = 3, scale = 1)
prefs     <- c(
    elev   = 0.5    # Preference towards elevation
  , dist   = -10.0   # Preference for distance to point of attraction
  , sl     = 0      # Preference for step length
  , log_sl = 0      # Preference for step length
  , cos_ta = 0      # Preference for turning angle
)

# # Prepare a function that resamples preferences from given distributions
# randomPrefs <- function(){
#   prefs <- c(
#       elev   = rnorm(n = 1, mean = 0.5, sd = 0.2)
#     , dist   = rnorm(n = 1, mean = -0.5, sd = 0.5)
#     , sl     = 0
#     , log_sl = 0
#     , cos_ta = 0
#   )
#   return(prefs)
# }
#
# # Try it
# randomPrefs()

# Simulate a test trajectory
sim <- move(
    xy       = matrix(runif(2, 50, 150), ncol = 2)
  , covars   = covars
  , ext      = ext
  , stop     = stop
  , n_rsteps = n_rsteps
  , n_steps  = n_steps
  , sl_dist  = sl_dist
  , prefs    = prefs
)

# Add an ID for the individual and an ID for each step
sim$ID <- 1
sim$step_id <- 1:nrow(sim)

# Visualize the simulation
plot(covars[[1]])
plot(ext, add = T, border = "red", lwd = 2)
points(sim$y ~ sim$x, type = "o", pch = 16, cex = 0.5)
plot(center, add = T, col = "red", pch = 20)

################################################################################
#### Multiple Trajectories
################################################################################
# Number of simulated individuals
ndisp <- 10

# Simulate track for each individual
sims <- pbmclapply(
    X                  = 1:ndisp
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x){

  # Simulate trajectory
  sim <- move(
      xy       = matrix(runif(2, 50, 150), ncol = 2)
    , covars   = covars
    , ext      = ext
    , stop     = stop
    , n_rsteps = n_rsteps
    , n_steps  = n_steps
    , sl_dist  = sl_dist
    , prefs    = prefs
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

# Visualize them
ids <- unique(sims$ID)
col <- rainbow(length(ids))
plot(covars[[1]])
plot(ext, add = T, border = "red", lwd = 2)
for (i in 1:length(unique(sims$ID))){
  sub <- subset(sims, ID == ids[i])
  points(sub$y ~ sub$x, col = col[i], pch = 16, cex = 0.4)
  lines(sub$y ~ sub$x, col = col[i], lwd = 0.3)
}

# In reality, we don't observe step lengths, turning angles etc., but we derive
# them from xy coordinates. Hence, let's that we only observed xy data + ID
obs <- dplyr::select(sims, x, y, ID, step_number, step_id)

# This is the type of data that we would observe in reality. Let's store it to
# file. Also store the simulated spatial layers to file
write_csv(obs, "03_Data/01_RawData/ObservedMovements.csv")
writeRaster(covars, "03_Data/01_RawData/CovariateLayers.tif")
