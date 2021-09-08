################################################################################
#### Simulating Movement Data with Known Preferences
################################################################################
# Description: Use ISSF analysis to simulate movement trajectories with known
# preferences

# Clear R's brain
rm(list = ls())

# Load required packages
library(NLMR)          # To simulate covariates
library(raster)        # To handle spatial data
library(Rcpp)          # For faster point interpolation
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(sf)            # For plotting spatial features
library(circular)      # For vonMises distribution

# Function to interpolate between spatial points
sourceCpp("/home/david/ownCloud/Dokumente/Bibliothek/Wissen/R-Scripts/interpolatepoints.cpp")

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

################################################################################
#### Simulate Covariates
################################################################################
# Set seed for reproducability
set.seed(12345)

# Specify resolution of covariates
n <- 300

# Simulate water cover
water <- nlm_randomcluster(
    ncol = n
  , nrow = n
  , p    = 0.5
  , ai   = c(0.8, 0.2)
)
plot(water)

# Simulate elevation
elev <- nlm_gaussianfield(
    ncol           = n
  , nrow           = n
  , autocorr_range = 10
  , mag_var        = 1
  , nug            = 0
)
plot(elev)

# Add a point of attraction at the center fo the study area
center <- coordinates(elev)
center <- colMeans(center)
center <- SpatialPoints(t(center))

# Calculate distance to center
dist <- distanceFromPoints(elev, center)

# Normalize covariates to range betwee 0 and 1 (water already is 0 to 1)
elev <- (elev - cellStats(elev, min)) / (cellStats(elev, max) - cellStats(elev, min))
dist <- (dist - cellStats(dist, min)) / (cellStats(dist, max) - cellStats(dist, min))

# Put covariate layers into a stack
covars <- stack(water, elev, dist)
names(covars) <- c("water", "elev", "dist")

# Let's also specify an extent on which the animals are allowed to move
ext <- extent(c(50, 250, 50, 250))
ext <- as(ext, "SpatialPolygons")

# Plot the covariates
as.data.frame(covars, xy = T) %>%
  gather(key = covariate, value = value, 3:5) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = st_as_sf(ext), inherit.aes = F, fill = NA, col = "white", lty = 2) +
    scale_fill_viridis_c(option = "viridis") +
    coord_sf() +
    theme_minimal() +
    facet_wrap("covariate") +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

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
    ta_new <- suppressWarnings(as.vector(rvonmises(n_rsteps
      , mu    = ta_dist$params$mu
      , kappa = ta_dist$params$kappa
    )))

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
#### Simulate Trajectory
################################################################################
# Simulation Parameters (you can play around with these as you wish, to see how
# results are influenced)
formula <- ~ water + elev + dist
prefs     <- c(-1, 0.5, -15)
sl_dist   <- list(name = "gamma", params = list(shape = 3, scale = 1))
ta_dist   <- list(name = "vonmises", params = list(kappa = 0, mu = 0))
n_rsteps  <- 25
n_steps   <- 500
stop      <- F

# Simulate a test trajectory
sim <- move(
    xy       = matrix(runif(2, 50, 250), ncol = 2)
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

# Add an ID for the individual and an ID for each step
sim$ID <- 1
sim$step_id <- 1:nrow(sim)

# Visualize the simulation
plot(covars[[1]])
plot(ext, add = T, border = "red", lwd = 2)
points(sim$y ~ sim$x, type = "o", pch = 16, cex = 0.2)
points(sim$y[1] ~ sim$x[1], type = "o", pch = 16, cex = 2, col = "green")
plot(center, add = T, col = "red", pch = 20)

################################################################################
#### Multiple Trajectories
################################################################################
# Number of simulated individuals
ndisp <- 100

# Simulate track for each individual
sims <- pbmclapply(
    X                  = 1:ndisp
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x){

  # Simulate trajectory
  sim <- move(
      xy       = matrix(runif(2, 50, 250), ncol = 2)
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
  points(sub$y ~ sub$x, col = col[i], pch = 16, cex = 0.2)
  lines(sub$y ~ sub$x, col = col[i], lwd = 0.3)
}

# In reality, we don't observe step lengths, turning angles etc., but we derive
# them from xy coordinates. Hence, let's that we only observed xy data + ID
obs <- dplyr::select(sims, x, y, ID, step_number, step_id)

# This is the type of data that we would observe in reality. Let's store it to
# file. Also store the simulated spatial layers to file
write_csv(obs, "03_Data/01_RawData/ObservedMovements.csv")
writeRaster(covars, "03_Data/01_RawData/CovariateLayers.grd", overwrite = T)
