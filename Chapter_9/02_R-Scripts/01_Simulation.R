################################################################################
#### Simulation of Movement using Two-State HMM & Step Selection Functions
################################################################################
# Clear R's brain
rm(list = ls())

# Suppress scientific notation
options(scipen = 999)

# Load required packages
library(tidyverse)     # For data wrangling
library(NLMR)          # To simulate covariates
library(raster)        # To handle spatial data
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

# Store them
write_rds(params_true, "03_Data/SimulationParameters.rds")

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

################################################################################
#### Verifying the Functions Work
################################################################################
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
  , progress  = T
)
sim$ID <- 1
sim$step_id <- 1:nrow(sim)

# Visualize
ggplot() +
  geom_raster(data = as.data.frame(covars, xy = T), aes(x = x, y = y, fill = elev)) +
  geom_path(data = sim, aes(x = x, y = y), size = 0.5) +
  geom_point(data = sim, aes(x = x, y = y, col = as.factor(state), group = ID), size = 1) +
  scale_fill_viridis_c() +
  coord_sf() +
  theme_minimal()

################################################################################
#### Simulate
################################################################################
# Design
design <- expand_grid(ID = 1:100)

# Indicate filenames for the simulated covariate layers and the simulated
# movement
design$FilenameCovariates <- paste0("03_Data/Simulations/Covariates_", sprintf("%02d", design$ID), ".grd")
design$FilenameMovement   <- paste0("03_Data/Simulations/Movement_", sprintf("%02d", design$ID), ".rds")

# Create the respective directory
dir.create("03_Data/Simulations", showWarnings = F)

# Write the design to file
write_rds(design, "03_Data/Design.rds")

# Loop through the design and simulate data
pbmclapply(
    X                  = 1:nrow(design)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {

  # Extract relevant information
  filename_covariates <- design$FilenameCovariates[x]
  filename_simulation <- design$FilenameMovement[x]

  # Simulate covariates
  if (file.exists(filename_covariates)) {
      covars <- stack(filename_covariates)
    } else {
      covars <- simCovars()
  }

  # Simulate movement
  if (!file.exists(filename_simulation)) {
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
    sim$ID <- x
    sim$step_id <- ((x - 1) * n_steps) + (1:n_steps)
    write_rds(sim, design$FilenameMovement[x])
  }

  # Extend the covariates
  covars <- extendRaster(covars, ext3)

  # Store the simulations to file
  writeRaster(covars, design$FilenameCovariates[x], overwrite = T)

  # Return nothing
  return(NULL)
})
