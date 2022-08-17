################################################################################
#### Simulating Movement Data with Known Preferences
################################################################################
# Description: Use integrated step selection functions to simulate movement
# trajectories with known preferences

# Clear R's brain
rm(list = ls())

# Load required packages
library(NLMR)          # To simulate covariates
library(raster)        # To handle spatial data
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(sf)            # For plotting spatial features

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")
# setwd("C:/Users/david/switchdrive/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Simulate Covariate Layers
################################################################################
# Specify resolution of covariates
n <- 300

# Generate a point of attraction at the center of the study area and compute the
# distance to it
r <- raster(nrows = n, ncols = n, xmn = 0, xmx = n, ymn = 0, ymx = n)
cent <- SpatialPoints(t(c(n / 2, n / 2)))
dist <- distanceFromPoints(r, cent)
dist <- (dist - cellStats(dist, min)) / (cellStats(dist, max) - cellStats(dist, min))
plot(dist)

# Function to simulate forest cover
simForest <- function(autocorr_range = 10, proportion = 0.5) {
  forest <- nlm_gaussianfield(
      ncol           = n
    , nrow           = n
    , autocorr_range = autocorr_range
    , mag_var        = 1
    , nug            = 0
  )
  cutoff <- quantile(forest, 1 - proportion)
  forest <- forest > cutoff
  forest <- focal(forest
    , w        = matrix(rep(1, 9), nrow = 3)
    , fun      = function(x) {sum(x) > 4.5}
    , pad      = T
    , padValue = 0
  )
  return(forest)
}

# Function to simulate elevation
simElev <- function(autocorr_range = 10) {
  elev <- nlm_gaussianfield(
      ncol           = n
    , nrow           = n
    , autocorr_range = autocorr_range
    , mag_var        = 1
    , nug            = 0
  )
  return(elev)
}

# Function to simulate the different covariates and return them as a stack
simCovars <- function(autocorr_range = 10) {

  # Simulate covariates
  forest <- simForest(autocorr_range)
  elev   <- simElev(autocorr_range)

  # Put them into a stack
  covars <- stack(forest, elev, dist)
  names(covars) <- c("forest", "elev", "dist")
  return(covars)
}

# Let's try it
covars <- simCovars(autocorr_range = 10)
plot(covars)

# Let's specify the extent on which animals are allowed to move
ext <- extent(0, n, 0, n)
ext <- as(ext, "SpatialPolygons")

# Let's also specify an extent within which individuals will be released
ext2 <- extent(ext) - 100
ext2 <- as(ext2, "SpatialPolygons")

# Plot the covariates
as.data.frame(covars, xy = T) %>%
  gather(key = covariate, value = value, 3:5) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = st_as_sf(ext2), inherit.aes = F, fill = NA, col = "white", lty = 2) +
    scale_fill_viridis_c(option = "viridis") +
    coord_sf() +
    theme_minimal() +
    facet_wrap("covariate") +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

################################################################################
#### Simulate Single Trajectory
################################################################################
# Simulation Parameters
formula <- ~ forest + elev + dist
prefs     <- c(-1, 0.5, -20)
sl_dist   <- list(name = "gamma", params = list(shape = 3, scale = 1))
ta_dist   <- list(name = "vonmises", params = list(kappa = 0.5, mu = 0))
n_rsteps  <- 10
n_steps   <- 100
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
ggplot() +
  geom_raster(data = as.data.frame(covars[[2]], xy = T), aes(x = x, y = y, fill = elev)) +
  geom_path(data = sim, aes(x = x, y = y), size = 0.2) +
  geom_point(data = sim, aes(x = x, y = y), size = 0.2) +
  geom_point(data = sim[1, ], aes(x = x, y = y), col = "green", size = 2) +
  geom_point(data = sim[nrow(sim), ], aes(x = x, y = y), col = "red", size = 2) +
  geom_sf(data = st_as_sf(ext2), col = "white", fill = NA) +
  scale_fill_viridis_c() +
  coord_sf() +
  theme_minimal()

################################################################################
#### Multiple Trajectories
################################################################################
# Number of simulated individuals per treatment
ndisp <- 1000

# Generate a design matrix containing the different treatments through which we
# will loop
sims <- expand_grid(
    AutocorrRange = c(1, 10, 100)
  , ID            = 1:ndisp
) %>% mutate(ID = 1:n())

# Simulate covariates for the different treatments. Unfortunately, we can't do
# this in parallel, as RFoptions must not be called within parallel code, yet is
# always called when simulating a covariate layer.
cat("Simulating covariate layers...\n")
pb <- txtProgressBar(min = 0, max = nrow(sims), style = 3)
sims$Covariates <- lapply(1:nrow(sims), function(x) {
  covs <- simCovars(autocorr_range = sims$AutocorrRange[x])
  setTxtProgressBar(pb, x)
  return(covs)
})

# Simulate track for each individual
cat("Simulating movement...\n")
sims$Simulations <- pbmclapply(
    X                  = 1:nrow(sims)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {

  # Simulate trajectory
  sim <- move(
      xy       = matrix(runif(2, xmin(ext2), xmax(ext2)), ncol = 2)
    , covars   = sims$Covariates[[x]]
    , formula  = formula
    , prefs    = prefs
    , sl_dist  = sl_dist
    , ta_dist  = ta_dist
    , ext      = ext
    , n_steps  = n_steps
    , n_rsteps = n_rsteps
    , stop     = stop
  )

  # Keep only variables that are observed in reality
  sim <- dplyr::select(sim, x, y, step_number)

  # Return the simulation
  return(sim)
})

# Store the data to file (we'll separate the movement data from the covariates)
write_rds(dplyr::select(sims, -Covariates), "03_Data/01_RawData/SimulatedMovement.rds")
write_rds(dplyr::select(sims, -Simulations), "03_Data/01_RawData/SimulatedCovariates.rds")

# Visualize simulations
sims %>%
  dplyr::select(-Covariates) %>%
  unnest(Simulations) %>%
  ggplot(aes(x = x, y = y, col = as.factor(ID))) +
    geom_path(size = 0.1) +
    geom_point(size = 0.1) +
    geom_sf(data = st_as_sf(ext2), col = "black", fill = NA, inherit.aes = F) +
    coord_sf(xlim = c(0, n), ylim = c(0, n)) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_color_viridis_d()
