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

# Load custom functions
source("02_R-Scripts/00_Functions.r")

################################################################################
#### Simulate Covariate Layers
################################################################################
# Set seed for reproducability
set.seed(12345)

# Specify resolution of covariates
n <- 300

# Simulate forest cover
forest <- nlm_randomcluster(
    ncol = n
  , nrow = n
  , p    = 0.5
  , ai   = c(0.8, 0.2)
)
plot(forest)

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

# Normalize covariates to range betwee 0 and 1 (forest already is 0 to 1)
elev <- (elev - cellStats(elev, min)) / (cellStats(elev, max) - cellStats(elev, min))
dist <- (dist - cellStats(dist, min)) / (cellStats(dist, max) - cellStats(dist, min))

# Put covariate layers into a stack
covars <- stack(forest, elev, dist)
names(covars) <- c("forest", "elev", "dist")

# Let's specify the extent on which animals are allowed to move
ext <- extent(covars)
ext <- as(ext, "SpatialPolygons")

# Let's also specify an extent within which individuals will be released
ext2 <- extent(c(50, 250, 50, 250))
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
#### Simulate Trajectory
################################################################################
# Simulation Parameters
formula <- ~ forest + elev + dist
prefs     <- c(-1, 0.5, -15)
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

# Visualize simulations
ggplot(sims, aes(x = x, y = y, col = as.factor(ID))) +
  geom_path(size = 0.1) +
  geom_point(size = 0.1) +
  geom_sf(data = st_as_sf(ext2), col = "red", fill = NA, inherit.aes = F) +
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_viridis_d()

# In reality, we don't observe step lengths, turning angles etc., but we derive
# them from xy coordinates. Hence, let's assume that we only observed xy data +
# ID
obs <- dplyr::select(sims, x, y, ID, step_number, step_id)

# This is the type of data that we would observe in reality. Let's store it to
# file. Also store the simulated spatial layers to file
write_csv(obs, "03_Data/01_RawData/ObservedMovements.csv")
writeRaster(covars, "03_Data/01_RawData/CovariateLayers.grd", overwrite = T)
