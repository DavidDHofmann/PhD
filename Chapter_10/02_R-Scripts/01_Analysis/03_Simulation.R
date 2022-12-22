################################################################################
#### Simulation
################################################################################
# Description: Simulation of dispersal

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_10")

# Load required packages
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(velox)          # For quick extraction
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(glmmTMB)        # To handle the movement model
library(Rcpp)           # To load C++ functions
library(pbmcapply)      # For multicore use with progress bar
library(mapview)        # For interactive plots

# Load custom functions
source("02_R-Scripts/00_Functions.R")
sourceCpp("02_R-Scripts/00_Functions.cpp")

################################################################################
#### Load Required Data
################################################################################
# Load the calibrated movement model
mod <- read_rds("03_Data/02_CleanData/MovementModel.rds")

# Load the scaling parameters
scaling <- read_rds("03_Data/02_CleanData/Scaling.rds")

# Load the fitted gamma distribution (for sampling step lengths)
sl_dist <- read_rds("03_Data/02_CleanData/GammaDistribution.rds")

# Load the required spatial layers layers
layers <- rast("03_Data/02_CleanData/SpatialLayers.tif")

# Extend the covariates artificially (fill added buffer with random values)
cat("Extending covariate layers...\n")
ext <- ext(layers) * 1.2
cov <- extendRaster(layers, ext)

# Apply our custom functions to make everything ready for the simulation
mod  <- prepareModel(mod)
cov  <- prepareCovars(layers)

# Define an extent and create a reference raster
ext <- ext(20.5, 26.5, -21.5, -17.5)
r <- rast(ext, resolution = 100 / 111000, crs = "epsg:4326")
# ext <- ext(23.2, 24.2, -20, -19)
# r <- rast(ext, resolution = 100 / 111000, crs = "epsg:4326")

# Load areas that we want to remove from the domain of "suitable" habitats
humans <- layers[["HumansBuff5000"]]
humans <- crop(humans, ext)
humans <- humans > 0
water  <- layers[["Water"]]
water  <- crop(water, ext)
major  <- vect("03_Data/02_CleanData/MajorWaters.shp")
water  <- mask(water, major, updatevalue = 1, inverse = T)

# Let's create a mask that we'll use to remove everything that is in one of the
# two layers
m <- water + humans
m <- m > 0
m <- as.polygons(m)
m <- m[m$Water == 1]
m <- buffer(m, width = 0)

################################################################################
#### Simulation Setup
################################################################################
# Number of simulated steps
n_steps <- 50

# How many dispersers do you want to simulate?
n_points <- 10

# How many random steps do you want to simulate per realized step?
n_rsteps <- 25

# Do you want to break the simulation of a track if it hits a boundary?
stop <- F

# What is the largest step possible in 4 hours (in meters)?
sl_max <- 35000

# Generate a layer of home ranges
hrs <- simHR(n = 400, m = m, r = r, buff = 10000, 50, 1800)
hrs <- as.polygons(hrs)
# mapview(as(hrs, "Spatial"))

# Function to sample locations from the different HRs
sampleLocations <- function(p) {
  locs <- lapply(1:length(p), function(x) {
    pts <- spatSample(p[x, ], size = 20)
    return(pts[1, ])
  })
  locs <- do.call(rbind, locs)
  return(locs)
}

################################################################################
#### Run Simulation
################################################################################
# Go through the design and run the simulation
cat("Simulating dispersal...\n")

# Generate 10 source points
pts <- spatSample(hrs, size = 20)
plot(hrs)
plot(pts, add = T, col = "red")

# Run the simulation for each source point
tracks <- pbmclapply(
    X                   = 1:length(pts)
  , mc.cores            = detectCores() - 1
  , ignore.interactive  = T
  , FUN                 = function(x) {
    sim <- suppressWarnings(
      disperse(
          source    = crds(pts[x, ])
        , covars    = cov
        , model     = mod
        , sl_dist   = sl_dist
        , sl_max    = sl_max
        , date      = as.POSIXct("2021-01-01 07:00:00", tz = "UTC")
        , n_rsteps  = n_rsteps
        , n_steps   = n_steps
        , scaling   = scaling
        , stop      = stop
      )
    )

    # Assign some more information
    sim$TrackID <- x

    # Create new columns containing the endpoints of each step
    sim$x_to <- lead(sim$x)
    sim$y_to <- lead(sim$y)
    sim <- dplyr::select(sim, x, y, x_to, y_to, everything())

    # Return the simulation
    return(sim)
})

# Identify closest track for each coordinate
tracks <- tracks %>%
  do.call(rbind, .) %>%
  subset(!is.na(x_to) & !is.na(y_to)) %>%
  nest(Data = -Timestamp) %>%
  mutate(Data = map(Data, function(x) {
      d       <- distance(cbind(x$x_to, x$y_to), lonlat = T)
      d       <- as.matrix(d)
      diag(d) <- NA
      closest <- lapply(1:ncol(d), function(z) {
        index <- which.min(d[, z])
        dista <- d[index, z]
        return(data.frame(ClosestTrack = index, Distance = dista))
      }) %>% do.call(rbind, .)
      return(cbind(x, closest))
  })) %>%
  unnest(Data) %>%
  arrange(TrackID, Timestamp)

# # For each step, identify the closest step from another group
# testing <- tracks %>%
#   do.call(rbind, .) %>%
#   nest(Data = -Timestamp) %>%
#   mutate(Data = map(Data, function(x) {
#       x <- testing$Data[[1]]
#       steps <- lapply(1:nrow(x), function(z) {
#         vect(rbind(c(x$x[z], x$y[z]), c(x$x_to[z], x$y_to[z])), type = "lines", crs = "epsg:4326")
#       }) %>% do.call(rbind, .)
#       distance(steps)
#
#       d       <- distance(cbind(x$x, x$y), lonlat = T)
#       d       <- as.matrix(d)
#       diag(d) <- NA
#       closest <- lapply(1:ncol(d), function(z) {
#         index <- which.min(d[, z])
#         dista <- d[index, z]
#         return(data.frame(ClosestTrack = index, Distance = dista))
#       }) %>% do.call(rbind, .)
#       return(cbind(x, closest))
#   })) %>%
#   unnest(Data) %>%
#   arrange(TrackID, Timestamp)


# Put tracks together
tracks <- do.call(rbind, tracks)

# Convert to tibble (this saves an amazing amount of space)
tracks <- as_tibble(tracks)

# Store the simulations
write_rds(tracks, design$Filename[i])

################################################################################
#### Combine Simulations
################################################################################
# Once the simulations are done, let's combine them into a single big dataframe
sims <- lapply(1:nrow(design), function(x) {
  dat <- read_rds(design$Filename[x])
  dat$SimID      <- x
  dat$FloodLevel <- design$FloodLevel[x]
  dat$SourceArea <- design$SourceArea[x]
  return(dat)
}) %>% do.call(rbind, .)

# Because in each simulation we start off with new IDs for trajectories, they
# are not unique across simulations. We thus combine ID and SimID to create an
# ID that is unqiue to each simulated path, across all simulations
sims <- sims %>%
  group_by(SimID, TrackID) %>%
  mutate(TrackID = cur_group_id())

# Make sure it worked
table(table(sims$TrackID))

# Collect garbage
gc()

# Let's also create a step counter, indicating the number of the step in its
# respective trajectory
sims <- sims %>%
  group_by(TrackID) %>%
  mutate(StepNumber = row_number())

# Check object size
format(object.size(sims), units = "Gb")

# Ungroup
sims <- ungroup(sims)

# Write to an rds
write_rds(sims, "03_Data/03_Results/DispersalSimulation.rds")

# Remove separate files to free some space
# file.remove(design$Filename)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/04_Simulation.rds")

# Print to terminal
cat("Done :)\n")
