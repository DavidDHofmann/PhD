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
library(sf)             # For plotting

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
hrs <- buffer(hrs, width = 0)
# mapview(as(hrs, "Spatial"))

# Function to sample locations from the different HRs
sampleLocations <- function(p, ID = NULL) {
  # p <- hrs
  # ID <- sample(hrs$ID, size = 10, replace = F)
  if (is.null(ID)) {
    p_sub <- p
  } else {
    p_sub <- p[p$ID %in% ID, ]
  }
  locs <- lapply(1:length(p_sub), function(x) {
    pts <- spatSample(p_sub[x, ], size = 20)
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

# Generate updated locations for residents
residents <- sampleLocations(hrs, ID = sample(hrs$ID, 50))
values(residents) <- data.frame(PackID = residents$ID)
hrs$Occupied <- hrs$ID %in% residents$PackID

# From some resident groups, we may see dispersers
dispersers     <- residents[sample(1:nrow(residents), size = 20, replace = F), ]
dispersers$Sex <- sample(c("M", "F"), size = length(dispersers), replace = T)
dispersers$Sex <- factor(dispersers$Sex, levels = c("M", "F"))
dispersers$CoalitionID <- 1:length(dispersers)

# Convert all to dataframes
dispersers <- as.data.frame(dispersers, geom = "XY")
residents  <- as.data.frame(residents, geom = "XY")

# Plot them
ggplot() +
  geom_sf(data = st_as_sf(hrs), aes(fill = Occupied), col = "white") +
  geom_point(data = residents, aes(x = x, y = y), col = "black", size = 3) +
  geom_point(data = dispersers, aes(x = x, y = y, col = Sex)) +
  scale_color_manual(values = c("cornflowerblue", "orange")) +
  scale_fill_manual(values = c("gray90", "gray70")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Run the simulation for each source point
tracks <- pbmclapply(
    X                   = 1:nrow(dispersers)
  , mc.cores            = detectCores() - 1
  , ignore.interactive  = T
  , FUN                 = function(x) {
    sim <- suppressWarnings(
      disperse(
          source    = cbind(dispersers$x[x], dispersers$y[x])
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

    # Assign some more information (We will need to make sure that each group
    # has a unique ID!!!)
    sim$PackID      <- dispersers$PackID[x]
    sim$CoalitionID <- dispersers$CoalitionID[x]
    sim$Sex         <- dispersers$Sex[x]

    # Create new columns containing the endpoints of each step
    sim$x_to <- lead(sim$x)
    sim$y_to <- lead(sim$y)
    sim <- dplyr::select(sim, x, y, x_to, y_to, everything())

    # Return the simulation
    return(sim)
})

# Bind all together
tracks <- tracks %>%
  do.call(rbind, .) %>%
  subset(!is.na(x_to) & !is.na(y_to))

# Identify closest group of (1) different sex dispersers from a different pack,
# or (2) different pack residents
tracks <- tracks %>%
  nest(Data = -Timestamp) %>%
  mutate(Data = map(Data, function(x) {

      # Get relevant coordiantes
      dis <- select(x, x = x_to, y = y_to, CoalitionID, PackID, Sex)
      res <- select(residents, x, y, PackID)

      # Find distances between dispersers and dispersers, as well as between
      # dispersers and residents
      d1 <- distance(cbind(dis$x, dis$y), cbind(dis$x, dis$y), lonlat = T)
      d2 <- distance(cbind(dis$x, dis$y), cbind(res$x, res$y), lonlat = T)

      # Filter to find closest other-sex, other-pack dispersal coalition
      filter1 <- outer(dis$CoalitionID, dis$CoalitionID, FUN = function(x, y) {x != y})
      filter2 <- outer(dis$Sex, dis$Sex, FUN = function(x, y) {x != y})

      # Filter to find closest other-pack resident coalition
      filter3 <- outer(dis$PackID, res$PackID, FUN = function(x, y) {x != y})

      # Apply filters
      d1[!filter1 | !filter2] <- NA
      d2[!filter3] <- NA

      # Extract IDs
      closest <- lapply(1:nrow(d1), function(z) {

        # Closest dispersers
        index_dispersers <- which.min(d1[z, ])
        if (length(index_dispersers) == 0) {
          close_dispersers <- NA
          dista_dispersers <- NA
        } else {
          close_dispersers <- dis$CoalitionID[index_dispersers]
          dista_dispersers <- d1[z, index_dispersers]
        }

        # Closest residents
        index_residents <- which.min(d2[z, ])
        if (length(index_residents) == 0) {
          close_residents <- NA
          dista_residents <- NA
        } else {
          close_residents <- res$PackID[index_residents]
          dista_residents <- d2[z, index_residents]
        }

        # Put all into a dataframe
        result <- data.frame(
            ClosestDispersers    = close_dispersers
          , DistanceToDispersers = dista_dispersers
          , ClosestResidents     = close_residents
          , DistanceToResidents  = dista_residents
        )
        return(result)
      }) %>% do.call(rbind, .)
      return(cbind(x, closest))
  })) %>%
  unnest(Data) %>%
  arrange(PackID, Timestamp)

# Plot them
ggplot() +
  geom_sf(data = st_as_sf(hrs), aes(fill = Occupied), col = "white") +
  geom_path(data = tracks, aes(x = x, y = y, col = Sex, group = CoalitionID)) +
  geom_point(data = residents, aes(x = x, y = y), col = "black") +
  scale_color_manual(values = c("cornflowerblue", "orange")) +
  scale_fill_manual(values = c("gray90", "gray70")) +
  theme_minimal() +
  theme(legend.position = "bottom")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/03_Simulation.rds")

# Print to terminal
cat("Done :)\n")
