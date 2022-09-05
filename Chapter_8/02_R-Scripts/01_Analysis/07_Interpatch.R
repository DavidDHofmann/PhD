################################################################################
#### Interpatch Connectivity
################################################################################
# Description: Computing interpatch connectivity from simulated dispersal paths

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(pbmcapply)      # To run stuff in parallel
library(igraph)         # For network analysis
library(rgeos)          # To manipulate spatial objects
library(ggnetwork)      # To plot network using ggplot
library(sf)             # To plot spatial features
library(ggspatial)      # To add scale bars etc to plots

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load reference raster
r <- rast("03_Data/02_CleanData/ReferenceRaster.tif")

# Load source areas
area <- vect("03_Data/02_CleanData/SourceAreas.shp")

# Rasterize them to the reference raster
area_r <- terra::rasterize(area, y = r, field = "ID")

# # Visualize them
# plot(area_r, main = "Source Areas")
# plot(area, add = T)
# text(area, "ID", cex = 0.5, halo = T)

################################################################################
#### Prepare Simulations
################################################################################
# Load dispersal simulations
sims <- read_rds("03_Data/03_Results/DispersalSimulation.rds")

# Keep only desired columns
sims <- sims[, c("x", "y", "TrackID", "StepNumber", "SourceArea", "FloodLevel")]

# Count number of simulations per source area
n_sample <- sims %>%
  select(TrackID, SourceArea, FloodLevel) %>%
  distinct() %>%
  count(SourceArea, FloodLevel) %>%
  pull(n) %>%
  unique()

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# Identify through which national parks the dispersers moved
visits <- data.frame(
    TrackID    = sims$TrackID
  , StepNumber = sims$StepNumber
  , FloodLevel = sims$FloodLevel
  , SourceArea = sims$SourceArea
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , Area       = raster::extract(raster(area_r), sims)
)

# Ignore any step in an "na" or "nan" area
visits <- subset(visits, !is.na(Area) & !is.nan(Area))

# Calculate for each step the distance to the first coordinate. We'll use this
# to determine how far an individual had to disperse before reaching another
# area
visits <- visits %>%
  nest(data = -TrackID) %>%
  mutate(data = pbmclapply(data
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x) {

      # Project coordinates
      coords <- reprojCoords(
          xy   = x[, c("x", "y")]
        , from = CRS("+init=epsg:4326")
        , to   = CRS("+init=epsg:32734")
      )

      # Compute distance to first coordinate
      first <- coords[1, ]
      distance <- sqrt((coords[, 1] - first[1]) ** 2 + (coords[, 2] - first[2]) ** 2)
      x$DistanceFromFirst <- distance

      # Return the resulting object
      return(x)
  })) %>%
  unnest(data)

# Identify how long it takes to reach the different areas
visits <- visits %>%
  rename(SourceArea = SourceArea, CurrentArea = Area) %>%
  group_by(TrackID, FloodLevel, SourceArea, CurrentArea) %>%
  summarize(
      StepNumber        = min(StepNumber)
    , DistanceFromFirst = min(DistanceFromFirst)
    , .groups           = "drop"
  ) %>%

  arrange(TrackID, StepNumber)

# Compute summary statistics by source area
summarizeVisits <- function(visits) {
  visits %>%
    group_by(FloodLevel, SourceArea, CurrentArea) %>%
    summarize(
        MeanStepNumber = mean(StepNumber)
      , SDStepNumber   = sd(StepNumber)
      , Frequency      = n()
      , .groups        = "drop"
    )
}

# Try it
summarizeVisits(visits)

# Generate bootstrap samples to compute standard deviation for the number of
# successful dispersals between source areas and the number of steps required
visits_nested <- nest(visits, Data = -c(TrackID, FloodLevel, SourceArea))
bootstrapped <- pbmclapply(1:1000, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  boot <- visits_nested %>%
    group_by(FloodLevel, SourceArea) %>%
    slice_sample(replace = T, n = n_sample) %>%
    unnest(Data) %>%
    summarizeVisits()
  boot$Bootstrap <- x
  return(boot)
})
bootstrapped <- do.call(rbind, bootstrapped)

# Bind and compute confidence intervals
visits_bootstrapped <- bootstrapped %>%
  group_by(FloodLevel, SourceArea, CurrentArea) %>%
  summarize(
      StepNumber   = mean(MeanStepNumber)
    , StepNumberSE = sd(MeanStepNumber)   # Bootstrapped sd is the metric's se
    , Freq         = mean(Frequency)
    , FreqSE       = sd(Frequency)        # Bootstrapped sd is the metric's se
    , .groups      = "drop"
  )

# Store visits to file
write_rds(visits_bootstrapped, "03_Data/03_Results/BootstrappedInterpatchConnectivity.rds")
