################################################################################
#### Betweenness
################################################################################
# Description: Computing betweenness from simulated dispersal paths

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

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load dispersal simulations
sims <- read_rds("03_Data/03_Results/DispersalSimulation.rds")
sims <- subset(sims, FloodLevel != "Mean")
# sims <- subset(sims, Area != 2)

# Keep only certain columns
sims <- sims[, c("x", "y", "TrackID", "StepNumber", "SourceArea", "FloodLevel")]

# Function to retrieve the visitation history from a sequence of values
visitHist <- function(x, singlecount = F) {
  transitions <- data.frame(from = lag(x), to = x) %>%
    group_by(from, to) %>%
    na.omit() %>%
    summarize(TotalConnections = n(), .groups = "drop")
  if (singlecount){
    transitions$TotalConnections <- 1
  }
  return(transitions)
}

# Function to compute betweenness
betweennessMap <- function(network = NULL, raster = NULL, tempfile = F) {
    betweenness <- raster
    values(betweenness) <- betweenness(network)
    names(betweenness) <- "betweenness"
    betweenness <- writeRaster(betweenness, tempfile())
    return(betweenness)
}

# Load the reference raster
r <- raster("03_Data/02_CleanData/ReferenceRaster.tif")

# Need to coarsen resolution
r  <- aggregate(r, fact = 2500 / 250, fun = max)

# Note that the resolution of the raster determines the number of vertices in
# the resulting network (each raster cell will be turned into a vertex). Let's
# therefore keep track of the IDs of each raster cell.
vertices <- 1:ncell(r)

# The resolution of the raster therefore also determines the layout of the
# network (i.e. the spatial location of each vertex). Let's keep track of the
# coordinates of each raster cell (i.e. each vertex) so that we can nicely plot
# the network graph afterwards.
lay  <- as.matrix(as.data.frame(r, xy = T)[, c(1, 2)])

# Fill the rasters with unique cell values (we will use the values as cell IDs)
values(r)  <- vertices

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# At each coordinate of the simulated trajectories we now extract the cell IDs
# from the different rasters
visits <- data.frame(
    TrackID    = sims$TrackID
  , StepNumber = sims$StepNumber
  , FloodLevel = sims$FloodLevel
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , r          = raster::extract(r, sims)
)

# # Remove simulations
# rm(sims)
# gc()

################################################################################
#### Individual Betweenness Maps
################################################################################
# Prepare design through which we want to loop
design <- expand_grid(
    Steps      = c(500, 1000, 2000)
  , SourceArea = unique(sims$SourceArea)
  , FloodLevel = unique(sims$FloodLevel)
)

################################################################################
#### Combined Betweenness Maps
################################################################################
# Prepare design through which we want to loop
design <- expand_grid(
    Steps       = c(500, 1000, 2000)
  , FloodLevel  = unique(sims$FloodLevel)
)

# Add a column for temporary but unique filename. Make sure the tempdir has
# plenty of storage.
design$filename <- tempfile(
    pattern = paste0(
        "Steps", design$Steps
      , "_FloodLevel", design$FloodLevel
      , "_"
    )
  , fileext = ".tif"
)

# Loop through the design and calculate betweenness map
maps <- list()
pb <- txtProgressBar(min = 0, max = nrow(design), style = 3)
for (i in 1:nrow(design)) {

  # Subset data to desired steps
  sub <- visits[
    visits$StepNumber <= design$Steps[i] &
    visits$FloodLevel == design$FloodLevel[i], ]

  # Nest tracks
  sub <- nest(sub, data = -TrackID)

  # Create visitation history
  cat("Getting visitation history...\n")
  history <- pbmclapply(
      X                  = sub$data
    , ignore.interactive = T
    , mc.cores           = 1
    , FUN                = function(y) {
        visitHist(y$r, singlecount = T)
    }) %>%
    do.call(rbind, .) %>%
    group_by(from, to) %>%
    summarize(TotalConnections = sum(TotalConnections), .groups = "drop") %>%
    ungroup() %>%
    mutate(weight = mean(TotalConnections) / TotalConnections)

  # Create network
  cat("Creating graph...\n")
  net <- graph_from_data_frame(history, vertices = vertices)

  # Calculate Betweenness
  cat("Calculating betweenness...\n")
  maps[[i]] <- betweennessMap(
      network  = net
    , raster   = r
  )

  # Print update
  setTxtProgressBar(pb, i)

}

# Combine maps
combined <- stack(maps)

# Store to file
writeRaster(combined, "03_Data/03_Results/BetweennessMaps.tif", overwrite = T)

# Add maps to the tibble
design <- mutate(design, betweenness = lapply(1:nlayers(combined), function(x) {
  combined[[x]]
}))
print(design)

# Store to file
write_rds(design, "03_Data/03_Results/Betweenness.rds")
