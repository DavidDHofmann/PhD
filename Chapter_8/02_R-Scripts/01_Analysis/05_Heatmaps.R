################################################################################
#### Heatmaps
################################################################################
# Description: Computing heatmaps from simulated dispersal paths

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
library(spatstat)       # To rasterize lines quickly
library(maptools)       # To rasterize lines quickly

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load dispersal simulations
sims <- read_rds("03_Data/03_Results/DispersalSimulation.rds")

# Keep only columns
sims <- sims[, c("x", "y", "TrackID", "StepNumber", "Area", "FloodLevel")]

# # Subsample
# sims <- subset(sims, FloodLevel == "Max")
# sims <- subset(sims, TrackID %in% sample(unique(sims$TrackID), size = 1000))

# Reproject coordinates to utm (required for spatstat)
sims[, c("x", "y")] <- reprojCoords(
    xy   = sims[, c("x", "y")]
  , from = CRS("+init=epsg:4326")
  , to   = CRS("+init=epsg:32734")
)

# Load the reference raster
r <- rast("03_Data/02_CleanData/ReferenceRaster.tif")

# Prepare extent that encompassess all coordinates + some buffer
ext <- extent(min(sims$x), max(sims$x), min(sims$y), max(sims$y)) +
  c(-1000, +1000, -1000, +1000)

# Span a raster with desired resolution
r <- raster(ext, res = 1000)
values(r) <- runif(ncell(r))
crs(r) <- CRS("+init=epsg:32734")

# Collect garbage
gc()

# Check out the number of rows
nrow(sims) / 1e6

# Create a dataframe with all source points and points in time at which we want
# to rasterize trajectories
rasterized <- expand_grid(
    Steps = c(68, 125, 250, 500, 1000, 2000)
  , Area  = unique(sims$Area)
  , Flood = unique(sims$FloodLevel)
)

# Add a column for temporary but unique filename. Make sure the tempdir has
# plenty of storage.
rasterized$filename <- tempfile(
    pattern = paste0(
        "Steps_", rasterized$Steps
      , "_Area_", rasterized$Area
      , "_Flood_", rasterized$Flood
      , "_"
    )
  , fileext = ".tif"
)

# Loop through the study design and reasterize trajectories
heatmaps <- list()
pb <- txtProgressBar(min = 0, max = nrow(rasterized), style = 3)
for (i in 1:nrow(rasterized)) {

  # Create heatmap
  heatmaps[[i]] <- rasterizeSims(
      simulations = sims
    , raster      = r
    , steps       = rasterized$Steps[i]
    , area        = rasterized$Area[i]
    , flood       = rasterized$Flood[i]
    , messages    = F
    , mc.cores    = detectCores() - 1
  )

  # Clean garbage
  gc()

  # Print update
  setTxtProgressBar(pb, i)

}

# Combine maps
combined <- stack(heatmaps)

# Reproject them
combined <- rast(combined)
combined <- terra::project(combined, CRS("+init=epsg:4326"), method = "bilinear")
combined <- stack(combined)

# Crop them to our reference raster
r <- raster("03_Data/02_CleanData/ReferenceRaster.tif")
combined <- crop(combined, r)

# Store to file
writeRaster(combined, "03_Data/03_Results/Heatmaps.tif", overwrite = T)

# Add maps to the tibble
rasterized <- mutate(rasterized, heatmap = lapply(1:nlayers(combined), function(x){
  combined[[x]]
}))

# Store to file
write_rds(rasterized, "03_Data/03_Results/Heatmaps.rds")
