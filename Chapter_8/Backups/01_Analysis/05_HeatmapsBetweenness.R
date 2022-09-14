################################################################################
#### Heatmaps & Betweenness Maps
################################################################################
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
library(igraph)         # For network analysis

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load reference raster
r <- raster("03_Data/02_CleanData/ReferenceRaster.tif")

# Load dispersal simulations
sims <- read_rds("03_Data/03_Results/DispersalSimulation.rds")

# Keep only columns
sims <- sims[, c("x", "y", "TrackID", "StepNumber", "SourceArea", "FloodLevel")]

# Reproject coordinates to utm (required for spatstat)
sims[, c("x", "y")] <- reprojCoords(
    xy   = sims[, c("x", "y")]
  , from = CRS("+init=epsg:4326")
  , to   = CRS("+init=epsg:32734")
)

# Generate a raster on which we will compute the heatmaps
ext <- extent(min(sims$x), max(sims$x), min(sims$y), max(sims$y)) +
  c(-1000, +1000, -1000, +1000)

# Span rasters with desired resolutions
r_heat <- raster(ext, res = 1000)
r_betw <- raster(ext, res = 2500)
values(r_heat) <- 1:ncell(r_heat)
values(r_betw) <- 1:ncell(r_betw)
crs(r_heat) <- crs(r_betw) <- CRS("+init=epsg:32734")

################################################################################
#### Local Metrics
################################################################################
# Folder into which we store the separate maps
dir.create("03_Data/03_Results/99_Heatmaps", showWarnings = F)
dir.create("03_Data/03_Results/99_Betweenness", showWarnings = F)

# Create a dataframe with all source points and number of steps for which we
# want to prepare said metrics
design <- expand_grid(
    Steps      = c(500, 1000, 2000)
  , SourceArea = unique(sims$SourceArea)
  , FloodLevel = unique(sims$FloodLevel)
)
design$FilenameHeatmap <- with(design, paste0("03_Data/03_Results/99_Heatmaps/Heatmap_Steps", Steps, "_SourceArea", SourceArea, "_FloodLevel", FloodLevel, ".tif"))
design$FilenameBetweenness <- with(design, paste0("03_Data/03_Results/99_Betweenness/Betweenness_Steps", Steps, "_SourceArea", SourceArea, "_FloodLevel", FloodLevel, ".tif"))

# Shuffle the design to get more reliable estimates of computation times
design <- design[sample(1:nrow(design), size = nrow(design), replace = F), ]

# Loop through the design and generate heatmaps
if (!file.exists("03_Data/03_Results/HeatmapsLocal.tif")) {
    heatmaps <- list()
    for (i in 1:nrow(design)) {
      cat("Computing heatmap", i, "out of", nrow(design), "\n")
      if (file.exists(design$FilenameHeatmap[i])) {
        heatmaps[[i]] <- raster(design$FilenameHeatmap[[i]])
      } else {
        heatmaps[[i]] <- rasterizeSims(
            simulations = sims
          , raster      = r_heat
          , steps       = design$Steps[i]
          , area        = design$SourceArea[i]
          , flood       = design$FloodLevel[i]
          , messages    = F
          , mc.cores    = detectCores() - 1
          , filename    = design$FilenameHeatmap[i]
        )
      }
    }

    # Combine maps, reproject them, crop, and store them
    combined <- stack(heatmaps)
    combined <- rast(combined)
    combined <- terra::project(combined, CRS("+init=epsg:4326"), method = "bilinear")
    combined <- stack(combined)
    combined <- crop(combined, r)
    writeRaster(combined, "03_Data/03_Results/HeatmapsLocal.tif", overwrite = T)
  } else {
    combined <- stack("03_Data/03_Results/HeatmapsLocal.tif")
}

# Add maps to the tibble
design <- mutate(design, Heatmap = lapply(1:nlayers(combined), function(x) {
  combined[[x]]
}))

# Loop through the design and generate betweenness maps
if (!file.exists("03_Data/03_Results/BetweennessLocal.tif")) {
    betweenness <- pbmclapply(1:nrow(design), ignore.interactive = T, mc.cores = detectCores() / 2, function(i) {
      if (file.exists(design$FilenameBetweenness[i])) {
        raster(design$FilenameBetweenness[[i]])
      } else {
        betweenSims(
            simulations = sims
          , raster      = r_betw
          , steps       = design$Steps[i]
          , area        = design$SourceArea[i]
          , flood       = design$FloodLevel[i]
          , messages    = T
          , mc.cores    = 1
          , filename    = design$FilenameBetweenness[i]
        )
      }
    })

    # Combine maps, reproject them, crop, and store them
    combined <- stack(betweenness)
    combined <- rast(combined)
    combined <- terra::project(combined, CRS("+init=epsg:4326"), method = "bilinear")
    combined <- stack(combined)
    combined <- crop(combined, r)
    writeRaster(combined, "03_Data/03_Results/BetweennessLocal.tif", overwrite = T)
  } else {
    combined <- stack("03_Data/03_Results/BetweennessLocal.tif")
}

# Add maps to the tibble
design <- mutate(design, Betweenness = lapply(1:nlayers(combined), function(x) {
  combined[[x]]
}))

# Store to file
write_rds(design, "03_Data/03_Results/HeatmapsBetweennessLocal.rds")

################################################################################
#### Global Metrics
################################################################################
# Create a dataframe with all source points and number of steps for which we
# want to prepare said metrics
design <- expand_grid(
    Steps               = c(500, 1000, 2000)
  , FloodLevel          = unique(sims$FloodLevel)
)
design$FilenameHeatmap <- with(design, paste0("03_Data/03_Results/99_Heatmaps/Heatmap_Steps", Steps, "_FloodLevel", FloodLevel, ".tif"))
design$FilenameBetweenness <- with(design, paste0("03_Data/03_Results/99_Betweenness/Betweenness_Steps", Steps, "_FloodLevel", FloodLevel, ".tif"))

# Shuffle the design to get more reliable estimates of computation times
design <- design[sample(1:nrow(design), size = nrow(design), replace = F), ]

# Loop through the design and generate heatmaps
if (!file.exists("03_Data/03_Results/HeatmapsGlobal.tif")) {
    heatmaps <- list()
    for (i in 1:nrow(design)) {
      cat("Computing heatmap", i, "out of", nrow(design), "\n")
      if (file.exists(design$FilenameHeatmap[i])) {
        heatmaps[[i]] <- raster(design$FilenameHeatmap[[i]])
      } else {
        heatmaps[[i]] <- rasterizeSims(
            simulations = sims
          , raster      = r_heat
          , steps       = design$Steps[i]
          , flood       = design$FloodLevel[i]
          , messages    = F
          , mc.cores    = detectCores() - 1
          , filename    = design$FilenameHeatmap[i]
        )
      }
    }

    # Combine maps, reproject them, crop, and store them
    combined <- stack(heatmaps)
    combined <- rast(combined)
    combined <- terra::project(combined, CRS("+init=epsg:4326"), method = "bilinear")
    combined <- stack(combined)
    combined <- crop(combined, r)
    writeRaster(combined, "03_Data/03_Results/HeatmapsGlobal.tif", overwrite = T)
  } else {
    combined <- stack("03_Data/03_Results/HeatmapsGlobal.tif")
}

# Add maps to the tibble
design <- mutate(design, Heatmap = lapply(1:nlayers(combined), function(x) {
  combined[[x]]
}))

# Loop through the design and generate betweenness maps
if (!file.exists("03_Data/03_Results/BetweennessGlobal.tif")) {
    betweenness <- pbmclapply(1:nrow(design), ignore.interactive = T, mc.cores = detectCores() / 4, function(i) {
      if (file.exists(design$FilenameBetweenness[i])) {
        raster(design$FilenameBetweenness[[i]])
      } else {
        betweenSims(
            simulations = sims
          , raster      = r_betw
          , steps       = design$Steps[i]
          , flood       = design$FloodLevel[i]
          , messages    = T
          , mc.cores    = 1
          , filename    = design$FilenameBetweenness[i]
        )
      }
    })

    # Combine maps, reproject them, crop, and store them
    combined <- stack(betweenness)
    combined <- rast(combined)
    combined <- terra::project(combined, CRS("+init=epsg:4326"), method = "bilinear")
    combined <- stack(combined)
    combined <- crop(combined, r)
    writeRaster(combined, "03_Data/03_Results/BetweennessGlobal.tif", overwrite = T)
  } else {
    combined <- stack("03_Data/03_Results/BetweennessGlobal.tif")
}

# Add maps to the tibble
design <- mutate(design, Betweenness = lapply(1:nlayers(combined), function(x) {
  combined[[x]]
}))

# Store to file
write_rds(design, "03_Data/03_Results/HeatmapsBetweennessGlobal.rds")
