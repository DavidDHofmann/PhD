################################################################################
#### Distance to Humans
################################################################################
# Description: Computing the distance to humans for each simulated step

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(terra)     # To handle spatial data
library(raster)    # To handle spatial data
library(tidyverse) # To wrangle data
library(maptools)  # To compute point density
library(spatstat)  # To compute point density

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load distance to humans layer
humans <- rast("03_Data/02_CleanData/DistanceToHumans.tif")

# Load dispersal simulations
sims <- read_rds("03_Data/03_Results/DispersalSimulation.rds")

# Keep only desired columns
sims <- sims[, c("x", "y", "sl_", "ta_", "TrackID", "StepNumber", "SourceArea", "FloodLevel")]

################################################################################
#### Distance to Humans
################################################################################
# Create directories into which we will store the output rasters
dir.create("03_Data/03_Results/99_Distances", showWarnings = F)

# Calculate the distance to humans for each step
sims$DistanceToHumans <- terra::extract(humans, cbind(sims$x, sims$y))[, 1]

# Design through which to loop
design <- expand_grid(
    FloodLevel = unique(sims$FloodLevel)
  , SourceArea = unique(sims$SourceArea)
)
design$Filename <- with(design, paste0("03_Data/03_Results/99_Distances/Distance_SourceArea", SourceArea, "_FloodLevel", FloodLevel, ".tif"))

# Loop throgh the design and generate density maps of the locations that are
# within a threshold distance to the nearest "human-influenced" pixel
if (!file.exists("03_Data/03_Results/Distance.tif")) {
    distmaps <- list()
    for (i in 1:nrow(design)) {
      cat("Computing distance map", i, "out of", nrow(design), "\n")
      if (file.exists(design$Filename[i])) {
        distmaps[[i]] <- rast(design$Filename[[i]])
      } else {
        closeto <- subset(sims
          , DistanceToHumans < 500 &
            FloodLevel == design$FloodLevel[i] &
            SourceArea == design$SourceArea[i]
        )
        closeto <- vect(x = cbind(closeto$x, closeto$y)
          , crs  = "epsg:4326"
          , atts = closeto
        )
        distmaps[[i]] <- terra::rasterize(closeto
          , humans
          , fun        = function(i) { length(i) }
          , background = 0
        )
        writeRaster(distmaps[[i]], design$Filename[i], overwrite = T)
      }
    }

    # Combine maps, reproject them, crop, and store them
    combined <- do.call(c, combined)
    writeRaster(combined, "03_Data/03_Results/Distance.tif", overwrite = T)
  } else {
    combined <- stack("03_Data/03_Results/Distance.tif")
}

# Add maps to the tibble
design <- mutate(design, Distance = lapply(1:nlyr(combined), function(x) {
  raster(combined[[x]])
}))

# Store the tibble
write_rds(design, "03_Data/03_Results/Distance.rds")
