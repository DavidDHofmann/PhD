################################################################################
#### Simulating Trajectories using the fitted Parameters
################################################################################
# Description: Simulating movement using the parameters that we estimated in the
# previous script
# Authors: David, Eva, & Dilsad
# Date: November 2020

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)       # To handle spatial data
library(rgeos)        # To manipulate spatial data
library(tidyverse)    # For data wrangling
library(pbmcapply)    # To run tasks in parallel
library(rgdal)        # To store spatial data

# Set the working directory
setwd("/home/david/ownCloud/Dokumente/Bibliothek/Wissen/R-Scripts/ProgrammingProject")

# Load custom functions
source("00_Functions.R")

# Load estimated parameters
dat <- read_rds("Data/Output/FittedDistributions.rds")

# Set seed for reproducability
set.seed(1234)

################################################################################
#### Simulate Movement
################################################################################
# We want to simulate movement in switzerland, so let's get the boundaries of
# switzerland for reference of size. We'll convert the data to utm so that we
# can work in meters
che <- getData("GADM", country = "CHE", level = 0)
che <- spTransform(che, CRS("+proj=utm +zone=32 ellps=WGS84"))

# Let's also draw a rectangular bounding box around switzerland. We'll use it as
# movement boundary.
area_move <- as(extent(che) + c(-50, +50, -50, +50) * 1000, "SpatialPolygons")
area_move$Description <- "MovementArea"
crs(area_move) <- crs(che)

# Let's assume that we want to set up our camera traps to cover an area of 10 x
# 10 kilometers. To select a location, let's randomly place the center point of
# the camera grid somewhere in switzerland.
center <- spsample(che, n = 1, type = "random")

# Now create the area for the camera trap around that point
area_cams <- as(extent(center) + c(-5, +5, -5, +5) * 1000, "SpatialPolygons")
area_cams$Description <- "CameraArea"
crs(area_cams) <- crs(che)

# Visualize the different regions
plot(area_move, border = "gray50", lty = 2)
plot(area_cams, add = T, border = "red")
plot(che, add = T)

# Simulate trajectories for each species. Note that all individuals are
# released from inside the "area_cams" polygon
dat$Simulations <- lapply(1:nrow(dat), function(x){

  # Simulate trajectories
  paths <- pbmclapply(1:2000
    , mc.cores           = detectCores() - 1
    , ignore.interactive = T
    , FUN                = function(y){
      move(
          n       = 365 * 24 / 4
        , start   = coordinates(spsample(area_cams, n = 1, type = "random"))
        , ext     = area_move
        , dist_sl = dat$dist_sl[[x]]
        , dist_ta = dat$dist_ta[[x]]
      )
  }) %>% do.call(rbind, .)

  # Assign ID and Species
  paths$ID <- 1:length(paths)
  paths$Species <- dat$Species[x]

  # Return the simulations
  return(paths)
})

# Put all simulations into a single object
paths <- do.call(rbind, dat$Simulations)

# Create unique IDs
paths$ID <- 1:nrow(paths)

# Assign correct crs
crs(paths) <- crs(area_cams)

# Visualize (some) simulations and highlight one of the paths
plot(area_move, border = "gray50", lty = 2)
plot(che, add = T)
plot(paths[1:100, ], add = T, col = "gray80", lwd = 0.1)
plot(paths[1, ], add = T, col = "red")
plot(area_cams, add = T, border = "purple")
legend("bottomright"
  , legend = c("swiss border", "movement border", "simulations", "camera area")
  , col    = c("black", "gray50", "gray80", "purple")
  , lty    = c(1, 2, 1, 1)
)

# Store the data to files
writeOGR(paths, "Data/Output/", "Simulations", "ESRI Shapefile", overwrite = T)
writeOGR(area_move, "Data/Output/", "MovementArea", "ESRI Shapefile", overwrite = T)
writeOGR(area_cams, "Data/Output/", "CameraArea", "ESRI Shapefile", overwrite = T)
writeOGR(che, "Data/Output/", "Switzerland", "ESRI Shapefile", overwrite = T)
