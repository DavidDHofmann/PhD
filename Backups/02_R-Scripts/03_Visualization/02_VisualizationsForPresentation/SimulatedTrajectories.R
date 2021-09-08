################################################################################
#### Plot of Simulated Trajectories
################################################################################
# Description: A simple plot of the Dispersal Durations

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)       # For spatial data
library(rgeos)        # For spatial data
library(davidoff)     # Custom functions
library(tidyverse)    # Data wrangling
library(parallel)     # Parallel computing

################################################################################
#### Simulated Trajectories
################################################################################
# Old directory
wd <- "/media/david/My Passport/Backups/WildDogs/15. PhD/00_WildDogs"
setwd(wd)

# Load required data
source_areas  <- shapefile("03_Data/03_Results/99_SourceAreas2.shp")
source_points <- shapefile("03_Data/03_Results/99_SourcePoints2.shp")

# Load some trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Only keep those from the static source point
sims <- subset(sims, PointSampling == "Static")

# Create animation
ani.options(interval = .1, ani.width = 1980, ani.height = 1980)
saveVideo({
  for (i in 1:2000){
    trajs <- sims2tracks(sims, steps = i)
    plot(source_areas, col = "gray50", border = "gray20", bg = "black")
    plot(trajs, add = T, col = "orange", lwd = 2.5)
  }
}, video.name = "99_Simulations99.mp4")
