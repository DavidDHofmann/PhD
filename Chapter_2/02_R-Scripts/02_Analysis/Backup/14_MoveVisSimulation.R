############################################################
#### Create Movement Animation of Simulation
############################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load packages
library(raster)       # For general raster manipulation
library(tidyverse)    # For data wrangling
library(moveVis)      # To animate movement
library(move)         # To coerce data to move objects
library(lubridate)    # To handle dates
library(viridis)      # For nice colors
library(tmap)         # To get a nice basemap
library(tmaptools)    # To get access to osm data
library(tictoc)       # To keep track of time

# Make use of multicore
beginCluster()

# Also tell moveVis to use multicore
use_multicore()

# Define the extent that we want to plot
extent <- "03_Data/02_CleanData/00_General_Raster250.tif" %>%
  raster() %>%
  extent()

# Load simulated trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")

# Coerce the data to a move object
move <- move(
    x       = sims$x
  , y       = sims$y
  , time    = sims$Timestamp
  , data    = sims
  , animal  = as.character(sims$ID)
  , proj    = CRS("+init=epsg:4326")
)

# We want white trajectories
move$colour <- "white"

# Prepare frames
tic()
frames <- frames_spatial(move
  , tail_colour   = "white"
  , trace_colour  = "white"
  , trace_show    = TRUE
  , path_legend   = FALSE
  , path_size     = 0.50
  , tail_size     = 0.25
  , map_service   = "mapbox"
  , map_type      = "satellite"
  , map_token     = "pk.eyJ1IjoiZG9keDkiLCJhIjoiY2p3dnltejJjMGR4YjN5bXp0ZjA2ZXBzMCJ9.4hirgQ-1SfJ2KHI7SR54cQ"
  , ext           = extent(move)
  ) %>%
  add_progress() %>%
  add_northarrow(colour = "white", size = 3) %>%
  add_scalebar(colour = "white", distance = 50) %>%
  add_labels(x = "Longitude", y = "Latitude")
toc()

# Look at a desired frame to make sure everything looks as desired.
frames[[50]]

# Store the animation
animate_frames(
    frames
  , width     = 1920
  , height    = 1080
  , out_file  = "04_Manuscript/SimulatedDispersal.mov"
  , overwrite = TRUE
)

# Terminate Cluster
endCluster()
