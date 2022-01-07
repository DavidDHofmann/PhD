############################################################
#### Comparing Accuracy of Simulated Paths and Least-Cost Paths
############################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load packages
library(tidyverse)    # For data wrangling
library(terra)        # For raster manipulation
library(raster)       # For raster manipulation
library(viridis)      # For nicer colors
library(rgeos)        # To manipulate spatial data

# Load custom functions
source("Functions.r")

# Load dispersal trajectories (csv)
disp <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv() %>%
  subset(State == "Disperser") %>%
  group_by(DogName) %>%
  nest() %>%

  # Create lines representing the observed dispersal tracks
  mutate(TrueTracks = map(data, function(x){
    coordinates(x) <- c("x", "y")
    x <- createSegments(x)
    x <- gLineMerge(x)
    crs(x) <- CRS("+init=epsg:4326")
    return(x)
  })) %>%

  # Create lines representing straight line dispersal track
  mutate(StraightTracks = map(data, function(x){
    x <- rbind(x[1, ], x[nrow(x), ])
    coordinates(x) <- c("x", "y")
    x <- createSegments(x)
    x <- gLineMerge(x)
    crs(x) <- CRS("+init=epsg:4326")
    return(x)
  }))

# Put tracks together
truetracks <- do.call(rbind, disp$TrueTracks)
straighttracks <- do.call(rbind, disp$StraightTracks)

# Load heatmaps
heat <- rast("03_Data/03_Results/99_Simulations.tif")
heat <- heat[[nlyr(heat)]]

# Crop to extent of dispersal events
heat <- crop(heat, truetracks, snap = "out")

# Visualize all
plot(sqrt(heat), col = viridis(50))
plot(truetracks, col = "red", add = T)
plot(straighttracks, col = "yellow", add = T)

# Extract values below all lines
true_values <- terra::extract(heat, vect(truetracks), touches = T) %>%
  lapply(unlist) %>%
  lapply(mean) %>%
  unlist()
straight_values <- terra::extract(heat, vect(straighttracks), touches = T) %>%
  lapply(unlist) %>%
  lapply(mean) %>%
  unlist()

# Put into single dataframe
values <- data.frame(True = true_values, Straight = straight_values) %>%
  gather(key = Type, value = Values)

# Visualize
ggplot(values, aes(x = Values, col = factor(Type))) + geom_boxplot()
