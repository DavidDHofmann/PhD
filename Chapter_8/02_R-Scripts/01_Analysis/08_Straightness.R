################################################################################
#### Straightness
################################################################################
# Description: Computing the straightness of the simulations

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(tidyverse)      # To wrangle data
library(reproj)         # To quickly reproject coordinates

# Load dispersal simulations
sims <- read_rds("03_Data/03_Results/DispersalSimulation.rds")

# Keep only desired columns
sims <- sims[, c("x", "y", "sl_", "ta_", "TrackID", "StepNumber", "SourceArea", "FloodLevel")]

# # Subset to a random selection of paths for now
# back <- sims
# ids <- sample(unique(sims$TrackID), size = 1000, replace = F)
# sims <- subset(sims, TrackID %in% ids)

# To compute straightness, the coordinates need to be projected
sims[, c("x", "y")] <- reproj_xy(
    x      = sims[, c("x", "y")]
  , source = 4326
  , target = 32734
  # , source = "epsg:4326"
  # , target = "epsg:32734"
)

# Calcualte straightness of each simulated track
sims <- sims %>%
  nest(Data = -c(TrackID, SourceArea, FloodLevel)) %>%
  mutate(Straightness = map_dbl(Data, function(x) {

      # Compute cumulative distance
      cum_dist <- sum(x$sl_, na.rm = T)

      # Compute euclidean distance
      dx <- x$x[1] - x$x[length(x$x)]
      dy <- x$y[1] - x$y[length(x$y)]
      euc_dist <- sqrt(dx ** 2 + dy ** 2)

      # Compute straightness
      straightness <- euc_dist / cum_dist

      # Return it
      return(straightness)

  }))

# Remove the original data
strai <- select(sims, -Data)

# Store the extracted data to file
write_rds(strai, "03_Data/03_Results/Straightness.rds")

# Visualize
ggplot(data = strai, aes(x = FloodLevel, y = Straightness)) +
  geom_boxplot() +
  facet_wrap(~ SourceArea) +
  theme_minimal()
