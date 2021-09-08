################################################################################
#### Analysing Simulated Data
################################################################################
# Description: Use the simulated data to generate connectivity metrics. We then
# compare the metrics resulting from the different landscapes, i.e. with or
# without buffer etc. To compare maps, we will use pearson's correlation
# coefficient, bhattacharyya's affinity, and the root mean squared error

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(sf)            # For plotting spatial features
library(rgeos)         # For manipulating spatial objects
library(ggpubr)        # To put plots together
library(spatstat)      # To rasterize quickly
library(maptools)      # To rasterize quickly

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_6")

# Load covariates and other necessary data from previous session
load("03_Data/Landscape.Rdata")
load("03_Data/Simulations.Rdata")

################################################################################
#### Useful Functions
################################################################################
# Function to quickly rasterize lines using spatstat
rasterizeSpatstat <- function(l, r){
  values(r) <- 0
  im <- as.im.RasterLayer(r)
  summed <- im
  for (y in 1:length(l)){
    line    <- as.psp(l[y, ], window = im)
    line    <- as.mask.psp(line)
    line_r  <- as.im.owin(line, na.replace = 0)
    summed  <- Reduce("+", list(summed, line_r))
  }
  return(raster(summed))
}

# Function to calculate the mean squared error
rmse <- function(x, y){sqrt(mean((x - y) ** 2))}

# Function to calculate bhattacharyya's distance
bhattacharyya <- function(x, y){
  m1 <- as.matrix(x)
  m2 <- as.matrix(y)
  mean_m1 <- mean(m1)
  mean_m2 <- mean(m2)
  mean_difference <- mean_m1 - mean_m2
  cov_m1 <- cov(m1)
  cov_m2 <- cov(m2)
  p <- (cov_m1 + cov_m2) / 2
  bh_distance <- 0.125 * t (mean_difference) * p ** (-1) * mean_difference +
    0.5 * log(det(p) / sqrt(det(cov_m1) * det(cov_m2)))
  bh_distance <- as.vector(bh_distance)
  return(bh_distance)
}

################################################################################
#### Generate Heatmap
################################################################################
# Generate heatmap for each set of simulations
design$Heatmap <- lapply(1:nrow(design), function(x){

  # Extract simulations
  sims <- design$Simulations[[x]]

  # Coerce simulated coordinates to tracks
  tracks <- lapply(unique(sims$ID), function(x){
    sub <- sims[sims$ID == x, ]
    coordinates(sub) <- c("x", "y")
    lines <- spLines(sub)
    return(lines)
  })
  tracks <- do.call(rbind, tracks)
  tracks$ID <- 1:length(tracks)

  # Rasterize and count the tracks
  heatmap <- rasterizeSpatstat(tracks, covars$Layers[[1]][[1]])

  # Return it
  return(heatmap)
})

################################################################################
#### Inter-Patch Connectivity
################################################################################
design$Interpatch <- lapply(1:nrow(design), function(x){

  # Extract simulations
  sims <- design$Simulations[[x]]

  # Coerce simulated coordinates to tracks
  tracks <- lapply(unique(sims$ID), function(x){
    sub <- sims[sims$ID == x, ]
    coordinates(sub) <- c("x", "y")
    lines <- spLines(sub)
    return(lines)
  })
  tracks <- do.call(rbind, tracks)
  tracks$ID <- 1:length(tracks)

  # Get first coordinate of each simulation
  first <- subset(sims, step_number == 1)
  coordinates(first) <- c("x", "y")

  # Determine with which national park each startpoint intersects (can only be
  # one)
  from <- gIntersects(nps, first, byid = T)
  from <- apply(from, 1, which)
  from <- as.vector(from)

  # Let's also check whith which national parks each trajectory intersects
  to <- gIntersects(nps, tracks, byid = T)
  to <- as.data.frame(to)
  names(to) <- 1:ncol(to)

  # Put all into a single dataframe
  inter <- cbind(first$ID, from, to)
  names(inter)[1:2] <- c("ID", "from")

  # Let's count the number of connections from one park to another
  conns <- inter %>%
    gather(key = to, value = reached, 3:ncol(.)) %>%
    subset(reached & from != to) %>%
    group_by(from, to) %>%
    summarize(total_connections = n(), .groups = "drop")

  # Return the connections
  return(conns)

})

################################################################################
#### Compare Heatmaps
################################################################################
# We need to compare maps of the sampe replicate so it makes sense to spread the
# tibble again
design <- design %>%
  dplyr::select(-c(Simulations, Interpatch)) %>%
  spread(Covariates, Heatmap)

# Specify the comparisons that we want to make
comps <- rbind(
    c("Real", "Real")
  , c("Randomized", "Real")
  , c("Cropped", "Real")
)

# Compare maps using the three different metrics
design$Comparisons <- lapply(1:nrow(design), function(x){
  comparisons <- lapply(1:nrow(comps), function(y){
    z1 <- values(unlist(design[x, ][, comps[y, 1]])[[1]])
    z2 <- values(unlist(design[x, ][, comps[y, 2]])[[1]])
    corr <- cor(z1, z2)
    bhat <- bhattacharyya(z1, z2)
    erro <- rmse(z1, z2)
    return(cbind(corr, bhat, erro))
  })
  comparisons <- do.call(rbind, comparisons)
  comparisons <- cbind(comps, comparisons)
  return(comparisons)
})
