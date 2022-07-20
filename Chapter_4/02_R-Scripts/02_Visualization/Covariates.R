################################################################################
#### Covariates
################################################################################
# Description: Plot of the covariates

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(tidyverse)     # For data wrangling
library(sf)            # For plotting
library(pbmcapply)     # For multicore abilities

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load observed movement data and covariates
obs <- read_csv("03_Data/01_RawData/ObservedMovements.csv")
covars <- stack("03_Data/01_RawData/CovariateLayers.grd")

# Extent in which animals were released
ext <- extent(c(50, 250, 50, 250))
ext <- as(ext, "SpatialPolygons")

################################################################################
#### Plot of Covariates and Simulated Trajectories
################################################################################
# Plot the covariates
p1 <- as.data.frame(covars, xy = T) %>%
  gather(key = covariate, value = value, 3:5) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = st_as_sf(ext), inherit.aes = F, fill = NA, col = "white", lty = 2) +
    scale_fill_viridis_c(option = "viridis") +
    coord_sf() +
    theme_minimal() +
    facet_wrap("covariate") +
    scale_x_continuous(breaks = seq(0, 300, by = 100)) +
    scale_y_continuous(breaks = seq(0, 300, by = 100)) +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
    theme(strip.background = element_rect(fill = "gray95", color = NA))

# Plot of trajectories
p2 <- ggplot(obs, aes(x = x, y = y, col = as.factor(ID))) +
  geom_sf(data = st_as_sf(ext), col = "cornflowerblue", fill = NA, inherit.aes = F, lwd = 1) +
  geom_path(size = 0.1) +
  geom_point(size = 0.1) +
  coord_sf(xlim = c(0, 300), ylim = c(0, 300)) +
  theme_minimal() +
  scale_color_viridis_d() +
  theme(
      legend.position = "none"
    , axis.title.y    = element_text(angle = 0, vjust = 0.5)
  )

# Store to file
ggsave(plot = p1, "04_Manuscript/99_Covariates.png", width = 6, height = 2.5, bg = "white", scale = 1.2)
ggsave(plot = p2, "04_Manuscript/99_Simulations.png", scale = 0.75, bg = "white")
