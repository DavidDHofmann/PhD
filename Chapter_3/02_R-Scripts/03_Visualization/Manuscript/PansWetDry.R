################################################################################
#### Plot of the same Area in the Wet and Dry Season
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(terra)
library(raster)
library(tidyverse)
library(tmaptools)
library(viridis)
library(ggpubr)

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Extent that I want to plot
exte <- ext(23.5, 23.8, -19.5, -19.2)

# Load "distance to pan" maps
dist_sep <- rast("03_Data/02_CleanData/00_Panmaps/DistanceTo/2019_09.tif")
dist_dec <- rast("03_Data/02_CleanData/00_Panmaps/DistanceTo/2019_12.tif")

# Crop them to the study extent
dist_sep <- crop(dist_sep, exte)
dist_dec <- crop(dist_dec, exte)

# Get a background satellite map
satmap <- read_osm(as(as.polygons(exte), "Spatial"), type = "bing", zoom = 12)
satmap <- as(satmap, "SpatRaster")
plotRGB(satmap)

# Pick a few pans to visualize
pans_sep <- rbind(
    c(23.56439, -19.29996)
  , c(23.5838, -19.35314)
  , c(23.74708, -19.34322)
) %>% as.data.frame() %>% setNames(c("x", "y")) %>% mutate(ID = 1:3) %>% mutate(Date = "2019-09")
pans_dec <- rbind(
    c(23.63163, -19.30857)
  , c(23.56863, -19.38096)
  , c(23.6665, -19.40605)
) %>% as.data.frame() %>% setNames(c("x", "y")) %>% mutate(ID = 4:6) %>% mutate(Date = "2019-12")
pans <- rbind(pans_sep, pans_dec)

# Get satellite images for those pans
sat_sep <- lapply(1:nrow(pans_sep), function(i) {
  sat <- read_osm(as(as.polygons(ext(vect(as.matrix(pans_sep[i, c("x", "y")]))) + 0.005), "Spatial"), type = "bing", zoom = 16)
  sat <- as(sat, "SpatRaster")
  sat <- project(sat, "epsg:4326", method = "near")
  return(sat)
})

sat_dec <- lapply(1:nrow(pans_dec), function(i) {
  sat <- read_osm(as(as.polygons(ext(vect(as.matrix(pans_dec[i, c("x", "y")]))) + 0.005), "Spatial"), type = "bing", zoom = 16)
  sat <- as(sat, "SpatRaster")
  sat <- project(sat, "epsg:4326", method = "near")
  return(sat)
})

# Convert stuff to dataframes
dist_sep <- as.data.frame(dist_sep, xy = T) %>% mutate(Date = "2019-09")
dist_dec <- as.data.frame(dist_dec, xy = T) %>% mutate(Date = "2019-12")

# Put everything together
dist <- rbind(dist_sep, dist_dec)

# Visualize
p_dist <- ggplot() +
  geom_raster(data = dist, aes(x = x, y = y, fill = class)) +
  geom_point(data = pans, aes(x = x, y = y), col = "red") +
  geom_text(data = pans, aes(x = x, y = y, label = ID), col = "red", nudge_y = 0.01) +
  scale_fill_viridis_c(name = "Distance (m)", trans = "sqrt") +
  theme_awesome() +
  facet_wrap(~ Date) +
  coord_sf()

p_sep <- lapply(1:length(sat_sep), function(i) {
  dat <- as.data.frame(sat_sep[[i]], xy = T)
  p <- ggplot(dat, aes(x = x, y = y, fill = rgb(red, green, blue, maxColorValue = 255))) +
    geom_raster() +
    scale_fill_identity() +
    coord_sf() +
    theme_void()
  return(p)
})

p_dec <- lapply(1:length(sat_sep), function(i) {
  dat <- as.data.frame(sat_dec[[i]], xy = T)
  p <- ggplot(dat, aes(x = x, y = y, fill = rgb(red, green, blue, maxColorValue = 255))) +
    geom_raster() +
    scale_fill_identity() +
    coord_sf() +
    theme_void()
  return(p)
})

# Put all together
p_sep <- ggarrange(p_sep[[1]], p_sep[[2]], p_sep[[3]], nrow = 1, labels = c("1", "2", "3"))
p_dec <- ggarrange(p_dec[[1]], p_dec[[2]], p_dec[[3]], nrow = 1, labels = c("4", "5", "6"))
ggarrange(p_sep, p_dist, p_dec, nrow = 3, heights = c(1, 2, 1))
