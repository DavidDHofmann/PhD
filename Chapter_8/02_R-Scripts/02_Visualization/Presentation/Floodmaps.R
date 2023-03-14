################################################################################
#### Dynamic Floodmaps
################################################################################
# Description: Plot of the two floodmaps used

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
setwd(wd)

# Load required packages
library(sf)           # To handle spatial data
library(tidyverse)    # For data.wrangling
library(terra)        # To handle spatial data
library(ggplot2)      # For nice plots
library(ggspatial)    # To add scale bars and north arrows

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load stuff that we would like to plot
africa <- read_sf("03_Data/02_CleanData/Africa.shp")
water  <- rast("03_Data/02_CleanData/WaterCover.tif")
roads  <- read_sf("03_Data/02_CleanData/Roads.shp")
vills  <- read_sf("03_Data/02_CleanData/Villages.shp")
vills  <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y)

# Reference raster
r <- rast("03_Data/02_CleanData/ReferenceRaster.tif")
r <- as.data.frame(r, xy = T)

# We only care about min/max extents
water  <- water[[c("min", "max")]]
names(water) <- c("Min", "Max")

# Crop to area of interest
water <- crop(water, ext(c(21.5, 24, -20.5, -18)))

# To compute the flood extent, we need to focus on the Okavango delta, thus,
# let's load a shapefile of it
oka <- vect("03_Data/02_CleanData/MajorWaters.shp")
oka <- oka[oka$name == "Okavango Delta", ]

# Buffer it slightly
oka <- buffer(oka, width = 10000)
oka <- fillHoles(oka)

# This looks good. We now mask anything outside that polygon
water_masked <- mask(water, oka)
water_masked <- trim(water_masked)
water_masked <- subst(water_masked, 0, NA)

# Convert shape of okavango delta to sf
crs(oka) <- "+init=epsg:4236"
oka            <- st_as_sf(oka)
oka            <- rbind(oka, oka)
oka$Area       <- round(expanse(water_masked, unit = "km"))
oka$FloodLevel <- names(water_masked)
oka$FloodLevel <- factor(oka$FloodLevel, levels = c("Min", "Max"))
oka$Label      <- paste0("Area~flooded:~", oka$Area, "~km^2")

# Create facet labels
lab <- data.frame(
    FloodLevel = factor(c("Min", "Max"), levels = c("Min", "Max"))
  , Label      = c("a", "b")
  , x = c(20.75, 20.75)
  , y = c(-17.75, -17.75)
)

# Convert to dataframe
water <- as.data.frame(water, xy = T)
water <- pivot_longer(water, Min:Max, names_to = "FloodLevel", values_to = "Flooded")
water$FloodLevel <- factor(water$FloodLevel, levels = c("Min", "Max"))

# Create country labels
labels_countries <- data.frame(
    x     = c(25.5, 26, 25.7, 21.5, 23.5)
  , y     = c(-19.3, -18.2, -17.6, -17.6, -17.8)
  , Label = c("BOTSWANA", "ZIMBABWE", "ZAMBIA", "ANGOLA", "NAMIBIA")
)

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.1, 25.6)
  , y     = c(-19.1, -18.2, -17.5, -20.7)
  , Label = c("Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans")
)

# Plot
p1 <- ggplot() +
  geom_raster(data = water, aes(x = x, y = y, fill = as.factor(Flooded))) +
  geom_sf(data = roads, col = "gray70", lwd = 0.2) +
  geom_sf(data = oka, col = "white", fill = NA, linewidth = 0.5, linetype = 2) +
  geom_point(
      data        = subset(vills, place == "City")
    , mapping     = aes(x = x, y = y)
    , col         = "gray80"
    , shape       = 15
    , size        = 2
    , show.legend = F
  ) +
  geom_text(
      data     = subset(vills, place == "City")
    , mapping  = aes(x = x, y = y, label = name)
    , col      = "gray80"
    , fontface = 3
    , size     = 3
    , nudge_y  = c(0.1, -0.1, 0.1)
  ) +
  scale_fill_manual(values = c("transparent", "cornflowerblue"), name = "") +
  theme(
      legend.position  = "none"
    , legend.key       = element_rect(color = NA)
    , axis.text        = element_blank()
    , axis.ticks       = element_blank()
    , panel.border     = element_blank()
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , panel.background = element_blank()
    , plot.background  = element_blank()
    , strip.background = element_blank()
    , strip.text       = element_blank()
  ) +
  facet_wrap(~ FloodLevel) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(water$x), max(water$x))
    , ylim   = c(min(water$y), max(water$y))
    , expand = F
  ) +
  xlab("") +
  ylab("")

p2 <- ggplot() +
  geom_raster(data = water, aes(x = x, y = y, fill = as.factor(Flooded))) +
  geom_sf(data = roads, col = "gray70", lwd = 0.2) +
  geom_sf(data = oka, col = "white", fill = "white", linewidth = 0.5, linetype = 2, alpha = 0.2) +
  geom_point(
      data        = subset(vills, place == "City")
    , mapping     = aes(x = x, y = y)
    , col         = "gray80"
    , shape       = 15
    , size        = 2
    , show.legend = F
  ) +
  geom_text(
      data     = subset(vills, place == "City")
    , mapping  = aes(x = x, y = y, label = name)
    , col      = "gray80"
    , fontface = 3
    , size     = 3
    , nudge_y  = c(0.1, -0.1, 0.1)
  ) +
  scale_fill_manual(values = c("transparent", "cornflowerblue"), name = "") +
  theme(
      legend.position  = "none"
    , legend.key       = element_rect(color = NA)
    , axis.text        = element_blank()
    , axis.ticks       = element_blank()
    , panel.border     = element_blank()
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , panel.background = element_blank()
    , plot.background  = element_blank()
    , strip.background = element_blank()
    , strip.text       = element_blank()
  ) +
  facet_wrap(~ FloodLevel) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(water$x), max(water$x))
    , ylim   = c(min(water$y), max(water$y))
    , expand = F
  ) +
  xlab("") +
  ylab("")

# Store the plot
ggsave("05_Presentation/99_FloodExtent1.png"
  , plot   = p1
  , bg     = "transparent"
  , width  = 8
  , height = 3.5
  , scale  = 1.2
)
ggsave("05_Presentation/99_FloodExtent2.png"
  , plot   = p2
  , bg     = "transparent"
  , width  = 8
  , height = 3.5
  , scale  = 1.2
)
