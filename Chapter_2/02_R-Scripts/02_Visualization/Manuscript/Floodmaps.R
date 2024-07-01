################################################################################
#### Dynamic Floodmaps
################################################################################
# Description: Plot of the two floodmaps used

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(plyr)         # To round
library(sf)           # To handle spatial data
library(tidyverse)    # For data.wrangling
library(terra)        # To handle spatial data
library(ggplot2)      # For nice plots
library(ggspatial)    # To add scale bars and north arrows
library(ggpubr)       # To put plots together
library(latex2exp)    # For latex expressions in the plot

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load stuff that we would like to plot
# africa <- read_sf("03_Data/02_CleanData/Africa.shp")
water  <- rast("03_Data/02_CleanData/WaterCover.tif")
roads  <- read_sf("03_Data/02_CleanData/Roads.shp")
vills  <- read_sf("03_Data/02_CleanData/Villages.shp")
vills  <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y)

# Reference raster
r <- rast("03_Data/02_CleanData/ReferenceRaster.tif")
r <- crop(r, ext(21, 24.5, -20.8, -17.8))

# We only care about min/max extents
water  <- water[[c("min", "max")]]
water  <- crop(water, r)
names(water) <- c("Min", "Max")

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
crs(oka)       <- "+init=epsg:4236"
oka            <- st_as_sf(oka)
oka            <- rbind(oka, oka)
oka$Area       <- round_any(expanse(water_masked, unit = "km")$area, 500)
oka$FloodLevel <- factor(c("Min", "Max"), levels = c("Min", "Max"))

# Create a label of the area flooded
oka$Label <- paste0(
    "Area flooded: "
  , format(oka$Area, big.mark = ",")
  , " km$^2$"
)

# Convert to dataframe
r     <- as.data.frame(r, xy = T)
water <- as.data.frame(water, xy = T)
water <- pivot_longer(water, 3:4, names_to = "FloodLevel", values_to = "Flooded")
water <- mutate(water, FloodLevel = factor(FloodLevel, levels = c("Min", "Max")))

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 22.7)
  , y     = c(-19.1, -18.2, -20.4)
  , Label = c("Okavango\nDelta", "Linyanti\nSwamp", "Lake\nNgami")
)
labels_waters <- expand_grid(labels_waters, FloodLevel = c("Min", "Max"))
labels_waters <- mutate(labels_waters, FloodLevel = factor(FloodLevel, levels = c("Min", "Max")))

# We don't want to display lake ngami at minimum flood
labels_waters <- subset(labels_waters,
  !(Label == "Lake\nNgami" & FloodLevel == "Min")
)

# Create facetted plot
p <- ggplot() +
  geom_raster(data = water, aes(x = x, y = y, fill = as.factor(Flooded))) +
  geom_sf(data = roads, col = "gray70", lwd = 0.2) +
  geom_sf_text(
      data    = oka
    , mapping = aes(label = TeX(Label, output = "character"))
    , nudge_y = 0.7
    , nudge_x = -0.25
    , size    = 3
    , parse   = T
    , col     = "red"
  ) +
  geom_point(
      data        = subset(vills, place == "City")
    , mapping     = aes(x = x, y = y)
    , col         = "gray50"
    , shape       = 15
    , size        = 2
    , show.legend = F
  ) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("cornflowerblue", 1.8)
    , fontface = 3
    , size     = 2.5
  ) +
  geom_text(
      data     = subset(vills, place == "City")
    , mapping  = aes(x = x, y = y, label = name)
    , col      = "gray50"
    , fontface = 3
    , size     = 3
    , nudge_y  = c(0.1, -0.1, 0.1)
  ) +
  scale_fill_manual(values = c("white", "cornflowerblue"), name = "") +
  theme(
      legend.position  = "none"
    , legend.box       = "vertical"
    , panel.background = element_blank()
    , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 12
      )
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  theme(strip.background = element_rect(fill = "gray95", color = "transparent")) +
  xlab("") +
  ylab("") +
  facet_wrap(~ FloodLevel)

# Store the plot
ggsave("04_Manuscript/Figures/FloodExtent.png"
  , plot   = p
  , bg     = "white"
  , width  = 7.1
  , height = 3.5
  , scale  = 1.3
  , device = png
)
