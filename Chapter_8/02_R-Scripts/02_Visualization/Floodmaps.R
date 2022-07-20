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
water <- rast("03_Data/02_CleanData/WaterCover.tif")
water <- water[[c(1, 3)]]

# Crop to dynamic area
flood <- dir(path = "03_Data/01_RawData/FLOODMAPS", pattern = ".tif$", full.names = T)[1]
flood <- rast(flood)
water <- crop(water, flood)

# Convert to dataframe
water <- as.data.frame(water, xy = T)
water <- pivot_longer(water, min:max, names_to = "FloodLevel", values_to = "Flooded")

# Plot
ggplot(water, aes(x = x, y = y, fill = as.factor(Flooded))) +
  geom_raster() +
  scale_fill_manual(values = c("white", "cornflowerblue"), name = "") +
  coord_equal() +
  theme(
      legend.position  = "none"
    , legend.box       = "vertical"
    , panel.background = element_blank()
    , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
    , axis.text.x      = element_text(angle = 45, hjust = 1)
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.4
    , line_width = 1.0
    , text_cex   = 1.0
    , height     = unit(0.1, "cm")
    # , bar_cols   = "black"
    , text_col   = "black"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.6, "cm"),
    , width    = unit(1.2, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 10
      )
  ) +
  facet_wrap(~ FloodLevel) +
  xlab("") +
  ylab("")
