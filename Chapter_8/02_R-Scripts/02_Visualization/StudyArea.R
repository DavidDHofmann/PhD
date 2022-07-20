################################################################################
#### Study Area
################################################################################
# Description: Plots of the study area

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
setwd(wd)

# Load required packages
library(sf)           # To handle spatial data
library(ggplot2)      # For nice plots
library(ggspatial)    # To add scale bars and north arrows

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load stuff that we would like to plot
africa <- read_sf("03_Data/02_CleanData/Africa.shp")
prot   <- read_sf("03_Data/02_CleanData/Protected.shp")
water  <- read_sf("03_Data/02_CleanData/MajorWaters.shp")
areas  <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
roads  <- read_sf("03_Data/02_CleanData/Roads.shp")

# Reorder the levels of national parks
prot$Desig <- factor(prot$Desig, levels = c("National Park", "Forest Reserve", "Protected"))

# Get an extent object to show on the map
extent <- st_as_sfc(st_bbox(prot))

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.2, 25.6)
  , y     = c(-19, -18.2, -17.6, -20.7)
  , Label = c("Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans")
)

# Create labels for some national parks and source areas
labels_nationalparks <- data.frame(
    x = c(26.56, 22.35, 23.71, 24.51, 24.7)
  , y = c(-19.08, -17.60, -19.29, -18.65, -20.4)
  , Label = paste0(c("Hwange", "Luengue-Luiana", "Moremi", "Chobe", "Nxai Pan"), "\nNP")
)

# Prepare plot of africa
p1 <- ggplot() +
  geom_sf(data = africa, fill = "gray90", col = "white", lwd = 0.1) +
  geom_sf(data = extent, fill = "red", col = "red", alpha = 0.2) +
  theme_minimal() +
  theme(
      panel.grid = element_blank()
    , axis.ticks = element_blank()
    , axis.text  = element_blank()
  )

# Prepare plot of study area
p2 <- ggplot() +
  geom_sf(data = prot, aes(fill = Desig), col = NA, alpha = 0.7) +
  geom_sf(data = water, fill = "cornflowerblue", col = NA) +
  geom_sf(data = areas, col = "orange", fill = "orange", alpha = 0.2) +
  geom_sf(data = roads, col = "gray70", lwd = 0.2) +
  geom_sf(data = extent, col = "red", fill = NA) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("cornflowerblue", 1.4)
    , fontface = 3
  ) +
  geom_text(
      data     = labels_nationalparks
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("green", 1.8)
    , fontface = 2
  ) +
  scale_fill_brewer(palette = "Greens", direction = -1, name = "") +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
    , pad_x      = unit(0.9, "cm")
    , pad_y      = unit(0.8, "cm")
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , pad_x    = unit(0.7, "cm")
    , pad_y    = unit(0.6, "cm")
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 12
      )
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  xlab("") +
  ylab("")

# Show the plots
p1
p2
