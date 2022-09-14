################################################################################
#### Heatmaps
################################################################################
# Description: Visualization of the heatmaps

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
setwd(wd)

# Load required packages
library(terra)          # To handle spatial data
library(raster)         # To handle spatial data
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(ggspatial)      # For scale bar and north arrow
library(rgdal)          # To handle spatial data
library(sf)             # To handle spatial data
library(scales)         # To squish oob values
library(latex2exp)      # For easy latex code
library(RColorBrewer)   # For custom colors

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load shapefiles that we want to plot
area  <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
roads <- read_sf("03_Data/02_CleanData/Roads.shp")
afric <- read_sf("03_Data/02_CleanData/Africa.shp")
vills <- read_sf("03_Data/02_CleanData/Villages.shp")
vills <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y)

# Prepare a custom color ramp (I don't like ggplots version of spectral)
spectral <- colorRampPalette(rev(brewer.pal(11, name = "Spectral")))

# Load the reference raster
r <- raster("03_Data/02_CleanData/ReferenceRaster.tif")
r <- as.data.frame(r, xy = T)

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

################################################################################
#### By Area
################################################################################
# Load the heatmaps and keep only desired columns
maps <- "03_Data/03_Results/HeatmapsBetweennessLocal.rds" %>%
  read_rds() %>%
  select(Steps, SourceArea, FloodLevel, Data = Heatmap) %>%
  subset(Steps == 2000)

# Convert the maps to dataframes
maps <- maps %>%
  mutate(Data = map(Data, function(x) {
    result <- as.data.frame(x, xy = T)
    names(result) <- c("x", "y", "Heat")
    return(result)
  })) %>% unnest(Data)

# Make sure the levels are correctly ordered
maps$FloodLevel <- factor(maps$FloodLevel, levels = c("Min", "Mean", "Max"))

# Plot
p1 <- ggplot() +
  geom_raster(
      data    = maps
    , mapping = aes(x = x, y = y, fill = Heat)
  ) +
  geom_sf(data = afric, lwd = 0.4, col = "black", fill = NA) +
  scale_fill_gradientn(
      colors  = spectral(100)
    , labels  = function(x){format(x, big.mark = "'")}
    , guide   = guide_colorbar(
      , title          = "#Traversing Trajectories"
      , show.limits    = T
      , title.position = "bottom"
      , title.hjust    = 0.5
      , ticks          = F
      , barheight      = unit(0.2, "cm")
      , barwidth       = unit(16.0, "cm")
    )
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
  ) +
  theme(
      legend.position    = "bottom"
    , legend.box         = "vertical"
    , panel.background   = element_blank()
    , panel.border       = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , text_cex   = 0.5
    , height     = unit(0.1, "cm")
    , bar_cols   = "white"
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(0.7, "cm"),
    , width    = unit(0.6, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("white", "white")
        , line_col  = NA
        , text_col  = "white"
        , text_size = 4
      )
  ) +
  facet_grid(FloodLevel ~ SourceArea)

################################################################################
#### By Steps
################################################################################
# Load the heatmaps and keep only desired columns
maps <- "03_Data/03_Results/HeatmapsBetweennessGlobal.rds" %>%
  read_rds() %>%
  select(Steps, FloodLevel, Data = Heatmap)

# Convert the maps to dataframes
maps <- maps %>%
  mutate(Data = map(Data, function(x) {
    result <- as.data.frame(x, xy = T)
    names(result) <- c("x", "y", "Heat")
    return(result)
  })) %>% unnest(Data)

# Make sure the levels are correctly ordered
maps$FloodLevel <- factor(maps$FloodLevel, levels = c("Min", "Mean", "Max"))

# Plot by steps
p2 <- ggplot() +
  geom_raster(
      data    = maps
    , mapping = aes(x = x, y = y, fill = Heat)
  ) +
  geom_sf(data = afric, lwd = 0.4, col = "black", fill = NA) +
  scale_fill_gradientn(
      colors  = spectral(100)
    , labels  = function(x){format(x, big.mark = "'")}
    , guide   = guide_colorbar(
      , title          = ""
      , show.limits    = T
      , title.position = "bottom"
      , title.hjust    = 0.5
      , ticks          = F
      , barheight      = unit(0.2, "cm")
      , barwidth       = unit(16.0, "cm")
    )
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
  ) +
  theme(
      legend.position    = "bottom"
    , legend.box         = "vertical"
    , panel.background   = element_blank()
    , panel.border       = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , text_cex   = 0.5
    , height     = unit(0.1, "cm")
    , bar_cols   = "white"
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(0.7, "cm"),
    , width    = unit(0.6, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("white", "white")
        , line_col  = NA
        , text_col  = "white"
        , text_size = 4
      )
  ) +  facet_grid(FloodLevel ~ Steps)

################################################################################
#### Store Plots
################################################################################
# Store the plot
ggsave("04_Manuscript/99_HeatmapsByAreas.png"
  , plot   = p1
  , bg     = "white"
  , height = 6
  , width  = 20
)

# Store the plot
ggsave("04_Manuscript/99_HeatmapsBySteps.png"
  , plot   = p2
  , bg     = "white"
  , width  = 8
  , height = 6
  , scale  = 1.2
)
