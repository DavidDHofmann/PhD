################################################################################
#### Betweenness
################################################################################
# Description: Visualization of betweenness

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

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load the betwenness maps
maps <- read_rds("03_Data/03_Results/Betweenness.rds")
maps <- subset(maps, Steps %in% c(500, 1000, 2000))

# Load shapefiles that we want to plot
area  <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
roads <- read_sf("03_Data/02_CleanData/Roads.shp")
afric <- read_sf("03_Data/02_CleanData/Africa.shp")
vills <- read_sf("03_Data/02_CleanData/Villages.shp")
vills <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y)

# Reference raster
r <- rast("03_Data/02_CleanData/ReferenceRaster.tif")
# r <- crop(r, extent(22, 25, -20.5, -18.5))

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

# # Apply focal filter to buffer/smooth maps
# maps <- mutate(maps, betweenness = map(betweenness, function(x) {
#   x <- focal(x, w = matrix(1, 3, 3), fun = mean)
#   return(x)
# }))

# Convert maps to dataframes
maps <- maps %>%
  mutate(betweenness = map(betweenness, function(x) {
    result <- as.data.frame(x, xy = T)
    names(result) <- c("x", "y", "Betweenness")
    return(result)
  })) %>% unnest(betweenness)

# Make sure the levels are correctly ordered
maps$FloodLevel <- factor(maps$FloodLevel, levels = c("Min", "Mean", "Max"))

# Convert the reference raster to a dataframe
r <- as.data.frame(r, xy = T)

# Visualize the betweenness
p1 <- ggplot() +
  geom_raster(
      data    = maps
    , mapping = aes(x = x, y = y, fill = Betweenness)
  ) +
  geom_sf(data = roads, col = "gray50", lwd = 0.1) +
  geom_sf(
      data        = area
    , col         = "white"
    , fill        = NA
    , lty         = 1
    , lwd         = 0.1
    , show.legend = F
    , alpha       = 0.15
  ) +
  geom_sf(data = afric, lwd = 0.4, col = "white", fill = NA) +
  geom_sf_text(
      data          = area
    , mapping       = aes(label = ID)
    , size          = 1.5
    , check_overlap = T
    , col           = "white"
  ) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "gray50"
    , fontface = 3
    , size     = 1.3
  ) +
  geom_text(
      data     = labels_countries
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "white"
    , size     = 2
    , fontface = 2
  ) +
  scale_fill_gradientn(
      colours = viridis::magma(100)
    , labels  = function(x){format(x, big.mark = "'")}
    , trans   = "sqrt"
    , guide   = guide_colorbar(
      , title          = "Betweenness"
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
      legend.position  = "bottom"
    , legend.box       = "vertical"
    , panel.background = element_blank()
    , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
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
  facet_grid(FloodLevel ~ Steps)

# Let's also generate a map where we only consider 2000 steps
maps_sub <- subset(maps, Steps == 2000 & FloodLevel != "Mean")
maps_sub$FloodLevel <- droplevels(maps_sub$FloodLevel)
p2 <- ggplot() +
  geom_raster(
      data    = maps_sub
    , mapping = aes(x = x, y = y, fill = Betweenness)
  ) +
  geom_sf(data = roads, col = "gray50", lwd = 0.1) +
  geom_point(
      data        = vills
    , mapping     = aes(x = x, y = y, size = place)
    , col         = "gray30"
    , shape       = 15
    , show.legend = F
  ) +
  geom_sf(
      data        = area
    , col         = "white"
    , fill        = "white"
    , lty         = 1
    , lwd         = 0.1
    , show.legend = F
    , alpha       = 0.15
  ) +
  geom_sf(data = afric, lwd = 0.4, col = "white", fill = NA) +
  geom_sf_text(
      data          = area
    , mapping       = aes(label = ID)
    , size          = 1.5
    , col           = "white"
  ) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "gray50"
    , fontface = 3
    , size     = 1.3
  ) +
  geom_text(
      data     = labels_countries
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "white"
    , size     = 2
    , fontface = 2
  ) +
  geom_text(
      data     = subset(vills, place == "City")
    , mapping  = aes(x = x, y = y, label = name)
    , col      = "gray30"
    , fontface = 3
    , size     = 2
    , nudge_y  = c(0.1, -0.1, 0.1)
  ) +
  scale_size_manual(values = c(1.0, 0.25)) +
  scale_fill_gradientn(
      colours = viridis::magma(100)
    , labels  = function(x){format(x, big.mark = "'")}
    , trans   = "sqrt"
    , guide   = guide_colorbar(
      , title          = "Betweenness"
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
      legend.position  = "bottom"
    , legend.box       = "vertical"
    , panel.background = element_blank()
    , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
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
  facet_grid(~ FloodLevel)

################################################################################
#### Store Plot
################################################################################
# Store the plot
ggsave("04_Manuscript/99_Betweenness.png"
  , plot   = p1
  , bg     = "white"
  , width  = 8
  , height = 5
  , scale  = 1
)

# Store the second plot
ggsave("04_Manuscript/99_Betweenness2000.png"
  , plot   = p2
  , bg     = "white"
  , width  = 8
  , height = 4
  , scale  = 1
)
