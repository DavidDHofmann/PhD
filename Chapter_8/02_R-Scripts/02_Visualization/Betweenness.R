################################################################################
#### Betweenness
################################################################################
# Description: Visualization of betweenness

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(terra)          # To handle spatial data
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(ggspatial)      # For scale bar and north arrow
library(rgdal)          # To handle spatial data
library(sf)             # To handle spatial data

# Load the betwenness maps
maps <- read_rds("03_Data/03_Results/Betweenness.rds")
maps <- subset(maps, Flood != "Mean" & Steps %in% c(500, 1000, 2000))

# Load shapefile of source areas
area <- read_sf("03_Data/02_CleanData/SourceAreas.shp")

# Reference raster
r <- rast("03_Data/02_CleanData/ReferenceRaster.tif")
# r <- crop(r, extent(22, 25, -20.5, -18.5))

# Merge heatmaps of the same number of steps and flood extent
maps <- maps %>%
  mutate(betweenness = map(betweenness, function(x) {
    result <- as.data.frame(x, xy = T)
    names(result) <- c("x", "y", "Betweenness")
    return(result)
  })) %>% unnest(betweenness)

# Convert the reference raster to a dataframe
r <- as.data.frame(r, xy = T)

# Visualize the maps
ggplot() +
  geom_raster(
      data    = maps
    , mapping = aes(x = x, y = y, fill = Betweenness)
  ) +
  geom_sf(
      data        = area
    , col         = "white"
    , fill        = NA
    , lty         = 1
    , lwd         = 0.1
    , show.legend = F
    , alpha       = 0.6
  ) +
  geom_sf_text(
      data          = area
    , mapping       = aes(label = Name)
    , size          = 1.4
    , check_overlap = T
    , col           = "white"
  ) +
  scale_fill_gradientn(
      colours = viridis::magma(100)
    , labels  = function(x){format(x, big.mark = "'")}
    # , trans   = "sqrt"
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
    , axis.text.x      = element_text(angle = 45, hjust = 1)
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
    , width    = unit(0.4, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("white", "white")
        , line_col  = NA
        , text_col  = "white"
        , text_size = 4
      )
  ) +
  facet_grid(Flood ~ Steps)
