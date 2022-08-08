################################################################################
#### Heatmaps
################################################################################
# Description: Visualization of the heatmaps

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
wd <- "C:/Users/david/switchdrive/University/15. PhD/Chapter_8"
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

# Load the heatmaps
maps <- read_rds("03_Data/03_Results/Heatmaps.rds")
maps <- subset(maps, Steps %in% c(500, 1000, 2000))

# Load shapefiles that we want to plot
area  <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
roads <- read_sf("03_Data/02_CleanData/Roads.shp")
afric <- read_sf("03_Data/02_CleanData/Africa.shp")

# Load the reference raster
r <- raster("03_Data/02_CleanData/ReferenceRaster.tif")
# r <- crop(r, extent(22, 25, -20.5, -18.5))

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.1, 25.6)
  , y     = c(-19, -18.2, -17.5, -20.7)
  , Label = c("Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans")
)

# Merge heatmaps of the same number of steps and flood extent
maps <- maps %>%
  nest(Data = -c(Steps, Flood)) %>%
  mutate(Data = map(Data, function(x) {
    result <- sum(stack(x$heatmap))
    result <- crop(result, r)
    return(result)
  }))

# Compute difference maps
diffs <- lapply(unique(maps$Steps), function(x) {
    submaps <- subset(maps, Steps == x & Flood != "Mean")
    submaps <- stack(submaps$Data)
    diff <- submaps[[2]] - submaps[[1]]
    names(diff) <- x
    return(diff)
  }) %>%
  stack() %>%
  as.data.frame(xy = T) %>%
  pivot_longer(X500:X2000, names_to = "Steps", values_to = "Difference") %>%
  mutate(Steps = as.numeric(substr(Steps, start = 2, stop = nchar(Steps))))

# Convert the maps to dataframes
maps <- maps %>%
  mutate(Data = map(Data, function(x) {
    result <- as.data.frame(x, xy = T)
    names(result) <- c("x", "y", "Heat")
    return(result)
  })) %>% unnest(Data)

# Make sure the levels are correctly ordered
maps$Flood <- factor(maps$Flood, levels = c("Min", "Mean", "Max"))

# Also convert the reference raster to a dataframe
r <- as.data.frame(r, xy = T)

# Visualize the difference maps
cols <- colorRampPalette(c("orange", "white", "cornflowerblue"))
p1 <- ggplot() +
  geom_raster(
      data    = diffs
    , mapping = aes(x = x, y = y, fill = Difference)
  ) +
  geom_sf(data = roads, col = "gray70", lwd = 0.1) +
  geom_sf(
      data        = area
    , col         = "black"
    , fill        = NA
    , lty         = 1
    , lwd         = 0.1
    , show.legend = F
    , alpha       = 0.6
  ) +
  geom_sf(data = afric, lwd = 0.2, col = "black", fill = NA) +
  geom_sf_text(
      data          = area
    , mapping       = aes(label = Name)
    , size          = 1.5
    , col           = "black"
    , nudge_y       = c(0.1, -0.1, -0.1, 0.2, 0)
  ) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("cornflowerblue", 1.4)
    , fontface = 3
    , size     = 1.3
  ) +
  scale_fill_gradientn(
      colours = cols(100)
    , limits  = c(-250, +250)
    , labels  = function(x){format(x, big.mark = "'")}
    , oob     = squish
    , guide   = guide_colorbar(
      , title          = TeX(r'(#$Trajectories_{Flood = min}$ - #$Trajectories_{Flood = max}$)')
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
    , bar_cols   = c("black", "white")
    , text_col   = "black"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(0.7, "cm"),
    , width    = unit(0.6, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 4
      )
  ) +
  facet_grid(~ Steps)

# Visualize the difference maps for 2000 steps
diffs_sub <- subset(diffs, Steps == 2000)
p2 <- ggplot() +
  geom_raster(
      data    = diffs_sub
    , mapping = aes(x = x, y = y, fill = Difference)
  ) +
  geom_sf(data = roads, col = "gray70", lwd = 0.1) +
  geom_sf(
      data        = area
    , col         = "black"
    , fill        = NA
    , lty         = 1
    , lwd         = 0.1
    , show.legend = F
    , alpha       = 0.6
  ) +
  geom_sf(data = afric, lwd = 0.2, col = "black", fill = NA) +
  geom_sf_text(
      data          = area
    , mapping       = aes(label = Name)
    , size          = 1.5
    , col           = "black"
    , nudge_y       = c(0.1, -0.1, -0.1, 0.2, 0)
  ) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("cornflowerblue", 1.4)
    , fontface = 3
    , size     = 1.3
  ) +
  scale_fill_gradientn(
      colours = cols(100)
    , limits  = c(-250, +250)
    , labels  = function(x){format(x, big.mark = "'")}
    , oob     = squish
    , guide   = guide_colorbar(
      , title          = TeX(r'(#$Trajectories_{Flood = min}$ - #$Trajectories_{Flood = max}$)')
      , show.limits    = T
      , title.position = "bottom"
      , title.hjust    = 0.5
      , ticks          = F
      , barheight      = unit(0.2, "cm")
      , barwidth       = unit(8.0, "cm")
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
    , bar_cols   = c("black", "white")
    , text_col   = "black"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(0.7, "cm"),
    , width    = unit(0.6, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 4
      )
  )

# Visualize the heatmaps
p3 <- ggplot() +
  geom_raster(
      data    = maps
    , mapping = aes(x = x, y = y, fill = Heat)
  ) +
  geom_sf(data = roads, col = "gray70", lwd = 0.1) +
  geom_sf(
      data        = area
    , col         = "white"
    , fill        = NA
    , lty         = 1
    , lwd         = 0.1
    , show.legend = F
    , alpha       = 0.6
  ) +
  geom_sf(data = afric, lwd = 0.2, col = "gray50", fill = NA) +
  geom_sf_text(
      data          = area
    , mapping       = aes(label = Name)
    , size          = 1.5
    , col           = "white"
    , nudge_y       = c(0.1, -0.1, -0.1, 0.2, 0)
  ) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "gray50"
    , fontface = 3
    , size     = 1.3
  ) +
  scale_fill_gradientn(
      colours = viridis::magma(100)
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
  facet_grid(Flood ~ Steps)

# Let's also generate a map where we only consider 2000 steps
maps_sub <- subset(maps, Steps == 2000 & Flood != "Mean")
maps_sub$Flood <- droplevels(maps_sub$Flood)
p4 <- ggplot() +
  geom_raster(
      data    = subset(maps_sub)
    , mapping = aes(x = x, y = y, fill = Heat)
  ) +
  geom_sf(data = roads, col = "gray70", lwd = 0.1) +
  geom_sf(
      data        = area
    , col         = "white"
    , fill        = NA
    , lty         = 1
    , lwd         = 0.1
    , show.legend = F
    , alpha       = 0.6
  ) +
  geom_sf(data = afric, lwd = 0.2, col = "gray50", fill = NA) +
  geom_sf_text(
      data          = area
    , mapping       = aes(label = Name)
    , size          = 1.5
    , col           = "white"
    , nudge_y       = c(0.1, -0.1, -0.1, 0.2, 0)
  ) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "gray50"
    , fontface = 3
    , size     = 1.3
  ) +
  scale_fill_gradientn(
      colours = viridis::magma(100)
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
  facet_grid(~ Flood)

################################################################################
#### Store Plots
################################################################################
# Store the first plot
ggsave("04_Manuscript/99_DifferenceHeatmaps.png"
  , plot   = p1
  , bg     = "white"
  , width  = 8
  , height = 3
  , scale  = 1.4
)

# Store the first plot
ggsave("04_Manuscript/99_DifferenceHeatmaps2000.png"
  , plot   = p2
  , bg     = "white"
  , width  = 4
  , height = 3
  , scale  = 1.4
)

# Store the second plot
ggsave("04_Manuscript/99_Heatmaps.png"
  , plot   = p3
  , bg     = "white"
  , width  = 8
  , height = 6
  , scale  = 1.2
)

# Store the second plot
ggsave("04_Manuscript/99_Heatmaps2000.png"
  , plot   = p4
  , bg     = "white"
  , width  = 8
  , height = 4
  , scale  = 1
)
