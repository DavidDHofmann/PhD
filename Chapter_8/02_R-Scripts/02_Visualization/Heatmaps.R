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
library(ggpubr)         # To arrange multiple plots

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load shapefiles that we want to plot
areas <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
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

# Generate area labels
labels_areas <- st_coordinates(st_point_on_surface(areas))
labels_areas <- cbind(labels_areas, st_drop_geometry(areas))

# Load the heatmaps and keep only desired columns
maps1 <- "03_Data/03_Results/HeatmapsBetweennessGlobal.rds" %>%
  read_rds() %>%
  mutate(SourceArea = NA, Level = "Global") %>%
  select(Steps, SourceArea, FloodLevel, Level, Data = Heatmap)
maps2 <- "03_Data/03_Results/HeatmapsBetweennessLocal.rds" %>%
  read_rds() %>%
  mutate(Level = "Local") %>%
  select(Steps, SourceArea, FloodLevel, Level, Data = Heatmap)
maps <- rbind(maps1, maps2)

# Convert the maps to dataframes
maps <- maps %>%
  mutate(Data = map(Data, function(x) {
    result <- as.data.frame(x, xy = T)
    names(result) <- c("x", "y", "Heat")
    return(result)
  })) %>% unnest(Data)

# Make sure the levels are correctly ordered
maps$FloodLevel <- factor(maps$FloodLevel, levels = c("Min", "Max"))

# Function to plot
plotHeatmap <- function(data, formula = ~FloodLevel, barwidth = 16) {
  ggplot() +
    geom_raster(
        data    = data
      , mapping = aes(x = x, y = y, fill = Heat)
    ) +
    geom_sf(data = roads, col = "gray90", lwd = 0.1) +
    geom_point(
        data        = subset(vills, place == "City")
      , mapping     = aes(x = x, y = y, size = place)
      , col         = "gray90"
      , shape       = 15
      , show.legend = F
      , size        = 1
    ) +
    geom_sf(
        data        = areas
      , col         = "black"
      , fill        = "black"
      , lty         = 1
      , lwd         = 0.1
      , show.legend = F
      , alpha       = 0.15
    ) +
    geom_sf(data = afric, lwd = 0.4, col = "black", fill = NA) +
    geom_text(
        data     = labels_waters
      , mapping  = aes(x = x, y = y, label = Label)
      , col      = "gray20"
      , fontface = 3
      , size     = 1.3
    ) +
    geom_text(
        data     = labels_countries
      , mapping  = aes(x = x, y = y, label = Label)
      , col      = "black"
      , size     = 2
      , fontface = 2
    ) +
    geom_text(
        data     = labels_areas
      , mapping  = aes(x = X, y = Y, label = ID)
      , col      = "black"
      , fontface = 3
      , size     = 2
    ) +
    geom_text(
        data     = subset(vills, place == "City")
      , mapping  = aes(x = x, y = y, label = name)
      , col      = "gray90"
      , fontface = 3
      , size     = 2
      , nudge_y  = c(0.1, -0.1, 0.1)
    ) +
    scale_size_manual(values = c(1.0, 0.25)) +
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
        , barwidth       = unit(barwidth, "cm")
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
      , bar_cols   = c("white", "white")
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
    facet_grid(formula)
}

# Apply it
p1 <- plotHeatmap(subset(maps, Level == "Global" & Steps == 2000))

# I also want to generate source-area specific plots
p2 <- list()
for (i in sort(unique(maps$SourceArea))) {
  p2[[i]] <- plotHeatmap(
      data     = subset(maps, Level == "Local" & Steps == 2000 & SourceArea == i)
    , barwidth = 12
  )
}

# Put the plots together
p3 <- ggarrange(p2[[1]], p2[[2]], p2[[3]], ncol = 1, labels = c("(1)", "(2)", "(3)"))
p4 <- ggarrange(p2[[4]], p2[[5]], p2[[6]], ncol = 1, labels = c("(4)", "(5)", "(6)"))
p5 <- ggarrange(p3, p4, ncol = 2)

################################################################################
#### Store the Maps
################################################################################
# Store the plot
ggsave("04_Manuscript/99_Heatmaps.png"
  , plot   = p1
  , bg     = "white"
  , width  = 8
  , height = 4
  , scale  = 1
)
ggsave("04_Manuscript/99_HeatmapsIndividual.png"
  , plot   = p5
  , bg     = "white"
  , height = 7
  , width  = 8
  , scale = 1.4
)
