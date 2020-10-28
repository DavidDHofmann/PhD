
################################################################################
#### Plot Source Areas and Source Points
################################################################################
# Description: Plot of the source points / source areas

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)  # For data wrangling
library(lubridate)  # To handle dates nicely
library(ggpubr)     # For nice plots

# Load required data
areas1  <- readOGR("03_Data/03_Results/99_SourceAreas.shp")
areas2  <- readOGR("03_Data/03_Results/99_SourceAreas2.shp")
points  <- readOGR("03_Data/03_Results/99_SourcePoints2.shp")
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")

# Randomize source points
head(points)
pointsr <- createPoints(areas = areas2, points = points, n = 100, randomize = T)

# Plot of static source points
p1 <- tm_shape(areas1) +
    tm_borders(
        col = "gray70"
      , lwd = 3
    ) +
  tm_shape(areas2) +
    tm_polygons(
        col         = "gray40"
      , border.col  = "black"
    ) +
  tm_shape(points) +
    tm_dots(
        col = "orange"
      , size = 0.1
    ) +
  tm_shape(kaza) +
    tm_borders(
        col = "white"
      , lwd = 1
      , lty = 2
    ) +
  tm_layout(
      frame       = "black"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.col      = "black"
  ) +
  tm_scale_bar(
        position    = "left"
      , text.size   = 0.5
      , width       = 0.125
      , text.color  = "white"
  ) +
  tm_compass(
      color.light = "white"
    , color.dark  = "white"
    , text.color  = "white"
  ) +
  tm_credits("(a)"
    , position  = c("left", "top")
    , size      = 1.5
    , col       = "white"
)

# Plot of random source points
p2 <- tm_shape(areas1) +
    tm_borders(
        col = "gray70"
      , lwd = 3
    ) +
  tm_shape(areas2) +
    tm_polygons(
        col         = "gray40"
      , border.col  = "black"
    ) +
  tm_shape(pointsr) +
    tm_dots(
        col = "orange"
      , size = 0.01
    ) +
  tm_shape(kaza) +
    tm_borders(
        col = "white"
      , lwd = 1
      , lty = 2
    ) +
  tm_layout(
      frame       = "black"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.col      = "black"
  ) +
  tm_scale_bar(
        position    = "left"
      , text.size   = 0.5
      , width       = 0.125
      , text.color  = "white"
  ) +
  tm_compass(
      color.light = "white"
    , color.dark  = "white"
    , text.color  = "white"
  ) +
  tm_credits("(b)"
    , position  = c("left", "top")
    , size      = 1.5
    , col       = "white"
)

# Put plots together
p <- tmap_arrange(p1, p2, ncol = 2)

# Store the plot
CairoPDF("05_Manuscript2/99_SourcePoints.pdf", width = 12, height = 5.25)
p
dev.off()
