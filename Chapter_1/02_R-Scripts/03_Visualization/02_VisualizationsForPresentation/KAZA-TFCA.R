################################################################################
#### Plot of KAZA-TFCA
################################################################################
# Description: Multiple plots depicting the KAZA-TFCA

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
# wd <- "C:/Users/david/switchdrive/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)  # For data wrangling
library(raster)     # To handle spatial data
library(rgeos)      # To manipulate spatial data
library(rgdal)      # To read and write spatial data
library(tmap)       # To plot nice maps
library(smoothr)    # To smooth spatial objects
library(davidoff)   # Custom functions

################################################################################
#### Data Preparation
################################################################################
# Load required data
africa <- "03_Data/02_CleanData/00_General_Africa_ESRI.shp" %>%
  readOGR()
africa2 <- "03_Data/02_CleanData/00_General_Africa_ESRI.shp" %>%
  readOGR()
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>%
  readOGR()
prot <- "03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp" %>%
  readOGR()
water <- "03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters_GEOFABRIK.shp" %>%
  readOGR(.)
disp <- "03_Data/02_CleanData/00_General_Dispersers_POPECOL.shp" %>%
  readOGR(.)

# Buffer slightly
africa <- gBuffer(africa, width = 0.01)

# Disaggregate
africa <- disaggregate(africa)

# Identify sizes of areas
africa$Size <- gArea(africa, byid = T)

# Keep only the largest two
africa <- subset(africa, Size %in% sort(africa$Size, decreasing = T)[1:2])

# Simplify and smoothen africa shape
africa <- gSimplify(africa, tol = 0.5)
africa <- smooth(africa, method = "ksmooth")

# Clip africa layer
africa <- gIntersection(africa2, africa, byid = T)

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Copy the africa layer
africa_copy <- africa

################################################################################
#### Plot
################################################################################
# Prepare a map of Africa
p1 <- tm_shape(africa_copy) +
    tm_borders(
        col = "gray50"
      , lwd = 5
    ) +
  tm_shape(africa) +
    tm_polygons(
        col = "gray10"
      , lwd = 0.7
      , border.col = "gray30"
    ) +
  tm_shape(kaza) +
    tm_polygons(
        col = "orange"
      , alpha = 0.5
      , border.col = "orange"
    ) +
  tm_shape(kaza_ext) +
    tm_borders(
        col = "white"
      , lty = 3
      , lwd = 2
  ) +
  tm_layout(
      bg.color = "transparent"
    , frame = F
)

# Plot of kaza only
p2 <- tm_shape(prot) +
    tm_polygons(
      , col           = "black"
      , border.col    = "black"
      , lwd           = 0
      , legend.show   = F
    ) +
  tm_shape(kaza, is.master = T) +
    tm_polygons(
        col        = "orange"
      , alpha      = 0.5
      , border.col = "orange"
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray70"
    ) +
  tm_scale_bar(
      position   = "left"
    , text.size  = 0.5
    , width      = 0.125
    , text.color = "white"
  ) +
  tm_compass(
      color.light = "white"
    , color.dark  = "white"
    , text.color  = "white"
  ) +
  tm_layout(
    bg.color = "black"
)

# Plot of kaza and protected areas
p3 <- tm_shape(prot) +
    tm_polygons(
      , col           = "gray40"
      , border.col    = NA
      , lwd           = 0
      , legend.show   = F
    ) +
  tm_shape(kaza, is.master = T) +
    tm_polygons(
        col = "orange"
      , alpha = 0.5
      , border.col = "orange"
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray70"
    ) +
  tm_scale_bar(
      position   = "left"
    , text.size  = 0.5
    , width      = 0.125
    , text.color = "white"
  ) +
  tm_compass(
      color.light = "white"
    , color.dark  = "white"
    , text.color  = "white"
  ) +
  tm_layout(bg.color = "black")

# Plot of kaza and protected areas and Hwange + Moremi
p4 <- tm_shape(prot) +
    tm_polygons(
        col           = "gray40"
      , border.col    = NA
      , lwd           = 0
      , legend.show   = F
    ) +
  tm_shape(kaza, is.master = T) +
    tm_polygons(
        col = "orange"
      , alpha = 0.5
      , border.col = "orange"
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray70"
    ) +
  tm_shape(subset(prot, Name %in% c("Moremi", "Hwange", "Kafue"))) +
    tm_polygons(
        col        = "white"
      , alpha      = 0.5
      , border.col = "white"
      , border.lwd = 3
    ) +
  tm_scale_bar(
      position   = "left"
    , text.size  = 0.5
    , width      = 0.125
    , text.color = "white"
  ) +
  tm_compass(
      color.light = "white"
    , color.dark  = "white"
    , text.color  = "white"
  ) +
  tm_layout(bg.color = "black")

# Plot of kaza only
p5 <- tm_shape(africa_copy) +
    tm_borders(
        col = "gray50"
      , lwd = 5
    ) +
  tm_shape(africa) +
    tm_polygons(
        col = "gray10"
      , lwd = 0.7
      , border.col = "gray30"
    ) +
  tm_shape(kaza) +
    tm_polygons(
        col = "white"
      , alpha = 0.8
      , border.col = "white"
    ) +
  tm_shape(kaza_ext) +
    tm_borders(
        col = "white"
      , lty = 3
      , lwd = 2
  ) +
  tm_layout(
      bg.color = "transparent"
    , frame = F
)

# Plot of kaza and protected areas and Dispersers
p6 <- tm_shape(prot) +
    tm_polygons(
        col           = "gray40"
      , border.col    = NA
      , lwd           = 0
      , legend.show   = F
    ) +
  tm_shape(kaza, is.master = T) +
    tm_polygons(
        col = "orange"
      , alpha = 0.5
      , border.col = "orange"
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray70"
    ) +
  tm_shape(disp) +
    tm_lines(
        col   = "white"
      , alpha = 0.5
      , lwd   = 1.5
    ) +
  tm_scale_bar(
      position   = "left"
    , text.size  = 0.5
    , width      = 0.125
    , text.color = "white"
  ) +
  tm_compass(
      color.light = "white"
    , color.dark  = "white"
    , text.color  = "white"
  ) +
  tm_layout(bg.color = "black")


# Store the plots
png("Plot1.png", width = 1080, height = 720, bg = "transparent")
p1
dev.off()
png("Plot2.png", width = 1080, height = 720, bg = "transparent", pointsize = 30)
p2
dev.off()
png("Plot3.png", width = 1080, height = 720, bg = "transparent", pointsize = 30)
p3
dev.off()
png("Plot4.png", width = 1080, height = 720, bg = "transparent", pointsize = 30)
p4
dev.off()
png("Plot5.png", width = 1080, height = 720, bg = "transparent", pointsize = 30)
p5
dev.off()
png("Plot6.png", width = 1080, height = 720, bg = "transparent", pointsize = 30)
p6
dev.off()

# Prepare raster
map <- raster(ncol = 8, nrow = 8, xmn = 0, xmx = 8, ymn = 0, ymx = 8)
map2 <- raster(ncol = 7, nrow = 7, xmn = 0.5, xmx = 7.5, ymn = 0.5, ymx = 7.5)
values(map) <- 1
map[1, 1:5] <- 2
map[2, 1:2] <- 2
map[3, 1:3] <- 2
map[3, 5] <- 2
map[4, 1:2] <- 2
map[5, 1:2] <- 2
map[1, 6:7] <- 3
map[2, 6:7] <- 3
map[6, 1:3] <- 3
map[7, 1:3] <- 3
map[8, 1:3] <- 3
map[1, 8] <- 4
map[2, 8] <- 4
map[3, 8] <- 4
map[4, 5:6] <- 4
map[5, 4:6] <- 4
map[5, 8] <- 4
map[6, 4:5] <- 4
map[6, 7:8] <- 4
map[7, 7:8] <- 4
map[8, 4:8] <- 4

# Set colors
colors <- c(
    rgb(146, 209, 79, max = 255)
  , rgb(255, 210, 67, max = 255)
  , rgb(0, 128, 1, max = 255)
  , "gray"
)

# Visualize
png("raster0.png", width = 1080, height = 1080, bg = "transparent")
plot(map, col = colors, axes = F, box = F, legend = F)
dev.off()

png("raster1.png", width = 1080, height = 1080, bg = "transparent")
plot(map, col = colors, axes = F, box = F, legend = F)
plot(as(map, "SpatialPointsDataFrame"), add = T, pch = 20, cex = 5)
dev.off()

png("raster2.png", width = 1080, height = 1080, bg = "transparent")
plot(map, col = colors, axes = F, box = F, legend = F)
plot(as(map, "SpatialPointsDataFrame"), add = T, pch = 20, cex = 5)
plot(as(map2, "SpatialPolygonsDataFrame"), add = T, pch = 20, lwd = 5)
dev.off()
