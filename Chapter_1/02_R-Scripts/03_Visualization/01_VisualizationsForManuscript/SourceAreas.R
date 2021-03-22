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
library(tidyverse)    # For data wrangling
library(lubridate)    # To handle dates nicely
library(rgdal)        # To load spatial data
library(rgeos)        # To manipulate spatial data
library(davidoff)     # Custom functions
library(tmap)         # For nice maps
library(raster)       # Manipulate rasters
library(RColorBrewer) # For nice color palettes
library(Cairo)        # To store the plots

################################################################################
#### Source Points
################################################################################
# Load required data
source_areas  <- readOGR("03_Data/03_Results/99_SourceAreas.shp")
buffer_areas  <- readOGR("03_Data/03_Results/99_BufferArea.shp")
prot   <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
kaza   <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
r      <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Extend the reference raster
r <- extendRaster(r, extent(r) + c(-1, 1, -1, 1) * metersToDegrees(100000))

# Identify protected areas outside our source areas
ints <- gIntersects(source_areas, prot, byid = T)
ints <- unname(rowSums(ints))
small_areas <- prot[!ints, ]
small_areas <- aggregate(small_areas)
small_areas <- as(small_areas, "SpatialPolygonsDataFrame")

# Create random source points
n_buffer  <- 2000
n_areas   <- 5000
points1 <- spsample(source_areas, n = n_areas, type = "random")
points2 <- spsample(buffer_areas, n = n_buffer, type = "random")
points <- rbind(points1, points2)

# Put buffer and main study area together
source_areas@data <- data.frame(ID = 1:nrow(source_areas), Area = "Main")
buffer_areas@data <- data.frame(ID = nrow(source_areas) + 1, Area = "Buffer")
small_areas@data  <- data.frame(ID = nrow(source_areas) + 2, Area = "Small")
areas <- rbind(buffer_areas, source_areas, small_areas)

# Visualize it
p1 <- tm_shape(r) +
    tm_raster(palette = "white", legend.show = F) +
  tm_shape(areas) +
    tm_polygons(col = "Area", palette = c("gray20", "gray40", "gray80"), lwd = 0) +
  tm_shape(kaza) +
    tm_borders(col = "black", lty = 2, lwd = 2) +
  tm_shape(points) +
    tm_dots(col = "orange", size = 0.01, alpha = 0.7, border.lwd = 0) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
    , frame                   = "gray20"
    , frame.lwd               = 3
  ) +
  tm_scale_bar(
        position  = c("right", "bottom")
      , text.size = 0.5
      , text.col  = "white"
      , width     = 0.125
  ) +
  tm_credits("a"
    , position = c("left", "top")
    , size     = 1.5
    , col      = "white"
    , fontface = "bold"
  ) +
  tm_compass(
      color.dark  = "white"
    , color.light = "white"
    , text.color  = "white"
    , position    = c("left", "bottom")
  ) +
  tm_layout(
      legend.position = c(0.08, 0.78)
    , legend.bg.color = "white"
)

# Store the plot
CairoPDF("04_Manuscript/99_SourcePoints.pdf", width = 7, height = 6.25)
p1
dev.off()

# Store plot as png too
png("04_Manuscript/99_SourcePoints.png", width = 1080, heigh = 1080, pointsize = 28)
p1
dev.off()
