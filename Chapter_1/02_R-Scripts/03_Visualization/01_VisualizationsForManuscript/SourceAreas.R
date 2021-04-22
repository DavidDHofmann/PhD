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
library(viridis)      # For nice colors

################################################################################
#### Extended Rasters & Source Points
################################################################################
# Load required data
main   <- readOGR("03_Data/03_Results/99_SourceAreas.shp")
buffer <- readOGR("03_Data/03_Results/99_BufferArea.shp")
prot   <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
kaza   <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
africa <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
r      <- raster("03_Data/02_CleanData/00_General_Raster.tif")
water  <- raster("03_Data/02_CleanData/01_LandCover_WaterCoverAveraged_MERGED.tif")

# Prepare country labels
labels_countries <- data.frame(
    x = c(20.39, 23.94, 20.07, 25.69, 28.22)
  , y = c(-15.28, -19.94, -19.39, -15.22, -18.9)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Extend the reference raster
r <- extendRaster(r, extent(r) + c(-1, 1, -1, 1) * metersToDegrees(100000))

# And crop it to the buffer
r <- crop(r, extent(buffer))

# Identify protected areas too small to be considered
ints <- gIntersects(main, prot, byid = T)
ints <- unname(rowSums(ints))
small <- prot[!ints, ]
small <- aggregate(small)
small <- as(small, "SpatialPolygonsDataFrame")

# Put buffer and main study area together
main@data <- data.frame(ID = 1:nrow(main), Area = "Source Areas (Main)")
buffer@data <- data.frame(ID = nrow(main) + 1, Area = "Source Areas (Buffer)")
areas <- rbind(main, buffer)

# Add labels for main and buffer area
labels_areas <- data.frame(
    x      = c(18.1, 18.1)
  , y      = c(-12.3, -12.9)
  , Label  = c("Buffer Area (n = 30'000)", "Main Area (n = 50'000)")
  , Label2 = c("Buffer Area", "Main Area")
)
coordinates(labels_areas) <- c("x", "y")
crs(labels_areas) <- CRS("+init=epsg:4326")

# Expand water
water_ext <- extendRaster(water, r)

# Plot
p1 <- tm_shape(r) +
    tm_raster(
        palette     = "white"
      , legend.show = F
    ) +
  tm_shape(water_ext) +
    tm_raster(
        palette = c("white", "gray30")
      , title   = ""
      , labels  = c("", "Water")
      , style   = "cat"
    ) +
  tm_shape(buffer) +
    tm_polygons(
        col          = "cornflowerblue"
      , border.col   = "cornflowerblue"
      , alpha        = 0.6
      , lwd          = 0
    ) +
  # tm_shape(labels_areas) +
  #   tm_text("Label2"
  #     , col             = "Label2"
  #     , palette         = c("cornflowerblue", "black")
  #     , fontface        = 3
  #     , size            = 0.7
  #     , just            = "left"
  #     , legend.col.show = F
  #   ) +
  tm_shape(kaza) +
    tm_borders(
        col = "black"
      , lty = 1
      , lwd = 2
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray20"
      , lwd = 0.5
    ) +
  tm_shape(labels_countries) +
    tm_text("Label"
      , col       = "gray20"
      , fontface  = 2
      , size      = 0.6
    ) +
  tm_graticules(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_compass(
      color.dark  = "black"
    , color.light = "black"
    , text.color  = "black"
    , position    = c(0.9, 0.08)
  ) +
  tm_scale_bar(
        position  = c(0.8, 0)
      , text.size = 0.5
      , text.col  = "black"
      , width     = 0.125
  ) +
  tm_credits("a"
    , position = c(0.01, 0.9)
    , size     = 1.5
    , col      = "black"
    , fontface = "bold"
  ) +
  tm_layout(
    , frame            = "gray20"
    , frame.lwd        = 3
    , legend.position  = c(0.068, 0.082)
    , legend.stack     = "vertical"
    , legend.text.size = 0.8
    , legend.bg.color  = "white"
  )

# Plot
p2 <- tm_shape(r) +
    tm_raster(
        palette     = "white"
      , legend.show = F
    ) +
  tm_shape(areas) +
    tm_polygons(
        col         = "Area"
      , title       = ""
      , palette     = c("#70ab70", "cornflowerblue")
      , lwd         = 0
      , border.col  = "#6ba36b"
      , legend.show = T
      , alpha       = 0.6
    ) +
  tm_shape(kaza) +
    tm_borders(
        col = "black"
      , lty = 1
      , lwd = 2
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray20"
      , lwd = 0.5
    ) +
  tm_shape(labels_countries) +
    tm_text("Label"
      , col       = "gray20"
      , fontface  = 2
      , size      = 0.6
    ) +
  # tm_shape(labels_areas) +
  #   tm_text("Label"
  #     , col             = "Label"
  #     , palette         = c("cornflowerblue", "black")
  #     , fontface        = 3
  #     , size            = 0.7
  #     , just            = "left"
  #     , legend.col.show = F
  #   ) +
  tm_graticules(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
    , frame            = "gray20"
    , frame.lwd        = 3
    , legend.position  = c(0.067, 0.082)
    , legend.stack     = "vertical"
    , legend.text.size = 0.8
    , legend.bg.color  = "white"
  ) +
  tm_compass(
      color.dark  = "black"
    , color.light = "black"
    , text.color  = "black"
    , position    = c(0.9, 0.08)
  ) +
  tm_scale_bar(
        position  = c(0.8, 0)
      , text.size = 0.5
      , text.col  = "black"
      , width     = 0.125
  ) +
  tm_credits("b"
    , position = c(0.01, 0.9)
    , size     = 1.5
    , col      = "black"
    , fontface = "bold"
)

# Put plots together
p <- tmap_arrange(p1, p2, nrow = 1)

# Store plot
tmap_save(p
  , filename = "04_Manuscript/99_SourceAreas.png"
  , width    = 9
  , height   = 3.75
)
