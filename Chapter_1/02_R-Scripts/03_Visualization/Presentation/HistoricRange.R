################################################################################
#### Plot of Historic Range
################################################################################
# Description: A cartoon like plot of the historic range of african wild dogs

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
library(sf)         # To plot stuff with sf

################################################################################
#### Data Preparation
################################################################################
# Load map of africa
africa    <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
africa2   <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
historic  <- readOGR("03_Data/01_RawData/DAVID/HistoricRange.shp")
dogs      <- readOGR("03_Data/02_CleanData/00_General_WildDogs_IUCN.shp")
kaza      <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")

# Buffer africa slightly
africa <- gBuffer(africa, width = 0.01)

# Disaggregate
africa <- disaggregate(africa)

# Identify area of each polygon
africa$Size <- gArea(africa, byid = T)

# Keep only the largest two
africa <- subset(africa, Size %in% sort(africa$Size, decreasing = T)[1:2])

# Simplify and smoothen africa shape
africa <- gSimplify(africa, tol = 0.5)
africa <- smooth(africa, method = "ksmooth")

# Smoothen and crop historic range
historic <- smooth(historic, method = "ksmooth")
historic <- gIntersection(historic, africa, byid = F)

# Simplify and smoothen dog distribution
dogs <- gSimplify(dogs, tol = 0.2)
dogs <- smooth(dogs, method = "ksmooth")

# Create a buffered polygon of africa
africa2 <- gBuffer(africa, width = 100 / 111000)

# Simplify and smoothen kaza
kaza <- gSimplify(kaza, tol = 0.1)
kaza <- smooth(kaza, method = "ksmooth")

# Identify wild dog strongholds
strong <- rbind(disaggregate(dogs)[c(1, 8, 2), ])

################################################################################
#### Plot
################################################################################
# Prepare a map of Africa (with KAZA)
p1 <- tm_shape(africa2) +
    tm_polygons(col = "gray30", border.col = "gray50", lwd = 5) +
  tm_shape(africa) +
    tm_polygons(
        col = "gray10"
      , lwd = 0.7
      , border.col = "gray30"
    ) +
  tm_shape(historic) +
    tm_polygons(
        col = "orange"
      , lwd = 0.7
      , border.col = "black"
      , alpha = 0.2
    ) +
  tm_shape(dogs) +
    tm_polygons(
        col           = "orange"
      , alpha         = 0.8
      , border.alpha  = 0
    ) +
  tm_shape(strong) +
    tm_polygons(
        col           = lighten("orange", 1.3)
      , border.alpha  = 0
    ) +
  # tm_shape(kaza) +
  #   tm_borders(
  #       col = "white"
  #     , lwd = 4
  #   ) +
  tm_layout(
      bg.color = "transparent"
    , frame = F
)
png("Historic.png", width = 1000, height = 1080, bg = "transparent", pointsize = 30)
p1
dev.off()

################################################################################
#### Using ggplot
################################################################################
# Keep the different shapes of the dog areas separated
dogs <- disaggregate(dogs)

# Convert everything to "SpatialDataFrames"
historic <- as(historic, "SpatialPolygonsDataFrame")
dogs     <- as(dogs, "SpatialPolygonsDataFrame")
strong   <- as(strong, "SpatialPolygonsDataFrame")

# Give the different shapes nice grouping names
historic$Group <- "historical"
dogs$Group     <- "present"
strong$Group   <- "stronghold"

# Put all into a single object
range <- rbind(historic, dogs, strong)

# Let's also add the image of a wild dog
dog <- image_read_svg("/home/david/ownCloud/University/15. PhD/General/Images/WildDog_Standing.svg")

# Plot
ggplot() +
  geom_sf(data = st_as_sf(africa2), col = "gray50", fill = "gray10", lwd = 1.5) +
  geom_sf(data = st_as_sf(range), aes(fill = Group, alpha = Group), lwd = 0) +
  scale_fill_manual(values = c("orange", "orange", lighten("orange", 1.3)), name = "Wild Dog Range") +
  scale_alpha_manual(values = c(0.2, 0.8, 1.0), name = "Wild Dog Range") +
  theme_void() +
  theme(legend.position = c(0.3, 0.3)) +
  annotation_raster(dog, xmin = -5, xmax = 5, ymin = -8, ymax = 2)

# Plot
ggplot() +
  geom_sf(data = st_as_sf(africa2), col = "gray40", fill = "gray93", lwd = 1.5) +
  geom_sf(data = st_as_sf(range), aes(alpha = Group), lwd = 0, fill = "sienna1") +
  scale_alpha_manual(values = c(0.1, 0.5, 1.0), name = "Wild Dog Range") +
  theme_void() +
  theme(legend.position = c(0.25, 0.2)) +
  annotation_raster(dog, xmin = -10, xmax = 5, ymin = -12, ymax = 2)

# Store
ggsave("HistoricRange.png", plot = last_plot(), scale = 0.75)
