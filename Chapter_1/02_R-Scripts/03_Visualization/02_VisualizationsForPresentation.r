############################################################
#### Preparation of all Plots for the Thesis
############################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)
library(raster)
library(rgeos)
library(smoothr)
library(tmap)
library(Cairo)
library(igraph)
library(rgdal)
library(davidoff)

############################################################
#### Random Network
############################################################
# Create random network
net <- erdos.renyi.game(n = 12, 20, type = "gnm")
plot(net, vertex.label = NA)

############################################################
#### Plot of Historic Range
############################################################
# Load map of africa
africa    <- shapefile("03_Data/02_CleanData/00_General_Africa")
africa2   <- shapefile("03_Data/02_CleanData/00_General_Africa")
historic  <- shapefile("03_Data/01_RawData/DAVID/HistoricRange")
dogs      <- shapefile("03_Data/02_CleanData/00_General_WildDogs_IUCN")
kaza      <- shapefile("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")

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

# Smoothen and crop historic range
historic <- smooth(historic, method = "ksmooth")
historic <- gIntersection(historic, africa, byid = F)

# Simplify and smoothen dog distribution
dogs <- gSimplify(dogs, tol = 0.2)
dogs <- smooth(dogs, method = "ksmooth")

# Create a buffered polygon of africa
africa2 <- gBuffer(africa, width = 100/111000)

# Simplify and smoothen kaza
kaza <- gSimplify(kaza, tol = 0.1)
kaza <- smooth(kaza, method = "ksmooth")
plot(dogs)

# Identify wild dog strongholds
strong <- rbind(disaggregate(dogs)[c(7, 19, 29), ])

# Prepare a map of Africa
p1 <- tm_shape(africa2) +
    tm_polygons(col = "gray70", border.col = "gray70", lwd = 2) +
  tm_shape(africa) +
    tm_polygons(
        col = "black"
      , lwd = 0.7
      , border.col = "gray70"
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
  tm_layout(
      asp         = 0.8
    , frame       = "white"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.color    = "white"
)

# Prepare a map of Africa where the kaza is added
p2 <- tm_shape(africa) +
    tm_polygons(col = "gray70", border.col = "gray70", lwd = 2) +
  tm_shape(africa) +
    tm_polygons(
        col = "black"
      , lwd = 0.7
      , border.col = "gray70"
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
  tm_shape(kaza) +
    tm_borders(
        col = "white"
      , lwd = 3
    ) +
  tm_layout(
      asp         = 0.8
    , frame       = "white"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.color    = "white"
)

# Store the plot
CairoPDF("Test", width = 5.25, height = 6)
p2
dev.off()

############################################################
#### Plot of KAZA-TFCA (Large)
############################################################
# Load required data
africa2 <- "03_Data/02_CleanData/00_General_Africa.shp" %>% readOGR()
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>% readOGR()
prot <- "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks.shp" %>% readOGR()
water <- "03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters_GEOFABRIK.shp" %>% readOGR(.)

# Clip africa layer
africa <- gIntersection(africa2, africa, byid = T)

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Create labels for countries
labels_countries <- data.frame(
    x = c(16, 24, 11, 38, 34)
  , y = c(-9, -30, -26, -13, -23)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Create lines pointing towards these countries
l1 <- rbind(c(17, -10), c(19, -12)) %>% spLines()
l2 <- rbind(c(24, -29), c(24, -23)) %>% spLines()
l3 <- rbind(c(13, -25), c(17, -22)) %>% spLines()
l4 <- rbind(c(35, -13), c(30, -14)) %>% spLines()
l5 <- rbind(c(33, -22), c(30, -19)) %>% spLines()
lines_countries <- rbind(l1, l2, l3, l4, l5)
crs(lines_countries) <- crs(labels_countries)

# Prepare a plot of the KAZA. Create labels for countries first.
labels_countries2 <- data.frame(
    x = c(20.39, 23.94, 20.07, 25.69, 28.22)
  , y = c(-15.28, -19.94, -19.39, -15.22, -18.9)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries2) <- c("x", "y")
crs(labels_countries2) <- CRS("+init=epsg:4326")

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.8, 25.6)
  , y     = c(-19, -18.2, -17, -20.7)
  , Label = c(
    "Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans"
  )
)
coordinates(labels_waters) <- c("x", "y")
crs(labels_waters) <- CRS("+init=epsg:4326")

# Create labels for some national parks
labels_nationalparks <- data.frame(
    x = c(26.56, 28.61, 21.15, 25.87, 20.38, 23.58, 23.71, 24.51, 20.78)
  , y = c(-19.08, -17.05, -17.26, -14.66, -16.08, -21.4, -19.29, -18.65, -18.81)
  , Label = paste0(c(
      "Hwange", "Matusadona", "Luengue-Luiana", "Kafue", "Mavinga"
    , "Central Kalahari", "Moremi", "Chobe", "Khaudum"
  ), "\nNP")
)
coordinates(labels_nationalparks) <- c("x", "y")
crs(labels_nationalparks) <- CRS("+init=epsg:4326")

# Prepare a map of Africa
p1 <- tm_shape(africa) +
    tm_polygons(
        col = "black"
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
      , lwd = 1.5
  )

# Store plot
CairoPDF("Test")
p1
dev.off()

# Prepare a map of Kaza
p2 <- tm_shape(prot) +
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
  tm_layout(bg.color = "black")

CairoPDF("Test")
p2
dev.off()
