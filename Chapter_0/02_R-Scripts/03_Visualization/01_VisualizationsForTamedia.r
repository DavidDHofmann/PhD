############################################################
#### Preparation of Plot for Tamedia
############################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Set a seed
set.seed(12345)

# Load required packages
library(rgdal)
library(raster)
library(tmap)
library(tidyverse)
library(Cairo)
library(viridis)
library(davidoff)
library(rgeos)
library(RColorBrewer)
library(grid)
library(gridExtra)

############################################################
#### Plot of Africa
############################################################
# Load required data
africa <- "03_Data/02_CleanData/00_General_Africa.shp" %>% readOGR()
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>% readOGR()

# Remove small islands for plotting
africa <- subset(africa, !(ID %in% c(27:41, 689, 690)))

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Prepare a map of Africa
p1 <- tm_shape(africa) +
    tm_polygons(
        col = "gray90"
      , border.col = "gray70"
      , lwd = 0.7
    ) +
  tm_shape(kaza) +
    tm_polygons(
        col        = "gray30"
      , col.border = "gray30"
      , lwd        = 1.5
      , alpha      = 0.5
    ) +
  tm_shape(kaza_ext) +
    tm_borders(
        col = "black"
      , lwd = 2
    ) +
  tm_layout(
      frame       = F
    , bg.color    = "transparent"
    , legend.show = FALSE
)

# Look at the plot
p1

############################################################
#### Plot of LCCs
############################################################
# Load required data
perm        <- raster("03_Data/03_Results/99_PermeabilityMap.tif")
corrs       <- raster("03_Data/03_Results/99_LeastCostCorridors2.tif")
kaza        <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp") %>% as(., "SpatialLines")
africa      <- readOGR("03_Data/02_CleanData/00_General_Africa.shp") %>% as(., "SpatialLines")
africa_crop <- readOGR("03_Data/02_CleanData/00_General_Africa.shp") %>% crop(., kaza)
prot        <- readOGR("03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(1Class).shp")
nati        <- readOGR("03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(3Classes).shp")

# Subset to national parks
nati <- subset(nati, Desig == "National Park")

# Subset to national parks that we want to plot
nati <- subset(nati, Name %in% c("Mavinga", "Luengue-Luiana", "Kafue"
  , "Hwange", "Central Kalahari", "Chobe", "Moremi", "Matusadona", "Khaudum"))

# There is a double entry for Kafue, get rid of the erronous one
nati$Area <- gArea(nati, byid = TRUE)
nati <- subset(nati, Area != min(Area))

# Create a separate shapefile for the text. We have to change some of the
# coordinates to make sure that they don't overlap
nati_text <- nati
nati_text$x <- coordinates(nati_text)[, 1]
nati_text$y <- coordinates(nati_text)[, 2]
nati_text <- nati_text@data
nati_text$y[nati_text$Name == "Kafue"] <-
  nati_text$y[nati_text$Name == "Kafue"] + 0.5
nati_text$y[nati_text$Name == "Chobe"] <-
  nati_text$y[nati_text$Name == "Chobe"] - 0.1
nati_text$y[nati_text$Name == "Matusadona"] <-
  nati_text$y[nati_text$Name == "Matusadona"] - 0.1
coordinates(nati_text) <- c("x", "y")
crs(nati_text) <- CRS("+init=epsg:4326")

# Check how they align
plot(nati)
points(nati_text)

# Add "NP" to the text (on a new line)
head(nati_text)
nati_text$Name <- paste0(nati_text$Name, "\nNP")

# We only keep the countries of interest in the cropped africa file
africa_crop <- subset(africa_crop, COUNTRY %in% c(
    "Angola"
  , "Namibia"
  , "Botswana"
  , "Zimbabwe"
  , "Zambia")
)

# Rename countries to german
# africa_crop$COUNTRY <- as.character(c("Angola", "Botswana", "Namibia", "Sambia", "Simbabwe"))

# Rescale between 0 and 1
corrs <- calc(corrs, fun = function(x){
  (x - min(x)) / (max(x) - min(x))
})

# Prepare the plot for LCCs
p2 <- tm_shape(corrs) +
    tm_raster(
        palette         = "-Spectral"
      , style           = "cont"
      # , title           = "Wildhund-Korridore"
      , title           = "Wild Dog Corridors"
      , breaks          = c(0, 0.5, 1)
      , labels          = c("unsuitable", "", "suitable")
      # , labels          = c("ungeeignet", "", "geeignet")
      , legend.reverse  = T
    ) +
  tm_shape(nati) +
    tm_borders(
        col   = "black"
      , alpha = 0.6
    ) +
  tm_shape(kaza) +
    tm_lines(
        col = "black"
      , lwd = 2
    ) +
  tm_shape(nati_text) +
    tm_text("Name"
      , col       = "black"
      , alpha     = 0.6
      , fontface  = 3
      , size      = 0.5
      , shadow    = TRUE
    ) +
  tm_shape(africa) +
    tm_lines(
        col = "black"
      , lwd = 1
      , lty = 2
    ) +
  tm_shape(africa_crop) +
    tm_text("COUNTRY"
      , col       = "black"
      , just      = "bottom"
      , fontface  = 2
    ) +
  tm_layout(
      legend.outside          = T
    , legend.outside.position = c("left", "top")
    , frame.lwd               = 2
  ) +
  tm_scale_bar(
      position    = "left"
    , text.size   = 0.5
    , text.color  = "white"
    , width       = 0.125
  ) +
  tm_compass(
      text.color   = "white"
    , color.dark   = "white"
    , color.light  = "white"
)

# Add common legend
p3 <- p2 + tm_add_legend(
      type    = "line"
    , labels  = c("Kavango-Zambezi", "Conservation Area", "Country Borders")
    , labels  = c("Kavango-Zambezi", "Schutzgebiet", "Landesgrenzen")
    , col     = c("black", "white", "black")
    , lty     = c(1, 1, 2)
    , lwd     = c(2, 2, 1.5)
  ) + tm_layout(
      legend.title.fontface = 2
    , legend.title.size     = 1.2
    , legend.text.size      = 1
    , legend.height         = 0.5
)

# Store the plot
tmap_save(p3
  , "Plot.pdf"
  , insets_tm = p1
  , insets_vp = viewport(0.164, 0.3, width = 0.56, height = 0.56)
  , width     = 9
  , height    = 5.25
  , scale     = 1.1
)
