################################################################################
#### Overview of Original Camera Trap Locations
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(readxl)
library(raster)
library(terra)
library(rgdal)
library(sf)
library(tidyverse)
library(ggspatial)
library(tmaptools)
library(ggnewscale)
library(rworldmap)
library(geodata)

# Load deployments
dat <- read_xlsx("/home/david/ownCloud/University/15. PhD/General/Cameratrapping/01_General/01_Deployments.xlsx")

# Here, we only want to consider the original deployments
dat <- dat[1:57, ]

# What area do they span in total?
pts <- vect(dat, geom = c("Longitude", "Latitude"), crs = "epsg:4326")
buf <- buffer(pts, 7500)
buf <- aggregate(buf)
buf <- fillHoles(buf)
area <- expanse(buf, unit = "km")

# Get some numbers on the number of cameras
nrow(dat)
count(dat, Habitat)
mean(dat$Height)
sd(dat$Height)
mean(as.numeric(dat$DistanceToRoad))
sd(as.numeric(dat$DistanceToRoad))

# Some more data for visual aid
study <- readOGR("/home/david/ownCloud/University/15. PhD/General/GoogleEarth/Botswana.kmz", layer = "Studyarea")
roads <- readOGR("/home/david/ownCloud/University/15. PhD/General/GoogleEarth/Botswana.kmz", layer = "Roads_Gabs")

# Load water
water <- vect("/home/david/ownCloud/University/15. PhD/Chapter_2/03_Data/02_CleanData/MajorWaters.gpkg")

# Define extent
ext <- extent(study)

# Crop stuff
roads <- crop(roads, ext)
study <- crop(study, ext)

# Download satellite image
sat <- read_osm(ext + 0.5, type = "bing", zoom = 11)
sat <- as(sat, "Raster")
sat <- projectRaster(sat, crs = CRS("+init=epsg:4326"), method = "ngb")
sat <- as.data.frame(sat, xy = T)
sat <- na.omit(sat)

# Location of dog camp
camp <- matrix(c(23.63718, -19.52126), ncol = 2)
camp <- SpatialPoints(camp, proj4string = CRS("+init=epsg:4326"))

# What area do they span in total?
cams <- vect(dat, geom = c("Longitude", "Latitude"), crs = "epsg:4326")
buff <- buffer(cams, 7500)
buff <- aggregate(buff)
buff <- fillHoles(buff)
area <- expanse(buff, unit = "km")
area <- round(area, -3)

################################################################################
#### Larger Study Area
################################################################################
# Get shape of africa
data(countriesCoarseLessIslands)
africa <- vect(countriesCoarseLessIslands)
africa <- buffer(africa, width = 0)
africa <- africa[africa$GEO3major == "Africa", ]

# # Get a digital elevation model (optional)
# ele <- elevation_global(path = tempdir(), res = 5)
# ele <- crop(ele, africa)
# ele <- mask(ele, africa)
# 
# # Compute hillshade
# ele_amplified <- ele ** 1.6
# slope  <- terrain(ele_amplified, v = "slope", unit = "radians")
# aspect <- terrain(ele_amplified, v = "aspect", unit = "radians")
# hill   <- shade(slope, aspect)

# Convert to sf
africa_sf <- st_as_sf(africa)
africa_sf <- st_set_crs(africa_sf, 4326)
ext_sf    <- st_as_sf(as(ext, "SpatialPolygons"))
ext_sf    <- st_set_crs(ext_sf, 4326)

# Visualize
p1 <- ggplot() +
#   geom_raster(data = as.data.frame(hill, xy = T), aes(x = x, y = y, fill = hillshade)) +
#   geom_raster(data = sat, aes(x = x, y = y, fill = rgb(red, green, blue, maxColorValue = 255))) +
  scale_fill_identity() +
  geom_sf(data = africa_sf, fill = "transparent") +
  geom_sf(data = st_as_sf(water)[1, ], fill = "cornflowerblue", color = "transparent", alpha = 0.5) +
  geom_sf(data = ext_sf, color = "orange", fill = "orange", alpha = 0.5) +
  geom_sf_text(data = africa_sf, aes(label = `ADMIN.1`)) +
#   scale_fill_gradient(low = "gray80", high = "white") +
  theme_void() +
  coord_sf(
      xlim = (ext + 5)[c(1, 2)] + 1
    , ylim = (ext + 5)[c(3, 4)]
  )

################################################################################
#### Plot of Cameras
################################################################################
# Prepare plot of the cameras
p2 <- ggplot() +
  geom_raster(data = sat, aes(x = x, y = y, fill = rgb(red, green, blue, maxColorValue = 255))) +
  scale_fill_identity() +
  new_scale_fill() +
  geom_sf(data = st_as_sf(roads), col = "gray80", linewidth = 0.1) +
  geom_sf(data = st_as_sf(study), fill = NA, col = "gray80", linewidth = 0.3) +
  geom_sf(data = st_as_sf(cams), aes(fill = as.factor(cams$Habitat)), pch = 21) +
  geom_sf(data = st_as_sf(camp), col = "red") +
  geom_sf_text(data = st_as_sf(cams), aes(label = cams$CameraID), col = "white", size = 2.5, nudge_y = 0.02) +
  scale_fill_manual(name = "Habitat", values = c("#ffff00", "#ffa500", "#00ff00")) +
  xlab("Longitude") +
  ylab("Latitude") +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
    , bar_cols   = c("black", "white")
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("white", "white")
        , line_col  = NA
        , text_col  = "white"
        , text_size = 12
      )
  ) +
  theme(
      legend.position   = "none"
    , plot.background   = element_rect(fill = "transparent", color = NA)
    , legend.key        = element_rect(fill = NA)
    , legend.background = element_blank()
  ) +
  coord_sf(
      xlim = ext[1:2]
    , ylim = ext[3:4]
  )

# Store the plots
ggsave("/home/david/Schreibtisch/Cameratraps_01"
  , plot   = p1
  , width  = 5
  , height = 6
  , scale  = 1
  , bg     = "white"
  , device = png
)

ggsave("/home/david/Schreibtisch/Cameratraps_02"
  , plot   = p2
  , width  = 12
  , height = 14
  , scale  = 0.5
  , bg     = "white"
  , device = png
)
