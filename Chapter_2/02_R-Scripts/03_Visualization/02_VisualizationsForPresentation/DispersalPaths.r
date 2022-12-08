################################################################################
#### Dispersal Paths
################################################################################
# Description: Visualization of all dispersal trajectories collected so far

# Clear R's brain
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(tidyverse)      # To wrangle data
library(ggpubr)         # To arrange multiple plots
library(tmaptools)      # To download satellite image
library(sf)             # To handle spatial data
library(rgdal)          # To handle spatial data
library(raster)         # To handle spatial data
library(ggspatial)      # For scale bar and north arrow
library(ggdark)         # To get access to dark themes

# Load dispersal data
dat <- read_csv("03_Data/02_CleanData/00_General_Dispersers.csv")

# Keep only tracks of dispersers
dat <- subset(dat, State == "Disperser")

# Let's download a satellite image for the extent of interest
ext <- extent(21, 28, -21, -17.5)
ext <- as(ext, "SpatialPolygons")
crs(ext) <- "+init=epsg:4326"
sat <- read_osm(ext, zoom = 9, type = "bing")

# Convert to raster and reproject
sat <- as(sat, "Raster")
sat <- projectRaster(sat, crs = crs(ext), method = "ngb")
sat <- trim(sat)

# Convert raster to dataframe
sat <- as.data.frame(sat, xy = T)

# Load protected areas
prot <- readOGR("03_Data/02_CleanData/02_LandUse_ProtectedAreas.shp")
prot <- subset(prot, Desig == "National Park")
prot <- st_as_sf(prot)

# Load country borders
afri <- readOGR("03_Data/02_CleanData/00_General_Africa.shp")
afri <- crop(afri, ext)
afri <- st_as_sf(afri)

# Prepare plot
p1 <- ggplot(sat, aes(x = x, y = y, fill = rgb(red, green, blue, maxColorValue = 255))) +
  geom_raster() +
  geom_point(data = dat, inherit.aes = F, aes(x = x, y = y, col = DogName), size = 0.4) +
  geom_path(data = dat, inherit.aes = F, aes(x = x, y = y, col = DogName), size = 0.3) +
  scale_fill_identity() +
  coord_sf(crs = "+init=epsg:4326") +
  theme_minimal() +
  scale_color_viridis_d(direction = -1, begin = 0.5) +
  theme(legend.position = "none") +
  xlab("Longitude") +
  ylab("Latitude") +
  annotation_scale(
      location   = "br"
    , width_hint = 0.2
    , line_width = 1
    , height     = unit(0.15, "cm")
    # , bar_cols   = c("white", "white")
    , text_col   = "white"
    , pad_x    = unit(2, "cm")
    , pad_y    = unit(1.2, "cm")
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
    , pad_x    = unit(0.8, "cm")
    , pad_y    = unit(1.0, "cm")
  )

# Plot of dispersal paths only
p2 <- ggplot(dat, aes(x = x, y = y, col = DogName)) +
  geom_point(size = 1) +
  geom_path(size = 0.5) +
  coord_sf(crs = "+init=epsg:4326") +
  theme_minimal() +
  scale_color_viridis_d(direction = -1, begin = 0.5) +
  dark_theme_classic() +
  theme(
      legend.position       = "none"
    , panel.border          = element_blank()
    , axis.line             = element_blank()
    , panel.background      = element_rect(fill   = "transparent")
    , plot.background       = element_rect(fill   = "transparent", color = NA)
    , panel.grid.major      = element_blank()
    , panel.grid.minor      = element_blank()
    , legend.background     = element_rect(fill   = "transparent")
    , legend.box.background = element_rect(fill   = "transparent", color = "transparent")
  ) +
  xlab("Longitude") +
  ylab("Latitude")

# Prepare plot with protected areas
p3 <- ggplot(sat, aes(x = x, y = y, fill = rgb(red, green, blue, maxColorValue = 255))) +
  geom_raster() +
  geom_point(data = dat, inherit.aes = F, aes(x = x, y = y, col = DogName), size = 0.4) +
  geom_path(data = dat, inherit.aes = F, aes(x = x, y = y, col = DogName), size = 0.3) +
  # geom_sf(data = subset(prot, Name %in% c("Hwange", "Moremi")), fill = "white", col = "white", alpha = 0.2, inherit.aes = F) +
  geom_sf(data = afri, fill = NA, col = "white", inherit.aes = F, lty = 1) +
  scale_fill_identity() +
  coord_sf(crs = "+init=epsg:4326") +
  theme_minimal() +
  scale_color_viridis_d(direction = -1, begin = 0.5) +
  theme(legend.position = "none") +
  xlab("Longitude") +
  ylab("Latitude") +
  annotation_scale(
      location   = "br"
    , width_hint = 0.2
    , line_width = 1
    , height     = unit(0.15, "cm")
    # , bar_cols   = c("white", "white")
    , text_col   = "white"
    , pad_x    = unit(2, "cm")
    , pad_y    = unit(1.2, "cm")
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
    , pad_x    = unit(0.8, "cm")
    , pad_y    = unit(1.0, "cm")
  )

# Store the plot
ggsave(plot = p1, "DispersalPaths1.png", width = 16, height = 9, scale = 0.6)
ggsave(plot = p2, "DispersalPaths2.png", width = 16, height = 9, scale = 0.6, bg = "transparent")
ggsave(plot = p3, "DispersalPaths3.png", width = 16, height = 9, scale = 0.6, bg = "transparent")
