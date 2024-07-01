################################################################################
#### Nice Plot of All Protected Areas in the World
################################################################################
# Description: Following this blog-post
# https://strimas.com/post/hexagonal-grids/
# https://stanford.edu/~vbauer/teaching/hillshade.html

# Clear R's brain
rm(list = ls())

# Load required packages
library(sf)
library(stars)
library(raster)
library(rgeos)
library(viridis)
library(rasterVis)
library(tidyverse)
library(velox)
library(rworldmap)
library(cleangeo)
library(rgdal)
library(devtools)
library(ggnewscale)
library(davidoff)

################################################################################
#### Preparation of Protected Areas
################################################################################
# Load protected areas
files <- c(
    "/home/david/Downloads/ProtectedAreas/Part1/WDPA_Aug2020-shapefile-polygons.shp"
  , "/home/david/Downloads/ProtectedAreas/Part2/WDPA_Aug2020-shapefile-polygons.shp"
  , "/home/david/Downloads/ProtectedAreas/Part3/WDPA_Aug2020-shapefile-polygons.shp"
)
prot1 <- read_sf(files[1])
prot2 <- read_sf(files[2])
prot3 <- read_sf(files[3])

# Put them together and remove non-terrestrial protected areas
prot <- rbind(prot1, prot2, prot3)
prot_land <- subset(prot, MARINE == 0)
prot_water <- subset(prot, MARINE != 0)

# Remove undesired stuff
rm(prot, prot1, prot2, prot3, files)
gc()

# Give a weight of one
prot_land$Weight <- 1
prot_water$Weight <- 1

# Rasterize protected areas onto a raster that covers the entire world
map <- getMap(resolution = "less islands")
map <- aggregate(map, dissolve = T)
r <- st_as_stars(raster(map, resolution = metersToDegrees(5000), vals = 0))
prot_land_rast <- st_rasterize(prot_land, r)
prot_water_rast <- st_rasterize(prot_water, r)

# Convert stars to raster
prot_land_rast <- as(prot_land_rast, "Raster")
prot_water_rast <- as(prot_water_rast, "Raster")

# Make raster binary
prot_land_rast <- prot_land_rast > 0
prot_water_rast <- prot_water_rast > 0

# Plot rasterized protected areas
plot(prot_land_rast)
plot(prot_water_rast)

# Store the rasters
writeRaster(prot_land_rast, "WDPA_Land.tif", overwrite = T)
writeRaster(prot_water_rast, "WDPA_Water.tif", overwrite = T)

################################################################################
#### Prepare Hexagons
################################################################################
# Reload them
prot_land_rast <- raster("WDPA_Land.tif")
prot_water_rast <- raster("WDPA_Water.tif")

# Get the map extent
extent <- as(extent(map), "SpatialPolygons")
crs(extent) <- crs(map)

# Create a mesh of hexagons
pts <- spsample(extent, type = "hexagonal", n = 10000)
hex <- HexPoints2SpatialPolygons(pts)
hex$ID <- 1:length(hex)

# # We only want to keep hexagons that cover land
# hex <- subset(hex, as.vector(gIntersects(map, hex, byid = T)))

# Identify how much of each hexagon is under protection
prot_land_velox <- velox(prot_land_rast)
prot_water_velox <- velox(prot_water_rast)
prot_land_ext <- prot_land_velox$extract(hex)
prot_water_ext <- prot_water_velox$extract(hex)
hex$ProtLand <- sapply(prot_land_ext, mean)
hex$ProtWater <- sapply(prot_water_ext, mean)
hex$ProtLand <- hex$ProtLand * 100
hex$ProtWater <- hex$ProtWater * 100

# Coerce to sf object
hex_terre <- subset(hex, as.vector(gIntersects(map, hex, byid = T)))
hex_land <- subset(hex, ProtLand > 0)
hex_water <- subset(hex, ProtWater > 0)

################################################################################
#### Plot
################################################################################
# Visualize Protected areas (Water + Land, 1 Color)
ggplot() +
  geom_sf(data = st_as_sf(hex_water), aes(fill = ProtWater), lwd = 0.1, col = "black") +
  scale_fill_gradient(low = "black", high = darken("darkolivegreen1"), trans = "sqrt") +
  geom_sf(data = st_as_sf(hex_terre), fill = "gray30", col = "gray29", lwd = 0.1) +
  new_scale_fill() +
  geom_sf(data = st_as_sf(hex_land), aes(fill = ProtLand), lwd = 0.1, col = "gray29") +
  scale_fill_gradient(low = "gray30", high = "darkolivegreen1", trans = "sqrt") +
  guides(fill = guide_colourbar(title.vjust = 2, title = "")) +
  theme(
      plot.background   = element_rect(fill = "transparent", color = NA)
    , panel.background  = element_rect(fill = "black", color = NA)
    , panel.grid.major  = element_blank()
    , legend.position   = "bottom"
    , legend.key.width  = unit(1,"cm")
    , legend.background = element_rect(fill = "transparent")
    , legend.text       = element_text(color = "white")
  )

# Store the plot
ggsave("ProtectedAreas0.png", scale = 1.4, device = "png", bg = "transparent")

# Visualize (Water + Land, 2 Colors)
ggplot() +
  geom_sf(data = st_as_sf(hex_water), aes(fill = ProtWater), lwd = 0.1, col = "black") +
  scale_fill_gradient(low = "black", high = darken("cornflowerblue"), trans = "sqrt") +
  geom_sf(data = st_as_sf(hex_terre), fill = "gray30", col = "gray29", lwd = 0.1) +
  new_scale_fill() +
  geom_sf(data = st_as_sf(hex_land), aes(fill = ProtLand), lwd = 0.1, col = "gray29") +
  scale_fill_gradient(low = "gray30", high = "orange", trans = "sqrt") +
  guides(fill = guide_colourbar(title.vjust = 2, title = "")) +
  theme(
      plot.background   = element_rect(fill = "transparent", color = NA)
    , panel.background  = element_rect(fill = "black", color = NA)
    , panel.grid.major  = element_blank()
    , legend.position   = "bottom"
    , legend.key.width  = unit(1,"cm")
    , legend.background = element_rect(fill = "transparent")
    , legend.text       = element_text(color = "white")
  )
ggsave("ProtectedAreas1.png", scale = 1.4, device = "png", bg = "transparent")

# Visualize (Land)
ggplot() +
  geom_sf(data = st_as_sf(hex_terre), fill = "gray30", col = "gray29", lwd = 0.1) +
  geom_sf(data = st_as_sf(hex_land), aes(fill = ProtLand), lwd = 0.1, col = "gray29") +
  scale_fill_gradient(low = "gray30", high = "orange", trans = "sqrt") +
  guides(fill = guide_colourbar(title.vjust = 2, title = "")) +
  theme(
      plot.background   = element_rect(fill = "transparent", color = NA)
    , panel.background  = element_rect(fill = "black", color = NA)
    , panel.grid.major  = element_blank()
    , legend.position   = "bottom"
    , legend.key.width  = unit(1,"cm")
    , legend.background = element_rect(fill = "transparent")
    , legend.text       = element_text(color = "white")
  )
ggsave("ProtectedAreas2.png", scale = 1.4, device = "png", bg = "transparent")
