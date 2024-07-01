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
library(viridis)      # For nice colors
library(sf)           # For plotting
library(ggdark)       # For dark themes
library(ggspatial)    # For scale bars and north arrows
library(ggpubr)       # To put plots together

################################################################################
#### Extended Rasters & Source Points
################################################################################
# Load required data
main   <- readOGR("03_Data/03_Results/99_SourceAreas.shp")
buffer <- readOGR("03_Data/03_Results/99_BufferArea.shp")
prot   <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
kaza   <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
africa <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")

# Calculate the extent of the entire study area (including the buffer)
ext1 <- as(extent(main), "SpatialPolygons")
ext2 <- as(extent(buffer), "SpatialPolygons")
crs(ext1) <- CRS("+init=epsg:4326")
crs(ext2) <- CRS("+init=epsg:4326")
gArea(spTransform(ext1, CRS("+init=epsg:32734"))) * 1e-6
gArea(spTransform(ext2, CRS("+init=epsg:32734"))) * 1e-6

# Load simulated trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")

# Identify the first coordinate of each trajectory
first <- sims %>%
  dplyr::select("x", "y", "TrackID", "Area") %>%
  group_by(TrackID) %>%
  slice(1)

# Remove the rest of the simulations
rm(sims)
gc()

# Prepare country labels
labels_countries <- data.frame(
    x     = c(20.39, 23.94, 20.07, 25.69, 28.22)
  , y     = c(-15.28, -19.94, -19.39, -15.22, -18.9)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Identify protected areas too small to be considered
ints <- gIntersects(main, prot, byid = T)
ints <- unname(rowSums(ints))
small <- prot[!ints, ]
small <- aggregate(small)
small <- as(small, "SpatialPolygonsDataFrame")

# Put buffer and main study area together
main@data <- data.frame(ID = 1:nrow(main), Area = "Main Source Areas")
buffer@data <- data.frame(ID = nrow(main) + 1, Area = "Buffer Source Areas")
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

# Make areas factorial
first$Area <- factor(first$Area, levels = c("Main", "Buffer"))

# Convert everything to sf
areas            <- st_as_sf(areas)
kaza             <- st_as_sf(kaza)
africa           <- st_as_sf(africa)
labels_countries <- st_as_sf(labels_countries)

# Plot of source patches
p1 <- ggplot() +
  geom_sf(
      data        = areas
    , aes(col = Area, fill = Area)
    , lty         = 1
    , lwd         = 0.1
    , show.legend = F
    , alpha       = 0.8
  ) +
  geom_sf_text(
      data     = labels_countries
    , mapping  = aes(label = Label)
    , col      = "white"
    , fontface = 2
    , size     = 5
  ) +
  geom_sf(
      data        = kaza
    , col         = "white"
    , fill        = NA
    , lty         = 1
    , lwd         = 1
    , show.legend = F
  ) +
  geom_sf(
      data        = africa
    , col         = "white"
    , fill        = NA
    , lty         = 2
    , lwd         = 0.5
    , show.legend = F
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(xmin(ext2), xmax(ext2))
    , ylim   = c(ymin(ext2), ymax(ext2))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    # , title    = "Betweenness"
    # , subtitle = "After 2'000 Steps"
  ) +
  dark_theme_classic() +
  scale_color_manual(values = c("orange", "cornflowerblue")) +
  scale_fill_manual(values = c("orange", "cornflowerblue")) +
  theme(
      legend.position       = "none"
    , panel.border          = element_rect(colour = "gray30", fill       = NA, size = 1)
    , axis.line             = element_line(colour = "gray30")
    , panel.background      = element_rect(fill   = "transparent")
    , plot.background       = element_rect(fill   = "transparent", color = NA)
    , panel.grid.major      = element_blank()
    , panel.grid.minor      = element_blank()
    , legend.background     = element_rect(fill   = "transparent")
    , legend.box.background = element_rect(fill   = "transparent", color = "transparent")
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 1
    , height     = unit(0.15, "cm")
    , bar_cols   = c("white", "white")
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
  )

# Plot of source points
p2 <- ggplot() +
  geom_point(
      data = first
    , aes(x = x, y = y, col = Area)
    , size = 0.1
    , alpha = 0.5
  ) +
  geom_sf_text(
      data     = labels_countries
    , mapping  = aes(label = Label)
    , col      = "white"
    , fontface = 2
    , size     = 5
  ) +
  geom_sf(
      data        = kaza
    , col         = "white"
    , fill        = NA
    , lty         = 1
    , lwd         = 1
    , show.legend = F
  ) +
  geom_sf(
      data        = africa
    , col         = "white"
    , fill        = NA
    , lty         = 2
    , lwd         = 0.5
    , show.legend = F
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(xmin(ext2), xmax(ext2))
    , ylim   = c(ymin(ext2), ymax(ext2))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    # , title    = "Betweenness"
    # , subtitle = "After 2'000 Steps"
  ) +
  dark_theme_classic() +
  scale_color_manual(values = c("orange", "cornflowerblue")) +
  theme(
      legend.position       = "none"
    , panel.border          = element_rect(colour = "gray30", fill       = NA, size = 1)
    , axis.line             = element_line(colour = "gray30")
    , panel.background      = element_rect(fill   = "transparent")
    , plot.background       = element_rect(fill   = "transparent", color = NA)
    , panel.grid.major      = element_blank()
    , panel.grid.minor      = element_blank()
    , legend.background     = element_rect(fill   = "transparent")
    , legend.box.background = element_rect(fill   = "transparent", color = "transparent")
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 1
    , height     = unit(0.15, "cm")
    , bar_cols   = c("white", "white")
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
  )


# Put plots together
p <- ggarrange(p1, p2)

# Store plot
ggsave(
    plot     = p
  , filename = "SourceAreas.png"
  , width    = 9
  , height   = 3.75
  , bg       = "transparent"
  , scale    = 1.5
)
