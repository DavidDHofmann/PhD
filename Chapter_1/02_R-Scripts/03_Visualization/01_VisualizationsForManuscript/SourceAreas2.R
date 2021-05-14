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
library(raster)       # Manipulate rasters
library(RColorBrewer) # For nice color palettes
library(viridis)      # For nice colors
library(sf)           # To plot spatial stuff with ggplot
library(ggspatial)    # To plot spatial stuff with ggplot

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

# Load simulated trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")

# Identify the first coordinate of each trajectory
first <- sims %>%
  dplyr::select("x", "y", "TrackID", "Area") %>%
  group_by(TrackID) %>%
  slice(1) %>%
  SpatialPointsDataFrame(
      coords      = cbind(.[["x"]], .[["y"]])
    , proj4string = CRS("+init=epsg:4326")
  )

# Replace area names
first$Area <- ifelse(first$Area == "Main", "Main Source Points", "Buffer Source Points")
first$Area <- factor(first$Area, levels = c("Main Source Points", "Buffer Source Points"))

# Remove the rest of the simulations
rm(sims)
gc()

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
main@data <- data.frame(ID = 1:nrow(main), Area = "Main Source Areas")
buffer@data <- data.frame(ID = nrow(main) + 1, Area = "Buffer Source Areas")
areas <- rbind(main, buffer)

# Coerce the reference raster to a dataframe
r <- as.data.frame(r, xy = T)

# Coerce spatial objects to class sf
areas   <- st_as_sf(areas)
first   <- st_as_sf(first)
kaza    <- st_as_sf(kaza)
africa  <- st_as_sf(africa)

# Plot source areas
p1 <- ggplot() +
  geom_sf(
      data        = areas
    , mapping     = aes(fill = Area)
    , col         = c(darken("#70ab70"), darken("cornflowerblue"))[areas$Area]
    , lty         = 1
    , lwd         = 0.2
    , show.legend = T
    , alpha       = 0.8
  ) +
  geom_sf(
      data        = kaza
    , col         = "black"
    , fill        = NA
    , lty         = 1
    , lwd         = 0.5
    , show.legend = F
  ) +
  geom_sf(
      data        = africa
    , col         = "black"
    , fill        = NA
    , lty         = 2
    , lwd         = 0.2
    , show.legend = F
  ) +
  scale_fill_manual(
    values = c("#70ab70", "cornflowerblue")
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    , title    = "Source Areas"
  ) +
  theme(
    , panel.background  = element_blank()
    , panel.border      = element_rect(colour = "black", fill = NA, size = 1)
    , legend.position   = c(0.23, 0.85)
    , legend.key        = element_rect(fill = "transparent", colour = "transparent")
    , legend.background = element_rect(fill = "transparent", colour = "transparent")
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 1
    , height     = unit(0.15, "cm")
    , bar_cols   = c("black", "white")
    , text_col   = "black"
  ) +
  annotation_north_arrow(
      location = "bl"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , pad_y = unit(0.5, "cm")
    , pad_x = unit(-0.02, "cm")
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 12
      )
)

# Plot source points
p2 <- ggplot() +
  geom_sf(
      data        = first
    , mapping     = aes(col = Area)
    , size        = 0.01
    , show.legend = T
    , alpha       = 0.8
  ) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  geom_sf(
      data        = kaza
    , col         = "black"
    , fill        = NA
    , lty         = 1
    , lwd         = 0.5
    , show.legend = F
  ) +
  geom_sf(
      data        = africa
    , col         = "black"
    , fill        = NA
    , lty         = 2
    , lwd         = 0.2
    , show.legend = F
  ) +
  scale_color_manual(
    values = c("#70ab70", "cornflowerblue")
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    , title    = "Source Points"
  ) +
  theme(
    , panel.background  = element_blank()
    , panel.border      = element_rect(colour = "black", fill = NA, size = 1)
    , legend.position   = c(0.23, 0.85)
    , legend.key        = element_rect(fill = "transparent", colour = "transparent")
    , legend.title      = element_blank()
    , legend.background = element_rect(fill = "transparent", colour = "transparent")
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 1
    , height     = unit(0.15, "cm")
    , bar_cols   = c("black", "white")
    , text_col   = "black"
  ) +
  annotation_north_arrow(
      location = "bl"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , pad_y = unit(0.5, "cm")
    , pad_x = unit(-0.02, "cm")
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 12
      )
)

# Put the plots together
p3 <- ggarrange(p1, p2, ncol = 2, labels = c("a", "b"))

# Store plot
ggsave("test.png", plot = p3, width = 6, height = 3, scale = 2)
