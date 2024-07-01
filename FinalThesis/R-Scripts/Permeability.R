################################################################################
#### Permeability Example
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/FinalThesis"
setwd(wd)

# Load required packages
library(terra)
library(tidyverse)
library(sf)
library(ggnewscale)

# Function to rotate data
rotateSF <- function(data, x_add = 0, y_add = 0, value = 0) {

  # Convert object to type sf
  if (inherits(data, "SpatVector")) {
      data <- st_as_sf(data)
    } else if (inherits(data, "SpatRaster")) {
      data <- st_as_stars(data)
      data <- st_as_sf(data)
  }

  # Create necessary matrices
  shear_matrix <- function() {
    rbind(
        c(2, 0)
      , c(1, 0.5)
    )
  }
  rotate_matrix <- function(value) {
    rbind(
        c(cos(value), -sin(value))
      , c(sin(value), cos(value))
    )
  }

  # Run transformation
  data <- mutate(data, geometry = geometry *
    shear_matrix() *
    rotate_matrix(value) +
    c(x_add, y_add)
  )

  # Return
  return(data)
}

# Recreate permeability surface from Diniz et al. 2019
perm <- rast(ncol = 8, nrow = 8, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
perm[]                <- 1
perm[4, 4:8]          <- 0
perm[5, 4:8]          <- 0
perm[6, 4:8]          <- 0
perm[7, 4:8]          <- 0
perm[8, 4:8]          <- 0
perm[1, 6:8]          <- 3
perm[2, 6:8]          <- 3
perm[3, 6:8]          <- 3
perm[6, 1:3]          <- 3
perm[7, 1:3]          <- 3
perm[8, 1:3]          <- 3
perm[11:13]           <- 2
perm[20]              <- 2
perm[c(27:28, 31:32)] <- 2
perm[c(35, 39)]       <- 2
perm[46]              <- 2
perm[52:54]           <- 2
names(perm) <- "Permeability"

# Split by values
perm_seg <- segregate(perm)

# Convert to polygon
perm_sf <- as.polygons(perm, aggregate = F)
perm_sf <- st_as_sf(perm_sf)

# Also derive a frame so we can smooth out the edges of the layers in the plot
frame <- as.polygons(ext(perm), crs = crs(perm))
frame <- st_as_sf(frame)

# Specify parmaters
framecolor <- "gray50"
spacing    <- 0.25

# Base plot
p0 <- ggplot()

# Add single layers on top of each other
for (i in 1:nlyr(perm_seg)) {

  # Prepare respective layer
  col             <- terrain.colors(4, rev = T)[i]
  perm_lyr        <- as.polygons(perm_seg[[i]], aggregate = F)
  names(perm_lyr) <- "Value"
  perm_lyr        <- st_as_sf(perm_lyr)

  # Add to stack
  p0 <- p0 +
    geom_sf(
        data        = rotateSF(perm_lyr, y_add = spacing * i)
      , mapping     = aes(fill = Value)
      , col         = "gray50"
      , show.legend = F
    ) +
    geom_sf(
        data        = rotateSF(frame, y_add = spacing * i)
      , fill        = "transparent"
      , col         = framecolor
      , linewidth   = 1
      , show.legend = F
    ) +
    scale_fill_gradient(low = "white", high = col) +
    new_scale_fill() +
    theme(legend.position = "none")
}

# Add resulting permeability map on top
p1 <- p0 +
  geom_sf(
      data    = rotateSF(perm_sf, y_add = 2)
    , mapping = aes(fill = Permeability)
    , color   = "black"
    , size    = 1
  ) +
  geom_sf(
      data        = rotateSF(frame, y_add = 2)
    , fill        = "transparent"
    , col         = framecolor
    , linewidth   = 1
    , show.legend = F
  ) +
  theme_void() +
  scale_fill_gradientn(colours = terrain.colors(100, rev = T)) +
  guides(
    fill = guide_colorbar(
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = F
      , barheight      = unit(0.2, "cm")
      , barwidth       = unit(6, "cm")
      , label          = F
    )
  ) +
  theme(legend.position = "top")

# Also a plot of the surface alone
p2 <- ggplot() +
  geom_sf(
      data    = perm_sf
    , mapping = aes(fill = Permeability)
    , color   = "black"
    , size    = 1
  ) +
  geom_sf(
      data        = frame
    , fill        = "transparent"
    , col         = framecolor
    , linewidth   = 1
    , show.legend = F
  ) +
  theme_void() +
  scale_fill_gradientn(colours = terrain.colors(100, rev = T)) +
  guides(
    fill = guide_colorbar(
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = F
      , barheight      = unit(0.2, "cm")
      , barwidth       = unit(6, "cm")
      , label          = F
    )
  ) +
  theme(legend.position = "none")

# Store it
ggsave(
    plot     = p1
  , filename = "Chapter_99/Figures/PermeabilityStack.png"
  , bg       = "transparent"
  , height   = 4
  , width    = 4.5
  , device   = png
  , scale    = 1
)

ggsave(
    plot     = p2
  , filename = "Chapter_99/Figures/PermeabilitySurface.png"
  , bg       = "transparent"
  , height   = 4
  , width    = 4
  , device   = png
  , scale    = 1
)
