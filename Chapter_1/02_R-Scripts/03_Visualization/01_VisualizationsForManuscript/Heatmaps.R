################################################################################
#### Plot of Heatmaps
################################################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)         # To handle raster data
library(rgdal)          # To load shapefiles
library(tidyverse)      # To wrangle data
library(davidoff)       # Custom functions
library(RColorBrewer)   # For colors
library(sf)             # To plot spatial objects with ggplot
library(ggspatial)      # For north arrow and scale bar
library(cowplot)        # To grab legends
library(ggpubr)         # To arrange ggplots

################################################################################
#### Load Required Data
################################################################################
# Load heatmaps
rasterized <- read_rds("03_Data/03_Results/99_Heatmaps.rds")

# Load shapefiles
buffer  <- readOGR("03_Data/03_Results/99_BufferArea.shp")
main    <- readOGR("03_Data/03_Results/99_SourceAreas.shp")
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
africa  <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
prot    <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Load reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Put heatmaps into a stack
heatmaps <- stack(rasterized$heatmap)

# Crop them to the extent of the main area
heatmaps <- crop(heatmaps, r)

# Prepare a merged map
merged <- heatmaps[[6]] + heatmaps[[12]]

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Convert to sf
kaza    <- st_as_sf(kaza)
africa  <- st_as_sf(africa)

# Convert heatmap to dataframe
merged <- as.data.frame(merged, xy = T)

# Also conver the reference raster to a dataframe
r <- as.data.frame(r, xy = T)

################################################################################
#### Plot
################################################################################
# Prepare color palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

# Main Plot
p1 <- ggplot() +
  geom_raster(
      data    = merged
    , mapping = aes(x = x, y = y, fill = layer)
  ) +
  scale_fill_gradientn(
      colours = myPalette(100)
    , guide   = guide_colorbar(
      , title          = "Number of Traversing Trajectories"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
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
    , title    = "Dispersal Heatmap"
    , subtitle = "After 2000 Steps"
  ) +
  theme(
      legend.position  = "bottom"
    , legend.box       = "vertical"
    , panel.background = element_blank()
    , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 12
      )
  )

# Plot for separate legend
p2 <- ggplot() +
  geom_raster(
      data        = merged
    , mapping     = aes(x = x, y = y, fill = layer)
    , show.legend = F
  ) +
  scale_fill_gradientn(
      colours = myPalette(100)
    , guide   = guide_colorbar(
        title          = "Number of Traversing Trajectories"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  geom_sf(
      data        = kaza
    , mapping     = aes(color = "KAZA-TFCA Borders")
    , fill        = NA
    , lty         = 1
    , lwd         = 0.5
    , show.legend = "line"
  ) +
  geom_sf(
      data        = africa
    , mapping     = aes(color = "Country Borders")
    , fill        = NA
    , lty         = 2
    , lwd         = 0.2
    , show.legend = "line"
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , col      = NULL
    , title    = "Dispersal Heatmap"
    , subtitle = "After 2000 Steps"
  ) +
  scale_color_manual(
      values = c("Country Borders" = "black", "KAZA-TFCA Borders" = "black")
    , guide  = guide_legend(
        override.aes = list(
          linetype = c(2, 1)
        , shape    = c(NA, NA)
        , lwd      = c(0.2, 0.5)
      )
    )
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  theme(
      legend.position       = c(0.20, 0.90)
    , legend.box            = "vertical"
    , legend.background     = element_blank()
    , legend.box.background = element_rect(fill  = "white")
    , legend.margin         = margin(0, 8, 2, 6)
    , legend.text           = element_text(color = "black")
    , legend.key            = element_blank()
    , legend.key.size       = unit(0.8, "lines")
    , legend.key.width      = unit(1.2, "lines")
    , panel.background      = element_blank()
  )

# Extract legend
legend <- get_legend(p2)

# Put into main plot
p3 <- p1 + annotation_custom(
      grob = legend
    , xmin = 18.75
    , xmax = 21
    , ymin = -13
    , ymax = -13.5
  )

# Show it
p3

# Store it to file
ggsave("04_Manuscript/99_Heatmap.png", plot = p3)

################################################################################
#### Function to Plot
################################################################################
# Write a function to plot a heatmap
plotHeatmap <- function(x, subtitle = NULL, legend = T, barwidth = 10){

  # Prepare dataframe
  x <- as.data.frame(x, xy = T)
  names(x) <- c("x", "y", "layer")

  # Main Plot
  p1 <- ggplot() +
    geom_raster(
        data    = x
      , mapping = aes(x = x, y = y, fill = layer)
    ) +
    scale_fill_gradientn(
        colours = myPalette(100)
      , guide   = guide_colorbar(
        , title          = "Number of Traversing Trajectories"
        , show.limits    = T
        , title.position = "top"
        , title.hjust    = 0.5
        , ticks          = T
        , barheight      = unit(0.6, "cm")
        , barwidth       = unit(barwidth, "cm")
      )
    ) +
    geom_sf(
        data        = kaza
      , mapping     = aes(color = "KAZA-TFCA Borders")
      , fill        = NA
      , lty         = 1
      , lwd         = 0.5
      , show.legend = F
    ) +
    geom_sf(
        data        = africa
      , mapping     = aes(color = "Country Borders")
      , fill        = NA
      , lty         = 2
      , lwd         = 0.2
      , show.legend = F
    ) +
    scale_color_manual(
        values = c("KAZA-TFCA Borders" = "black", "Country Borders" = "black")
      , guide  = guide_legend(
          override.aes = list(
            linetype = c(1, 2)
          , shape    = c(NA, NA)
          , lwd     = c(0.5, 0.2)
        )
      )
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
      , title    = "Dispersal Heatmap"
      , subtitle = subtitle
    ) +
    theme(
        legend.position      = "top"
      , legend.justification = "right"
      , legend.box           = "vertical"
      , legend.box.margin    = margin(-10, 0, -10, 0)
      , panel.background     = element_blank()
      , panel.border         = element_rect(colour = "black", fill = NA, size = 1)
      , plot.subtitle        = element_text(margin = margin(t = 0, b = -30))
    ) +
    annotation_scale(
        location   = "bl"
      , width_hint = 0.2
      , line_width = 0.5
      , height     = unit(0.15, "cm")
    ) +
    annotation_north_arrow(
        location = "br"
      , height   = unit(1.5, "cm"),
      , width    = unit(1.2, "cm"),
      , style    = north_arrow_fancy_orienteering(
            fill      = c("black", "black")
          , line_col  = NA
          , text_col  = "black"
          , text_size = 12
        )
    )

  if (legend){

    # Plot for separate legend
    p2 <- ggplot() +
      geom_raster(
          data        = x
        , mapping     = aes(x = x, y = y, fill = layer)
        , show.legend = F
      ) +
      geom_sf(
          data        = kaza
        , mapping     = aes(color = "KAZA-TFCA Borders")
        , fill        = NA
        , lty         = 1
        , lwd         = 0.5
        , show.legend = "line"
      ) +
      geom_sf(
          data        = africa
        , mapping     = aes(color = "Country Borders")
        , fill        = NA
        , lty         = 2
        , lwd         = 0.2
        , show.legend = "line"
      ) +
      scale_color_manual(
          values = c("black", "black")
        , guide  = guide_legend(
            override.aes = list(
              linetype = c(2, 1)
            , shape    = c(NA, NA)
            , lwd      = c(0.2, 0.5)
          )
        )
      ) +
      coord_sf(
          crs    = 4326
        , xlim   = c(min(r$x), max(r$x))
        , ylim   = c(min(r$y), max(r$y))
        , expand = F
      ) +
      theme(
          legend.position       = c(0.20, 0.90)
        , legend.box            = "vertical"
        , legend.background     = element_blank()
        , legend.box.background = element_rect(fill  = "white")
        , legend.margin         = margin(0, 8, 2, 6)
        , legend.text           = element_text(color = "black")
        , legend.key            = element_blank()
        , legend.key.size       = unit(0.8, "lines")
        , legend.key.width      = unit(1.2, "lines")
        , panel.background      = element_blank()
      )

    # Extract legend
    legend <- get_legend(p2)

    # Put into main plot
    p3 <- p1 + annotation_custom(
          grob = legend
        , xmin = 18.75
        , xmax = 21
        , ymin = -13
        , ymax = -13.5
      )
  } else {
    p3 <- p1
  }

  # Return it
  return(p3)

}

# Let's apply the function to get all desired plots
maps <- lapply(1:nlayers(heatmaps), function(x){
  subtitle <- paste0("After ", rasterized$steps[x], " Steps")
  map <- plotHeatmap(
      x        = heatmaps[[x]]
    , subtitle = subtitle
    , legend   = F
    , barwidth = 7
  )
  return(map)
})

# Arrange plots nicely
p <- ggarrange(maps[[2]], maps[[4]], maps[[6]], maps[[8]], maps[[10]], maps[[12]])

# Store the arranged plot
ggsave("04_Manuscript/99_HeatmapsIndividual.png"
  , plot   = p
  , scale  = 2
  , height = 6
  , width  = 9
)
