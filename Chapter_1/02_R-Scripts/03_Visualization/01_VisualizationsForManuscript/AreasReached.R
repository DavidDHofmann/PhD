################################################################################
#### Visualize Connections Between National Parks
################################################################################
# Description: In this script, we visualize the interpatch connectivity

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(igraph)         # For network analysis
library(rgdal)          # To load spatial data
library(sf)             # To plot spatial stuff with ggplot
library(ggspatial)      # To add scale bars etc to plots
library(rgeos)          # To calculate areas
library(ggnetwork)      # To plot network using ggplot
library(viridis)        # For nice colors
library(davidoff)       # Custom functions
library(cowplot)        # Grab ggplot legend
library(gtable)         # To combine ggplots
library(grid)           # To combine ggplots
library(gridExtra)      # To combine ggplots

################################################################################
#### Prepare Network
################################################################################
# Reload data on areas reached
visits <- read_rds("03_Data/03_Results/99_AreasReached.rds")

# Load protected areas
prot    <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Keep only national parks
prot <- subset(prot, Desig == "National Park")

# Assign a unique ID to each of the areas
prot$ID <- 1:nrow(prot)

# Identify the area of each protected area
prot$Area <- gArea(spTransform(prot, CRS("+init=epsg:32734")), byid = T)

# Coerce the visitations to a graph
net <- graph_from_data_frame(
    d        = visits
  , vertices = unique(prot$ID)
  , directed = T
)

# Add area as vertex information
V(net)$Area <- prot$Area

# Prepare layout
lay <- coordinates(gCentroid(prot, byid = T))

# Prepare plot for ggplotting
net_p <- ggnetwork(net, layout = lay, arrow.gap = 0.1, scale = F)

################################################################################
#### Prepare Additional Plotting Data
################################################################################
# Load additional shapefiles and raster for the background
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
africa  <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
prot    <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
r       <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Simplify Protection zones
prot$Desig <- as.character(prot$Desig)
prot$Desig[prot$Desig == "Forest Reserve"] <- "Protected"

# Make nicer names
prot$Desig[prot$Desig == "National Park"] <- "National Parks"
prot$Desig[prot$Desig == "Protected"] <- "Protected Areas"
prot$Desig <- as.factor(as.character(prot$Desig))

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Rename
kaza$Name <- "KAZA-TFCA Borders"

# Convert objects to sf
kaza                  <- st_as_sf(kaza)
africa                <- st_as_sf(africa)
prot                  <- st_as_sf(prot)

# Convert heatmap to dataframe
r <- as.data.frame(r, xy = T)

################################################################################
#### Plot
################################################################################
# Prepare color palette
pal <- colorRampPalette(plasma(100, begin = 0.8, end = 0))

# Main Plot
p1 <- ggplot() +
  geom_sf(
      data        = prot
    , mapping     = aes(fill = Desig)
    , col         = "#6ba36b"
    , lwd         = 0
    , show.legend = F
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
  geom_edges(
      data      = net_p
    , mapping   = aes(x = x, y = y, xend = xend, yend = yend
      , size = RelFrequency, col = MeanStepNumber)
    , curvature = 0.2
    , arrow     = arrow(length = unit(6, "pt"), type = "closed", angle = 10)
  ) +
  scale_size_area(
      name     = "Visitation Frequency"
    , max_size = 1
  ) +
  scale_color_gradientn(
      colors  = pal(100)
    , guide   = guide_colorbar(
        title          = "Number of Steps"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(3.0, "cm")
    )
  ) +
  scale_fill_manual(
    values = c("#70ab70", "#d9f0d3")
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
    , title    = "Areas Reached and Visitation Frequency"
    , subtitle = "In Relation to Number of Steps"
  ) +
  guides(
    size = guide_legend(title.position = "top")
  ) +
  theme(
      legend.position      = "bottom"
    , legend.box           = "horizontal"
    , legend.title.align   = 0.5
    , panel.background     = element_blank()
    , panel.border         = element_rect(colour = "black", fill = NA, size = 1)
    , legend.title         = element_text(size = 10),
    , legend.text          = element_text(size = 8)
    , legend.margin        = margin(c(0, 0, 0, 0))
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
  geom_sf(
      data        = prot
    , mapping     = aes(fill = Desig)
    , col         = "#6ba36b"
    , lwd         = 0
  ) +
  geom_sf(
      data        = kaza
    , mapping     = aes(col = "KAZA-TFCA Borders")
    , fill        = NA
    , show.legend = "line"
    , lwd         = 0.5
  ) +
  geom_sf(
      data        = africa
    , mapping     = aes(col = "Country Borders")
    , fill        = NA
    , show.legend = "line"
    , lwd         = 0.2
  ) +
  scale_fill_manual(
    values = c("National Parks" = "#70ab70", "Protected Areas" = "#d9f0d3")
    , guide = guide_legend(
      override.aes = list(
          linetype = c("blank", "blank")
        , shape    = c(NA, NA)
      )
    )
  ) +
  scale_color_manual(
      values = c("Country Borders" = "black", "KAZA-TFCA Borders" = "black")
    , guide = guide_legend(
        override.aes = list(
            linetype = c(2, 1)
          , lwd      = c(0.2, 0.5)
        )
      )
  ) +
  theme(
      legend.title          = element_blank()
    , legend.spacing.y      = unit(0, "cm")
    , legend.box            = "vertical"
    , legend.background     = element_blank()
    , legend.box.background = element_rect(fill  = "white")
    , legend.margin         = margin(2, 8, 2, 6)
    , legend.text           = element_text(color = "black")
    , legend.key            = element_blank()
    , legend.key.size       = unit(0.8, "lines")
    , legend.key.width      = unit(1.2, "lines")
    , panel.background      = element_blank()
  )

# Create plot with dots of different size for yet another legend
p3 <- ggplot() +
  geom_point(
      data    = net_p
    , mapping = aes(x = x, y = y)
    , col     = "black"
    , size    = 0.1
  ) +
  geom_point(
      data    = net_p
    , mapping = aes(x = x, y = y, size = Simulations)
    , col     = "orange"
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
    , title    = "Areas Reached and Visitation Frequency"
    , subtitle = "In Relation to Number of Steps"
  ) +
  guides(
    size = guide_legend(title.position = "top")
  ) +
  theme(
      legend.position      = "bottom"
    , legend.box           = "horizontal"
    , legend.title.align   = 0.5
    , panel.background     = element_rect(fill = "transparent", colour = NA)
    , panel.grid           = element_blank()
    , legend.margin        = margin(c(0, 0, 13, 0))
    , legend.title         = element_text(size = 10),
    , legend.text          = element_text(size = 8)
  )

# Extract legends from all plots
legend1 <- get_legend(p1)
legend2 <- get_legend(p2)
legend3 <- get_legend(p3)

# Put legend with protected areas into main plot
p4 <- p1 + annotation_custom(
      grob = legend2
    , xmin = 18.75
    , xmax = 21
    , ymin = -13
    , ymax = -14
  )

# Remove the other legend from main plot
p5 <- p4 + theme(legend.position = "none")

# Add circles to main plot
g1 <- ggplotGrob(p5)
g2 <- ggplotGrob(p3)
g2 <- gtable_filter(g2, "panel")
pos <- c(subset(g1$layout, grepl("panel", g1$layout$name), select = t:r))
p5 <- gtable_add_grob(g1, g2, t = pos$t, l = pos$l)
p5 <- ggplotify::as.ggplot(p5)

# Put legends together
legends <- grid.arrange(legend1, legend3, nrow = 1, widths = c(2, 1))
legends <- gtable_add_padding(legends, unit(c(0, 1, 0, 0.3), "cm"))

# Put them below the first plot
p6 <- arrangeGrob(p5, legends, heights = c(10, 1))
p6 <- ggplotify::as.ggplot(p6)

# Visualize
p6

# Store
ggsave("04_Manuscript/99_AreasReached.png", plot = p6)
