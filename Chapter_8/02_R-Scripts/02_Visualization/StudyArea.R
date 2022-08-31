################################################################################
#### Study Area
################################################################################
# Description: Plots of the study area

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
setwd(wd)

# Load required packages
library(tidyverse)    # For data wrangling
library(raster)       # To handle raster files
library(sf)           # To handle spatial data
library(ggplot2)      # For nice plots
library(ggspatial)    # To add scale bars and north arrows
library(ggpubr)       # To extract plot legends
library(grid)         # To arrange plots
library(gridExtra)    # To arrange plots
library(gtable)       # To arrange plots
library(RColorBrewer) # For colors
library(ggpattern)    # For hatched fills

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Get the reference raster and derive an extent from it
r <- raster("03_Data/02_CleanData/ReferenceRaster.tif")
extent <- as(extent(r), "SpatialPolygons")
crs(extent) <- CRS("+init=epsg:4326")
extent <- st_as_sf(extent)
r <- as.data.frame(r, xy = T)

# Load stuff that we would like to plot
africa <- read_sf("03_Data/02_CleanData/Africa.shp")
prot   <- read_sf("03_Data/02_CleanData/Protected.shp")
water  <- read_sf("03_Data/02_CleanData/MajorWaters.shp")
river  <- read_sf("03_Data/02_CleanData/MajorRivers.shp")
areas  <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
roads  <- read_sf("03_Data/02_CleanData/Roads.shp")
fault  <- read_sf("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/01_RawData/DAVID/Faults.shp")
vills  <- read_sf("03_Data/02_CleanData/Villages.shp")
vills  <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y)

# I also want to add the Okavango river
kava <- raster("03_Data/01_RawData/MERIT/Rivers.tif")
kava <- crop(kava, extent(c(20.5, 22, -18.2, -17.5)))
kava <- rasterToPolygons(kava, fun = function(x) {x == 1})
kava <- aggregate(kava)
kava <- buffer(kava, width = 0.1 / 111, dissolve = T)
kava <- st_as_sf(kava)

# Keep only faults of interest
fault <- fault[1:2, ]

# Reorder the levels of national parks
prot$Desig <- factor(prot$Desig, levels = c("National Park", "Forest Reserve", "Protected"))

# Create country labels
labels_countries <- data.frame(
    x     = c(25.5, 26, 25.7, 21.5, 23.5)
  , y     = c(-19.3, -18.2, -17.6, -17.6, -17.8)
  , Label = c("BOTSWANA", "ZIMBABWE", "ZAMBIA", "ANGOLA", "NAMIBIA")
)

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 22.8, 25.6)
  , y     = c(-19, -18.2, -20.5, -20.7)
  , Label = c("Okavango\nDelta", "Linyanti\nSwamp", "Lake\nNgami", "Makgadikgadi\nPans")
)

# Create labels for some national parks and source areas
labels_nationalparks <- data.frame(
    x     = c(26.3, 22.35, 23.67, 24.51, 24.7, 23.35)
  , y     = c(-19.08, -17.70, -19.35, -18.65, -20.4, -21.3)
  , Label = paste0(c("Hwange", "Luengue-Luiana", "Moremi", "Chobe", "Nxai Pan", "Central Kalahari"), "\nNP")
)

# Create labels for rivers
labels_rivers <- data.frame(
    x     = c(20.65, 23.18, 23.15, 23.9)
  , y     = c(-17.85, -18.65, -20.35, -20.25)
  , Label = c("Okavango", "Selinda", "Nhabe", "Boteti")
)

# Add fault labels
labels_faults <- data.frame(
    x     = c(21.8, 23.78)
  , y     = c(-20.72, -19.7)
  , Label = c("Kunyere Fault", "Thamalakane Fault")
)

# Add albels for source areas
labels_source <- data.frame(
    x     = c(23.4, 22, 23.5, 22.7, 23.1, 23.4, 21.55, 24.1, 22.6)
  , y     = c(-20.3, -19.4, -18.9, -18.3, -19.4, -20.5, -19.5, -19.5, -18)
  , Label = 1:9
)

# Prepare plot of africa
p1 <- ggplot() +
  geom_sf(data = africa, fill = "gray90", col = "white", lwd = 0.1) +
  geom_sf(data = extent, fill = "gray30", col = "gray30", alpha = 0.2) +
  theme_minimal() +
  coord_sf(xlim = c(-18, 48), ylim = c(-32, 35)) +
  theme(
      panel.grid = element_blank()
    , axis.ticks = element_blank()
    , axis.text  = element_blank()
  )

# Prepare plot of study area
p2 <- ggplot() +
  geom_sf(data = prot, aes(fill = Desig), col = NA, alpha = 0.7) +
  geom_sf(data = water, fill = "cornflowerblue", col = NA) +
  geom_sf(data = river, col = "cornflowerblue") +
  geom_sf(data = kava, col = "cornflowerblue") +
  geom_sf(data = roads, col = "gray70", lwd = 0.2) +
  geom_sf(data = africa, col = "black", fill = NA, lwd = 0.3) +
  geom_sf(data = areas, col = "orange", fill = "orange", alpha = 0.2) +
  # geom_sf(data = areas, col = "black", fill = "black", alpha = 0.2) +
  geom_sf_pattern(data = areas,
    aes(pattern = as.factor(ID), pattern_spacing = as.factor(ID), pattern_angle = as.factor(ID))
    , col = "orange"
    , fill = NA
    , pattern_color = "orange"
    , pattern_fill = NA
    , pattern_density = 0.01
    , pattern_spacing = 0.01
    , pattern_size = 0.2
  ) +
  geom_sf(data = fault, col = "gray30", lty = 2) +
  geom_point(
      data        = vills
    , mapping     = aes(x = x, y = y, size = place)
    , col         = "gray50"
    , shape       = 15
    , show.legend = F
  ) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("cornflowerblue", 1.4)
    , fontface = 3
    , size     = 3
  ) +
  geom_text(
      data     = labels_rivers
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("cornflowerblue", 1.4)
    , fontface = 3
    , size     = 1.5
    , angle    = c(-30, 40, 20, -32)
  ) +
  geom_text(
      data     = labels_nationalparks
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("green", 1.8)
    , fontface = 3
    , size     = 3
  ) +
  geom_text(
      data     = labels_countries
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "black"
    , fontface = 2
  ) +
  geom_text(
      data     = subset(vills, place == "City")
    , mapping  = aes(x = x, y = y, label = name)
    , col      = "gray50"
    , fontface = 3
    , size     = 3
    , nudge_y  = c(0.1, -0.1, 0.1)
  ) +
  geom_text(
      data     = labels_faults
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "gray30"
    , fontface = 3
    , size     = 1.5
    , angle = c(-5, 52)
  ) +
  geom_point(
      data     = labels_source
    , mapping  = aes(x = x, y = y)
    , col      = "orange"
    , size     = 4
  ) +
  geom_text(
      data     = labels_source
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "white"
    , fontface = 3
    , size     = 2
  ) +
  scale_fill_brewer(palette = "Greens", direction = -1, name = "") +
  scale_size_manual(values = c(2.0, 0.5)) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
    , pad_x      = unit(0.9, "cm")
    , pad_y      = unit(0.8, "cm")
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , pad_x    = unit(0.7, "cm")
    , pad_y    = unit(0.6, "cm")
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 12
      )
  ) +
  theme_minimal() +
  theme(legend.position = "none", panel.border = element_rect(colour = "gray30", fill = NA, size = 1)) +
  xlab("") +
  ylab("") +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  )

################################################################################
#### Legends
################################################################################
# Create a new plot for the legend
df <- prot[1:5, ]
df <- dplyr::select(df, geometry)
df <- cbind(
    df
  , Name  = factor(c("National Parks", "Forest Reserves", "Protected", "Major Waters", "Source Areas")
    , levels = c("National Parks", "Forest Reserves", "Protected", "Major Waters", "Source Areas"))
  , Color = c(rev(brewer.pal(n = 3, name = "Greens")), "cornflowerblue", "orange")
  , Alpha = c(0.7, 0.7, 0.7, 1, 0.2)
)

# Plot
l1 <- ggplot() +
  geom_sf(data = df, aes(fill = Name, alpha = Name, col = Name)) +
  scale_fill_manual(values = df$Color, name = "Legend") +
  scale_color_manual(values = c("white", "white", "white", "white", "orange"), name = "Legend") +
  scale_alpha_manual(values = c(0.7, 0.7, 0.7, 1.0, 0.2), name = "Legend") +
  theme_minimal()
l2 <- ggplot() +
  geom_point(data = vills, aes(x = x, y = y, size = place), shape = 15, col = "gray50") +
  geom_line(data = data.frame(x = c(1, 2, 1, 2), y = c(1, 2, 1, 2), Class = c("Road", "Fault")), aes(x = x, y = y, col = Class, lty = Class)) +
  scale_color_manual(values = c("gray30", "gray70"), name = "") +
  scale_size_manual(values = c(2, 0.5), name = "") +
  scale_linetype_manual(values = c(2, 1), name = "") +
  theme_minimal() +
  theme(legend.spacing.y = unit(-0.38, "cm"))

# Grab legends
l1 <- get_legend(l1)
l2 <- get_legend(l2)

################################################################################
#### Combine Plots
################################################################################
# Arrange the legends
legend <- grid.arrange(l1, l2, nrow = 1)
legend <- gtable_add_padding(legend, unit(c(-2, -1, 0, 0), "cm"))

# Arrange legend with plots
p3 <- ggarrange(p1, legend, nrow = 2, heights = c(2, 1))
p4 <- ggarrange(p3, p2, nrow = 1, widths = c(2, 5))

# Store the plot
ggsave("04_Manuscript/99_StudyArea.png"
  , plot   = p4
  , bg     = "white"
  , width  = 8
  , height = 4
  , scale  = 1.4
)
