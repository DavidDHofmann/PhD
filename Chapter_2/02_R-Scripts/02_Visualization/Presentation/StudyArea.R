################################################################################
#### Study Area
################################################################################
# Description: Plots of the study area

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
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

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Generate directory
dir.create("05_Presentation", showWarnings = F)
dir.create("05_Presentation/Figures", showWarnings = F)

# Get the reference raster and derive an extent from it
r <- raster("03_Data/02_CleanData/ReferenceRaster.tif")
extent <- as(extent(r), "SpatialPolygons")
crs(extent) <- CRS("+init=epsg:4326")

# Expanded study area (buffer zone)
buffer <- as(extent(r) * 1.2, "SpatialPolygons")
crs(buffer) <- CRS("+init=epsg:4326")
buffer <- buffer - extent

# Convert to sf
buffer <- st_as_sf(buffer)
extent <- st_as_sf(extent)
r <- as.data.frame(r, xy = T)

library(rgeos)
africa <- shapefile("03_Data/02_CleanData/Africa.shp")
africa2 <- shapefile("03_Data/02_CleanData/Africa.shp")

# Buffer slightly
africa <- gBuffer(africa, width = 0.01)

# Disaggregate
africa <- disaggregate(africa)

# Identify sizes of areas
africa$Size <- gArea(africa, byid = T)

# Keep only the largest two
africa <- subset(africa, Size %in% sort(africa$Size, decreasing = T)[1:2])

# Simplify and smoothen africa shape
africa <- gSimplify(africa, tol = 0.5)
library(smoothr)
africa <- smooth(africa, method = "ksmooth")

# Clip africa layer
africa <- gIntersection(africa2, africa, byid = T)
africa <- st_as_sf(africa)

# Load stuff that we would like to plot
prot   <- read_sf("03_Data/02_CleanData/Protected.shp")
water  <- read_sf("03_Data/02_CleanData/MajorWaters.shp")
river  <- read_sf("03_Data/02_CleanData/MajorRivers.shp")
areas  <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
roads  <- read_sf("03_Data/02_CleanData/Roads.shp")
fault  <- read_sf("03_Data/02_CleanData/Faults.shp")
vills  <- read_sf("03_Data/02_CleanData/Villages.shp")
cutli  <- read_sf("03_Data/02_CleanData/Cutlines.shp")
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
    x     = c(22.6, 23.7, 22.7, 25.6)
  , y     = c(-19, -18.1, -20.4, -20.7)
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

# Generate area labels
labels_areas <- st_coordinates(st_point_on_surface(areas))
labels_areas <- cbind(labels_areas, st_drop_geometry(areas))

# Load dispersal trajectories
disp <- read_csv("03_Data/02_CleanData/DispersersUpdated.csv")
disp <- subset(disp, State == "Disperser")

# Prepare plot of africa
p1 <- ggplot() +
  geom_sf(data = africa, fill = "white", col = "white", linewidth = 1.5) +
  geom_sf(data = africa, fill = "gray10", col = "gray50", lwd = 0.7) +
  geom_sf(data = buffer, fill = "red", col = "red", alpha = 0.4) +
  geom_sf(data = extent, fill = "gray30", col = "gray30", alpha = 0.4) +
  theme_minimal() +
  coord_sf(xlim = c(-18, 48), ylim = c(-32, 35), clip = "off") +
  theme(
      panel.grid = element_blank()
    , axis.ticks = element_blank()
    , axis.text  = element_blank()
    , plot.background = element_blank()
  )

# Prepare plot of study area
p2 <- ggplot() +
  geom_sf(data = prot, aes(fill = Desig), col = NA, alpha = 1) +
  geom_sf(data = water, fill = "cornflowerblue", col = NA) +
  geom_sf(data = river, col = "cornflowerblue") +
  geom_sf(data = kava, col = "cornflowerblue") +
  geom_sf(data = roads, col = "gray70", linewidth = 0.2) +
  geom_sf(data = africa, col = "white", fill = NA, linewidth = 1) +
#   geom_sf(data = fault, col = "red", lty = 4) +
#   geom_point(
#       data        = vills
#     , mapping     = aes(x = x, y = y, size = place)
#     , col         = "gray80"
#     , shape       = 15
#     , show.legend = F
#   ) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("cornflowerblue", 1.8)
    , fontface = 3
    , size     = 3
  ) +
#   geom_text(
#       data     = labels_rivers
#     , mapping  = aes(x = x, y = y, label = Label)
#     , col      = darken("cornflowerblue", 1.8)
#     , fontface = 3
#     , size     = 1.5
#     , angle    = c(-30, 40, 20, -32)
#   ) +
  geom_text(
      data     = subset(labels_nationalparks, Label != "Moremi\nNP")
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("green", 2.5)
    , fontface = 3
    , size     = 3
  ) +
  geom_text(
      data     = labels_countries
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "white"
    , fontface = 2
  ) +
  geom_text(
      data     = subset(vills, place == "City")
    , mapping  = aes(x = x, y = y, label = name)
    , col      = "gray80"
    , fontface = 3
    , size     = 3
    , nudge_y  = c(0.1, -0.1, 0.1)
  ) +
#   geom_text(
#       data     = labels_faults
#     , mapping  = aes(x = x, y = y, label = Label)
#     , col      = "red"
#     , fontface = 3
#     , size     = 1.5
#     , angle = c(-5, 52)
#   ) +
  scale_fill_brewer(palette = "Greens", direction = -1, name = "") +
  scale_size_manual(values = c(2.0, 0.5)) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
    , pad_x      = unit(0.9, "cm")
    , pad_y      = unit(0.8, "cm")
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , pad_x    = unit(0.7, "cm")
    , pad_y    = unit(0.6, "cm")
    , style    = north_arrow_fancy_orienteering(
          fill      = c("white", "white")
        , line_col  = NA
        , text_col  = "white"
        , text_size = 12
      )
  ) +
  theme(
      legend.position  = "none"
    , legend.key       = element_rect(color = NA)
    , axis.text        = element_text(color = "white")
    , axis.ticks       = element_line(color = "white")
    , panel.border     = element_rect(color = "white", fill = NA, size = 1)
    , panel.grid.major = element_line(color = NA)
    , panel.grid.minor = element_line(color = NA)
    , panel.background = element_rect(fill   = "transparent")
    , plot.background  = element_rect(fill   = "transparent", color = NA)
  ) +
  xlab("") +
  ylab("") +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  )

# Add source areas
p3 <- p2 + geom_sf(
      data  = subset(areas, Type == "Main")
    , col   = "orange"
    , fill  = "orange"
    , alpha = 0.8
  ) +
  geom_text(
      data     = subset(labels_areas, Type == "Main")
    , mapping  = aes(x = X, y = Y, label = ID)
    , col      = "white"
    , fontface = 3
    , size     = 2
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  )

# Add emigration zones
p4 <- p3 +
  geom_sf(data = cutli, lty = 2, col = "blue") +
  geom_sf(
      data  = subset(areas, Type == "Buffer")
    , col   = "purple"
    , fill  = "purple"
    , alpha = 0.8
  ) +
  geom_point(
      data     = subset(labels_areas, Type == "Buffer")
    , mapping  = aes(x = X, y = Y)
    , col      = "purple"
    , size     = 4
  ) +
  geom_text(
      data     = subset(labels_areas, Type == "Buffer")
    , mapping  = aes(x = X, y = Y, label = ID)
    , col      = "white"
    , fontface = 3
    , size     = 2
  ) +
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
df <- prot[1:6, ]
df <- dplyr::select(df, geometry)
df <- cbind(
    df
  , Name  = factor(c("National Parks", "Forest Reserves", "Protected", "Major Waters", "Source Areas", "Emigration Zones")
    , levels = c("National Parks", "Forest Reserves", "Protected", "Major Waters", "Source Areas", "Emigration Zones"))
  , Color = c(rev(brewer.pal(n = 3, name = "Greens")), "cornflowerblue", "orange", "purple")
)
df2 <- data.frame(
    x     = c(1, 1)
  , y     = c(1, 1)
  , Class = c("Cutlines", "Roads")
)

# Plot
l1 <- ggplot() +
  geom_sf(data = df, aes(fill = Name, alpha = Name, col = Name)) +
  scale_fill_manual(values = df$Color, name = "Legend") +
  scale_color_manual(values = c("transparent", "transparent", "transparent", "transparent", "orange", "purple"), name = "Legend") +
  scale_alpha_manual(values = c(1.0, 1.0, 1.0, 1.0, 0.8, 0.8), name = "Legend") +
  theme_minimal() +
  theme(
      legend.text       = element_text(color = "gray50")
    , legend.title      = element_text(color = "gray50")
    , legend.key        = element_rect(color = NA, fill = NA)
    , legend.background = element_rect(fill = "transparent", color = NA)
  )
l2 <- ggplot() +
#   geom_point(data = vills, aes(x = x, y = y, size = place), shape = 15, col = "gray80") +
  geom_line(data = df2, aes(x = x, y = y, col = Class, lty = Class)) +
  scale_color_manual(values = c("blue", "gray70"), name = "") +
  scale_size_manual(values = c(1, 0.5, 2), name = "") +
  scale_linetype_manual(values = c(2, 4, 1), name = "") +
  theme_minimal() +
  theme(
      legend.text  = element_text(color = "gray50")
    , legend.title = element_text(color = "gray50")
  )

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
p1 <- ggarrange(p1, legend, nrow = 2, heights = c(2, 1))
p2 <- ggarrange(p1, p2, nrow = 1, widths = c(2, 5), labels = "auto", font.label = list(color = "gray50"))
p3 <- ggarrange(p1, p3, nrow = 1, widths = c(2, 5), labels = "auto", font.label = list(color = "gray50"))
p4 <- ggarrange(p1, p4, nrow = 1, widths = c(2, 5), labels = "auto", font.label = list(color = "gray50"))

# Store the plots
ggsave("05_Presentation/Figures/StudyArea1.png"
  , plot   = p2
  , bg     = "transparent"
  , width  = 8
  , height = 4
  , scale  = 1.4
)
ggsave("05_Presentation/Figures/StudyArea2.png"
  , plot   = p3
  , bg     = "transparent"
  , width  = 8
  , height = 4
  , scale  = 1.4
)
ggsave("05_Presentation/Figures/StudyArea3.png"
  , plot   = p4
  , bg     = "transparent"
  , width  = 8
  , height = 4
  , scale  = 1.4
)
