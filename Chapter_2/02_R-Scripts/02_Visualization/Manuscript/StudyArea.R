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
library(terra)        # To handle raster files
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

# Get the reference raster and derive an extent from it
r <- raster("03_Data/02_CleanData/ReferenceRaster.tif")
extent <- as(extent(r), "SpatialPolygons")
crs(extent) <- CRS("+init=epsg:4326")

# Expanded study area (buffer zone)
buffer <- as(extent(r) * 1.2, "SpatialPolygons")
crs(buffer) <- CRS("+init=epsg:4326")
buffer <- buffer - extent

# Define map borders
xlims <- extent(buffer)[1:2]
ylims <- extent(buffer)[3:4]

# Convert to sf
buffer <- st_as_sf(buffer)
extent <- st_as_sf(extent)
r <- as.data.frame(r, xy = T)

# Load stuff that we would like to plot
africa <- read_sf("03_Data/02_CleanData/Africa.shp")
prot   <- read_sf("03_Data/02_CleanData/Protected.shp")
water  <- read_sf("03_Data/02_CleanData/MajorWaters.shp")
river  <- read_sf("03_Data/02_CleanData/MajorRivers.shp")
areas  <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
roads  <- read_sf("03_Data/02_CleanData/Roads.shp")
fault  <- read_sf("03_Data/02_CleanData/Faults.shp")
mabab  <- read_sf("03_Data/02_CleanData/MababeDepression.shp") %>% st_cast("LINESTRING")
kaza   <- read_sf("03_Data/02_CleanData/KAZA.shp")
vills  <- read_sf("03_Data/02_CleanData/Villages.shp")
vills  <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y) %>%
  mutate(place = gsub(place, pattern = "Village", replacement = "Villages")) %>%
  mutate(place = gsub(place, pattern = "City", replacement = "Cities"))

# Crop cutlines
cutli <- vect("03_Data/02_CleanData/Cutlines.shp")
cutli_main <- cutli[c(4, 8), ]
cutli_rest <- areas %>%
  subset(Type == "Buffer") %>%
  as("Spatial") %>%
  vect() %>%
  buffer(width = 12000) %>%
  aggregate() %>%
  intersect(cutli[-c(4, 8)], .)
cutli <- rbind(cutli_main, cutli_rest)
cutli <- st_as_sf(cutli)

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
    x     = c(26.3, 22.35, 23.67, 24.51, 24.8, 23.35)
  , y     = c(-19.08, -17.70, -19.35, -18.65, -20.0, -21.3)
  , Label = paste0(c("Hwange", "Luengue-Luiana", "Moremi", "Chobe", "Nxai Pan", "Central Kalahari"), "\nNP")
)

# Create labels for rivers
labels_rivers <- data.frame(
    x     = c(20.65, 23.18, 23.15, 24, 24.1, 22.40)
  , y     = c(-17.85, -18.65, -20.35, -20.20, -19.21, -19.55)
  , Label = c("Okavango", "Selinda", "Nhabe", "Boteti", "Khwai", "Thaoge")
)

# Add fault labels
labels_faults <- data.frame(
    x     = c(21.8, 23.78, 24.17)
  , y     = c(-20.72, -19.7, -18.88)
  , Label = c("Kunyere Fault", "Thamalakane Fault", "Mababe Depression")
)

# Generate area labels
labels_areas <- st_coordinates(st_point_on_surface(areas))
labels_areas <- cbind(labels_areas, st_drop_geometry(areas))

# Load (updated) dispersal trajectories
# disp <- read_csv("03_Data/02_CleanData/DispersersUpdated.csv")
# disp <- subset(disp, !is.na(x))
# disp <- mutate(disp, Dispersing = State %in% c("Explorer", "Disperser"))
# disp <- subset(disp, Dispersing)

# Prepare plot of africa
p1 <- ggplot() +
  geom_sf(data = africa, fill = "gray90", col = "white", linewidth = 0.2) +
  geom_sf(data = kaza, fill = "#D62828", col = "#D62828", alpha = 0.2) +
  geom_sf(data = st_union(buffer, extent), fill = "black", col = "black", alpha = 0.3) +
  coord_sf(xlim = c(12, 40), ylim = c(-35, -5), clip = "on") +
  theme_void()

# Prepare plot of study area
p2 <- ggplot() +
  geom_sf(data = buffer, fill = "gray80", col = "gray80", alpha = 0.4, lwd = 0) +
  geom_sf(data = prot, aes(fill = Desig), col = NA, alpha = 0.7) +
  geom_sf(data = water, fill = "cornflowerblue", col = NA) +
  geom_sf(data = river, col = "cornflowerblue") +
  geom_sf(data = kava, col = "cornflowerblue") +
  geom_sf(data = roads, col = "gray70", lwd = 0.2) +
  geom_sf(data = st_crop(africa, extent), col = "black", fill = NA, lwd = 0.3) +
  geom_sf(data = fault, col = "gray30", lty = 4) +
  geom_sf(data = mabab, col = "gray30", fill = NA, lty = 4) +
  geom_sf(data = cutli, lty = 2, col = "purple") +
  geom_point(
      data        = vills
    , mapping     = aes(x = x, y = y, size = place)
    , col         = "gray50"
    , shape       = 15
    , show.legend = F
  ) +
  # geom_path(
  #     data    = disp
  #   , mapping = aes(x = x, y = y, group = as.factor(DogName))
  #   , col     = viridis::viridis(20)[20]
  #   , size    = 0.3
  # ) +
  geom_sf(
      data  = subset(areas, Type == "Main")
    , col   = "orange"
    , fill  = "orange"
    , alpha = 0.2
  ) +
  geom_sf(
      data  = subset(areas, Type == "Buffer")
    , col   = "purple"
    , fill  = "purple"
    , alpha = 0.2
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
    , size     = 2.25
    , angle    = c(-30, 40, 20, 45, 45, 80)
  ) +
  geom_text(
      data     = subset(labels_nationalparks, Label != "Moremi\nNP")
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
      data     = subset(vills, place == "Cities")
    , mapping  = aes(x = x, y = y, label = name)
    , col      = "gray50"
    , fontface = 3
    , size     = 3
    , nudge_y  = c(0.1, -0.1, 0.1)
  ) +
  geom_text(
      data     = subset(vills, name == "Nokaneng")
    , mapping  = aes(x = x, y = y, label = name)
    , col      = "gray50"
    , fontface = 3
    , size     = 1.8
    , nudge_y  = 0.05
    , nudge_x  = -0.05
  ) +
  geom_text(
      data     = labels_faults
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "gray30"
    , fontface = 3
    , size     = 1.5
    , angle = c(-5, 52, 45)
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
  geom_text(
      data     = subset(labels_areas, Type == "Main")
    , mapping  = aes(x = X, y = Y, label = ID)
    , col      = "black"
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
    , pad_y      = unit(0.5, "cm")
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , pad_x    = unit(0.3, "cm")
    , pad_y    = unit(0.3, "cm")
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
  ylab("")  +
  coord_sf(
      crs    = 4326
    , xlim   = xlims
    , ylim   = ylims
    , expand = F
  )

# +
#   coord_sf(
#       crs    = 4326
#     , xlim   = c(min(r$x), max(r$x))
#     , ylim   = c(min(r$y), max(r$y))
#     , expand = F
#   )

################################################################################
#### Legends
################################################################################
# Create a new plot for the legend
df <- prot[1:7, ]
df <- dplyr::select(df, geometry)
df <- cbind(
    df
  , Name  = factor(c("National Parks", "Forest Reserves", "Other Protected Areas", "Major Waters", "Source Areas", "Egression Zones", "KAZA-TFCA")
    , levels = c("National Parks", "Forest Reserves", "Other Protected Areas", "Major Waters", "Source Areas", "Egression Zones", "KAZA-TFCA"))
    , Color = c(rev(brewer.pal(n = 3, name = "Greens")), "cornflowerblue", "orange", "purple", "#D62828")
    , Alpha = c(0.7, 0.7, 0.7, 1, 0.2, 0.2, 0.2)
)
df2 <- data.frame(
    x     = c(1, 1, 2, 1)
  , y     = c(1, 1, 2, 1)
  , Class = as.factor(c("Dispersers", "Cutlines", "Faults / Depressions", "Roads"))
)

# Plot
l1 <- ggplot() +
  geom_sf(data = df, aes(fill = Name, alpha = Name, col = Name)) +
  scale_fill_manual(values = df$Color, name = "Legend") +
  scale_color_manual(values = c("white", "white", "white", "white", "orange", "purple", "#D62828"), name = "Legend") +
  scale_alpha_manual(values = c(0.7, 0.7, 0.7, 1.0, 0.2, 0.2, 0.2), name = "Legend") +
  theme_minimal()
l2 <- ggplot() +
  geom_point(data = vills, aes(x = x, y = y, size = place), shape = 15, col = "gray50") +
  geom_line(data = df2, aes(x = x, y = y, col = Class, lty = Class)) +
  scale_color_manual(values = c("purple", viridis::viridis(20)[20], "gray30", "gray70"), name = "") +
  scale_size_manual(values = c(1, 0.3, 0.5, 2), name = "") +
  scale_linetype_manual(values = c(2, 1, 4, 1), name = "") +
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
p4 <- ggarrange(p3, p2, nrow = 1, widths = c(2, 5), labels = "auto")

# Store the plot
ggsave("04_Manuscript/Figures/StudyArea.png"
  , plot   = p4
  , bg     = "white"
  , width  = 8
  , height = 4
  , scale  = 1.6
  , device = png
)
