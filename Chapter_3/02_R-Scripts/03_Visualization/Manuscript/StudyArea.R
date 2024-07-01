################################################################################
#### Study Area
################################################################################
# Description: Plots of the study area

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(terra)        # To handle raster files
library(sf)           # To handle spatial data
library(tidyverse)    # For data wrangling
library(ggspatial)    # To add scale bars and north arrows
library(colorspace)   # To lighten or darken colors
library(RColorBrewer) # For colors
library(ggpubr)       # To arrange plots
library(magick)       # To plot jpg
library(cowplot)      # To plot jpg

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load stuff that we would like to plot
africa <- read_sf("03_Data/02_CleanData/Africa.gpkg")
prot   <- read_sf("03_Data/02_CleanData/Protected.gpkg")
water  <- read_sf("03_Data/02_CleanData/MajorWaters.gpkg")
river  <- read_sf("03_Data/02_CleanData/MajorRivers.gpkg")
roads  <- read_sf("03_Data/02_CleanData/Roads.gpkg")
areas  <- read_sf("03_Data/02_CleanData/Sources.gpkg")
vills  <- read_sf("03_Data/02_CleanData/Villages.gpkg")
vills  <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y) %>%
  mutate(place = gsub(place, pattern = "Village", replacement = "Villages")) %>%
  mutate(place = gsub(place, pattern = "City", replacement = "Cities"))

# Load dispersal trajectories and derive extent from them
# disp <- read_csv("03_Data/02_CleanData/DispersersSubsampled.csv")
disp <- "03_Data/02_CleanData/Dispersers.csv" %>%
  read_csv() %>%
  nest(GPS = -ID) %>%
  mutate(GPS = map(GPS, function(x) {
    dis <- subset(x, State == "Disperser" & !is.na(x))
    dis <- resampleFixes(dis, hours = 1, start = 1, tol = 0.5)
    return(dis)
  })) %>%
  unnest(GPS)
exte <- ext(c(22, 27, -20.7, -17.9))

# Load moving windows
wind <- read_rds("03_Data/02_CleanData/Windows.rds")
wind <- vect(wind$Window)
wind <- aggregate(wind, dissolve = T)
wind <- as.polygons(exte) - wind
crs(wind) <- "+init=epsg:4326"
wind <- st_as_sf(wind)

# Calculate a centerpoint
cent <- data.frame(x = mean(exte[1:2]), y = mean(exte[3:4]))

# Reorder the levels of national parks
prot$Desig <- factor(prot$Desig, levels = c("National Park", "Forest Reserve", "Protected"))

# Create country labels
labels_countries <- data.frame(
    x     = c(25.2, 26.4, 25.7, 21.5, 22.5)
  , y     = c(-19.3, -18.2, -17.6, -17.6, -18.0)
  , Label = c("BOTSWANA", "ZIMBABWE", "ZAMBIA", "ANGOLA", "NAMIBIA")
)

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 22.7, 25.3)
  , y     = c(-19, -18.1, -20.4, -20.5)
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
    x     = c(20.65, 23.18, 23.15, 23.9)
  , y     = c(-17.85, -18.65, -20.35, -20.25)
  , Label = c("Okavango", "Selinda", "Nhabe", "Boteti")
)

################################################################################
#### Africa
################################################################################
# Prepare plot of africa
p1 <- ggplot() +
  geom_sf(data = africa, fill = "gray90", col = "gray90", lwd = 0.5) +
  geom_point(data = cent, aes(x = x, y = y), col = "darkgreen", alpha = 0.4, size = 5) +
  theme_minimal() +
  coord_sf(xlim = c(-18, 48), ylim = c(-32, 35)) +
  theme(
      panel.grid = element_blank()
    , axis.ticks = element_blank()
    , axis.text  = element_blank()
    , axis.title = element_blank()
  )

################################################################################
#### Study Area
################################################################################
# Prepare plot of study area
p2 <- ggplot() +
  geom_sf(data = prot, aes(fill = Desig), col = NA, alpha = 0.7) +
  geom_sf(data = water, fill = "cornflowerblue", col = NA) +
  geom_sf(data = river, col = "cornflowerblue") +
  geom_sf(data = roads, col = "gray70", lwd = 0.2) +
  geom_sf(data = africa, col = "black", fill = NA, lwd = 0.3) +
  # geom_sf(data = wind, col = "red", fill = "red", alpha = 0.2) +
  geom_point(data = disp, mapping = aes(x = x, y = y), size = 0.1, col = "gray30") +
  geom_path(data = disp, mapping = aes(x = x, y = y, group = ID), linewidth = 0.1, col = "gray30") +
  # geom_sf(data = areas, fill = NA, col = "black", linewidth = 0.5) +
  # geom_sf_text(data = areas, aes(label = Name), col = "black") +
  geom_point(
      data        = subset(vills, place == "Cities")
    , mapping     = aes(x = x, y = y)
    , col         = "orange"
    , shape       = 17
    , show.legend = F
  ) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("cornflowerblue", 0.4)
    , fontface = 3
    , size     = 3
  ) +
  geom_text(
      data     = labels_rivers
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("cornflowerblue", 0.4)
    , fontface = 3
    , size     = 1.5
    , angle    = c(-30, 40, 20, -32)
  ) +
  geom_text(
      data     = subset(labels_nationalparks, Label != "Moremi\nNP")
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("green", 0.8)
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
      data     = labels_countries
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "black"
    , fontface = 2
  ) +
  geom_text(
      data     = subset(vills, place == "Cities")
    , mapping  = aes(x = x, y = y, label = name)
    , col      = "orange"
    , fontface = 2
    , size     = 3
    , nudge_y  = c(0.1, -0.1, 0.1)
  ) +
  scale_fill_brewer(palette = "Greens", direction = -1, name = "") +
  annotation_scale(
      location   = "br"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
    # , pad_x      = unit(0.9, "cm")
    # , pad_y      = unit(0.8, "cm")
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    # , pad_x    = unit(0.7, "cm")
    , pad_y    = unit(0.8, "cm")
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
    , xlim   = c(xmin(exte), xmax(exte))
    , ylim   = c(ymin(exte), ymax(exte))
    , expand = F
  )

################################################################################
#### Legends
################################################################################
# Create a new plot for the legend
df <- prot[1:4, ]
df <- dplyr::select(df, geom)
df <- cbind(
    df
  , Name   = factor(c("National Parks", "Forest Reserves", "Protected Areas", "Water")
  , levels = c("National Parks", "Forest Reserves", "Protected Areas", "Water"))
  , Color  = c(rev(brewer.pal(n = 3, name = "Greens")), "cornflowerblue")
  , Alpha  = c(0.7, 0.7, 0.7, 1)
)
df2 <- data.frame(
    x     = c(1, 2, 1)
  , y     = c(1, 2, 1)
  , Class = as.factor(c("Cutlines", "Faults", "Roads"))
)

# Plot
l1 <- ggplot() +
  geom_sf(data = df, aes(fill = Name, alpha = Name, col = Name)) +
  scale_fill_manual(values = df$Color, name = "Legend") +
  scale_color_manual(values = c("white", "white", "white", "white", "orange", "purple"), name = "Legend") +
  scale_alpha_manual(values = c(0.7, 0.7, 0.7, 1.0, 0.2, 0.2), name = "Legend") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
l2 <- ggplot(data = data.frame(x = c(0, 1), y = c(0, 1), ID = c("Dispersal Trajectories", "Dispersal Trajectories")), aes(x = x, y = y, col = ID)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("gray30"), name = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())

# Grab legends
l1 <- get_legend(l1)
l2 <- get_legend(l2)

################################################################################
#### Wet vs Dry Season
################################################################################
# Create plots for the wet and dry season
wet <- image_read("04_Manuscript/Figures/2023-03-10.jpg")
wet <- image_crop(wet, "800x1274") %>% image_border(color = "black")
dry <- image_read("04_Manuscript/Figures/2022-09-16.jpg")
dry <- image_crop(dry, "800x1274") %>% image_border(color = "black")

# Plot
p3 <- ggdraw() +
  draw_image(wet)
p4 <- ggdraw() +
  draw_image(dry)

################################################################################
#### Store Plots
################################################################################
ggsave("04_Manuscript/Figures/StudyArea1.png", bg = "white", plot = p1, width = 3, height = 3, scale = 2, device = png)
ggsave("04_Manuscript/Figures/StudyArea2.png", bg = "white", plot = p2, width = 5.2, height = 3, scale = 1.5, device = png)
ggsave("04_Manuscript/Figures/StudyArea3.png", bg = "transparent", plot = l1, width = 2.5, height = 0.3, scale = 2, device = png)
ggsave("04_Manuscript/Figures/StudyArea4.png", bg = "transparent", plot = l2, width = 1.0, height = 0.2, scale = 2, device = png)
