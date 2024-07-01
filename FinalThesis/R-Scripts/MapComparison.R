################################################################################
#### Comparison of Connectivity Maps
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(terra)        # To handle spatial data
library(raster)       # To handle spatial data
library(tidyverse)    # To wrangle data
library(RColorBrewer) # Color palettes
library(sf)           # To plot spatial features
library(ggspatial)    # For scale bar and north arrow
library(scales)       # To squish colorscale

# Function to normalize a raster
normalizeRaster <- function(x) {
  (x - min(x[], na.rm = T)) / (max(x[], na.rm = T) - min(x[], na.rm = T))
}

# Load connectivity maps from previous chapters
map_ch0 <- "/home/david/ownCloud/University/15. PhD/Chapter_0/03_Data/03_Results/99_LeastCostCorridors1.tif" %>%
    rast() %>%
    setNames("Least-Cost Corridors")

# Prepare map from simulation chapter
map_ch1 <- "/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/03_Results/99_Heatmaps.rds" %>%
    read_rds() %>%
    subset(steps == 2000) %>%
    pull(heatmap) %>%
    stack() %>%
    sum() %>%
    rast() %>%
    normalizeRaster() %>%
    setNames("Simulated Heatmap")

# Prepare maps from flood extreme chapter
map_ch2 <- "/home/david/ownCloud/University/15. PhD/Chapter_2/03_Data/03_Results/HeatmapsBetweennessGlobal.rds" %>%
    read_rds() %>%
    subset(Steps == 2000) %>%
    pull(Heatmap) %>%
    stack() %>%
    rast() %>%
    setNames(c("Minimum Flood", "Maximum Flood")) %>%
    normalizeRaster()

# Prepare maps from seasonality chapter
map_ch3 <- "/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/03_Results/Connectivity.rds" %>%
    read_rds() %>%
    mutate(Filename = file.path("/home/david/ownCloud/University/15. PhD/Chapter_3", Filename)) %>%
    mutate(Data = map(Filename, function(x) {
      heatmap <- read_rds(x)
      heatmap <- heatmap$Heatmap
      heatmap <- unwrap(heatmap)
      return(heatmap)
    })) %>%
    subset(Formula == "Full") %>%
    dplyr::select(Source, ModelCode, Heatmap = Data) %>%
    nest(Heatmaps = -ModelCode) %>%
    mutate(Heatmaps = map(Heatmaps, function(x) {
        sum(rast(x$Heatmap))
    })) %>%
    pull(Heatmaps) %>%
    do.call(c, .) %>%
    disagg(fact = 4, method = "bilinear") %>%
    focal(w = matrix(1, 5, 5), fun = "mean") %>%
    aggregate(fact = 2, fun = "mean") %>%
    setNames(c("Static Configuration", "Dynamic Configuration")) %>%
    normalizeRaster()

# Define extent
exte <- ext(21, 26, -21, -17.5)

# Prepare reference raster
r <- rast(exte, resolution = 1000 / 111000, crs = "epsg:4326")

# Reproject
map_ch0 <- resample(map_ch0, r, method = "bilinear")
map_ch1 <- resample(map_ch1, r, method = "bilinear")
map_ch2 <- resample(map_ch2, r, method = "bilinear")
map_ch3 <- resample(map_ch3, r, method = "bilinear")

# Crop connectivity maps accordingly
map_ch0 <- crop(map_ch0, exte)
map_ch1 <- crop(map_ch1, exte)
map_ch2 <- crop(map_ch2, exte)
map_ch3 <- crop(map_ch3, exte)

# Put them all into a stack
maps <- c(map_ch0, map_ch1, map_ch2, map_ch3) %>%
  as.data.frame(xy = T) %>%
  pivot_longer(3:ncol(.), names_to = "Map", values_to = "Value") %>%
  mutate(Map = factor(Map, levels = c("Least-Cost Corridors", "Simulated Heatmap", "Minimum Flood", "Maximum Flood", "Static Configuration", "Dynamic Configuration")))

# Load source areas
src_ch0 <- "/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/02_CleanData/Protected.gpkg" %>%
    read_sf() %>%
    dplyr::select(Name, geom) %>%
    mutate(Map = "Least-Cost Corridors")

src_ch1 <- "/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/02_CleanData/Protected.gpkg" %>%
    read_sf() %>%
    dplyr::select(Name, geom) %>%
    mutate(Map = "Simulated Heatmap")

src_ch2 <- "/home/david/ownCloud/University/15. PhD/Chapter_2/03_Data/02_CleanData/SourceAreas.shp" %>%
    read_sf() %>%
    subset(Type == "Main") %>%
    dplyr::select(Name = ID, geom = geometry) %>%
    expand_grid(., Map = c("Minimum Flood", "Maximum Flood")) %>%
    st_as_sf()

src_ch3 <- "/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/02_CleanData/Sources.gpkg" %>%
    read_sf() %>%
    dplyr::select(Name, geom) %>%
    expand_grid(., Map = c("Static Configuration", "Dynamic Configuration")) %>%
    st_as_sf()

# Put them together
src <- rbind(src_ch0, src_ch1, src_ch2, src_ch3) %>%
  mutate(Map = factor(Map, levels = c("Least-Cost Corridors", "Simulated Heatmap", "Minimum Flood", "Maximum Flood", "Static Configuration", "Dynamic Configuration")))


# Load additional stuff to plot
roads <- read_sf("/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/02_CleanData/Roads.gpkg")
afric <- read_sf("/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/02_CleanData/Africa.gpkg")
vills <- read_sf("/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/02_CleanData/Villages.gpkg")
vills <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y)

# Create country labels
labels_countries <- data.frame(
    x     = c(24.5, 25.5, 25.7, 21.5, 23.5)
  , y     = c(-19.3, -18.2, -17.6, -17.6, -17.8)
  , Label = c("BOTSWANA", "ZIMBABWE", "ZAMBIA", "ANGOLA", "NAMIBIA")
)

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.1, 25.2)
  , y     = c(-19.1, -18.2, -17.5, -20.7)
  , Label = c("Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans")
)

# Prepare plot
spectral <- colorRampPalette(rev(brewer.pal(11, name = "Spectral")))
colorscheme <- "dark"
p <- ggplot() +
  geom_raster(data = maps, aes(x = x, y = y, fill = Value)) +
  geom_sf(
      data        = src
    , color       = "white"
    , fill        = "white"
    , lty         = 1
    , linewidth   = 0.2
    , show.legend = F
    , alpha       = 0.15
  ) +
  geom_sf(
      data      = roads
    , col       = ifelse(colorscheme == "dark", "gray90", "gray30")
    , linewidth = 0.1
  ) +
  geom_point(
      data        = subset(vills, place == "City")
    , mapping     = aes(x = x, y = y, size = place)
    , col         = ifelse(colorscheme == "dark", "gray90", "gray30")
    , shape       = 15
    , show.legend = F
    , size        = 1
  ) +
  geom_sf(
      data      = afric
    , linewidth = 0.4
    , col       = "black"
    , fill      = NA
  ) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = ifelse(colorscheme == "dark", "gray80", "gray30")
    , fontface = 3
    , size     = 1.3
  ) +
  geom_text(
      data     = labels_countries
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "black"
    , size     = 2
    , fontface = 2
  ) +
  geom_text(
      data     = subset(vills, place == "City")
    , mapping  = aes(x = x, y = y, label = name)
    , col      = ifelse(colorscheme == "dark", "gray90", "gray30")
    , fontface = 3
    , size     = 2
    , nudge_y  = c(0.1, -0.1, 0.1)
  ) +
  facet_wrap(~ Map, dir = "v", nrow = 2) +
  scale_x_continuous(breaks = seq(21.5, 26.5, by = 1)) +
  scale_y_continuous(breaks = seq(-20.5, -17.5, by = 1)) +
  coord_sf(
      crs    = 4326
    , xlim   = exte[1:2]
    , ylim   = exte[3:4]
    , expand = F
  ) +
  theme_minimal() +
  scale_fill_gradientn(
      colors  = spectral(100)
    , breaks  = c(0, 0.45, 0.9)
    , labels  = c("Low", "Medium", "High")
    , oob     = squish
    , guide   = guide_colorbar(
      , title          = "Connectivity"
      , show.limits    = T
      , title.position = "bottom"
      , title.hjust    = 0.5
      , ticks          = F
      , barheight      = unit(0.2, "cm")
      , barwidth       = unit(4, "cm")
    )
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
  ) +
  theme(
      legend.position  = "bottom"
    , legend.box       = "vertical"
    , panel.background = element_blank()
    , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
    , strip.background = element_rect(fill = "gray95", color = "transparent")
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , text_cex   = 0.5
    , height     = unit(0.1, "cm")
    , bar_cols   = c("white", "black")
    , text_col   = ifelse(colorscheme == "dark", "white", "black")
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(0.7, "cm"),
    , width    = unit(0.6, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = ifelse(colorscheme == "dark", "white", "black")
        , line_col  = NA
        , text_col  = ifelse(colorscheme == "dark", "white", "black")
        , text_size = 4
      )
  )

# Store plot to file
ggsave("/home/david/ownCloud/University/15. PhD/FinalThesis/Chapter_06/Figures/MapComparison.png"
  , plot   = p
  , bg     = "white"
  , height = 3.8
  , width  = 6
  , scale  = 1.6
  , device = png
)
