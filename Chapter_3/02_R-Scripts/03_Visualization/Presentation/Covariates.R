################################################################################
#### Plot of Covariates
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(terra)
library(tidyverse)
library(ggspatial)

# Crop to extent of delta
ext <- ext(21.74909, 24.30089, -20.34901, -18.14901)

# Load covariates that we want to visualize and crop them
water <- rast("03_Data/02_CleanData/WaterStatic.tif") %>% crop(ext) %>% as.data.frame(xy = T) %>% setNames(c("x", "y", "layer"))
dista <- rast("03_Data/02_CleanData/DistanceToWaterStatic.tif") %>% crop(ext) %>% as.data.frame(xy = T) %>% setNames(c("x", "y", "layer"))
trees <- rast("/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/02_CleanData/00_Vegmaps/2022-03-06.tif")[[1]] %>% crop(ext) %>% as.data.frame(xy = T) %>% setNames(c("x", "y", "layer"))
shrub <- rast("/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/02_CleanData/00_Vegmaps/2022-03-06.tif")[[2]] %>% crop(ext) %>% as.data.frame(xy = T) %>% setNames(c("x", "y", "layer"))
human <- rast("03_Data/02_CleanData/Humans.tif") %>% crop(ext) %>% as.data.frame(xy = T) %>% setNames(c("x", "y", "layer"))

# Generate plots
p1 <- ggplot() +
  geom_raster(data = water, aes(x = x, y = y, fill = as.factor(layer))) +
  scale_fill_manual(values = c("white", "cornflowerblue")) +
  coord_sf(xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), crs = 4326, expand = F) +
  theme_void() +
  theme(legend.position = "none") +
  annotation_scale(
      location   = "bl"
    , bar_cols   = c("white", "black")
    , text_col   = "black"
  ) +
  annotation_north_arrow(
      location = "br"
    , style    = north_arrow_fancy_orienteering(
          fill      = "black"
        , line_col  = NA
        , text_col  = "black"
      )
  )

p2 <- ggplot() +
  geom_raster(data = dista, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(direction = 1, trans = "sqrt") +
  coord_sf(xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), crs = 4326, expand = F) +
  theme_void() +
  theme(legend.position = "none") +
  annotation_scale(
      location   = "bl"
    , bar_cols   = c("white", "black")
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , style    = north_arrow_fancy_orienteering(
          fill      = "white"
        , line_col  = NA
        , text_col  = "white"
      )
  )

p3 <- ggplot() +
  geom_raster(data = trees, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(direction = 1, trans = "sqrt") +
  coord_sf(xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), crs = 4326, expand = F) +
  theme_void() +
  theme(legend.position = "none") +
  annotation_scale(
      location   = "bl"
    , bar_cols   = c("white", "black")
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , style    = north_arrow_fancy_orienteering(
          fill      = "white"
        , line_col  = NA
        , text_col  = "white"
      )
  )

p4 <- ggplot() +
  geom_raster(data = shrub, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(direction = 1, trans = "sqrt") +
  coord_sf(xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), crs = 4326, expand = F) +
  theme_void() +
  theme(legend.position = "none") +
  annotation_scale(
      location   = "bl"
    , bar_cols   = c("white", "black")
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , style    = north_arrow_fancy_orienteering(
          fill      = "white"
        , line_col  = NA
        , text_col  = "white"
      )
  )

p5 <- ggplot() +
  geom_raster(data = human, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(direction = 1, trans = "sqrt") +
  coord_sf(xlim = c(xmin(ext), xmax(ext)), ylim = c(ymin(ext), ymax(ext)), crs = 4326, expand = F) +
  theme_void() +
  theme(legend.position = "none") +
  annotation_scale(
      location   = "bl"
    , bar_cols   = c("white", "black")
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , style    = north_arrow_fancy_orienteering(
          fill      = "white"
        , line_col  = NA
        , text_col  = "white"
      )
  )

# Store the plots
ggsave("05_Presentation/Water.png", plot = p1, device = png)
ggsave("05_Presentation/DistanceToWater.png", plot = p2, device = png)
ggsave("05_Presentation/Woodland.png", plot = p3, device = png)
ggsave("05_Presentation/Shrubs.png", plot = p4, device = png)
ggsave("05_Presentation/Humans.png", plot = p5, device = png)
