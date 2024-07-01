################################################################################
#### Egression
################################################################################
# Description: Visualization of egressions into different egression zones

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(raster)         # To handle spatial data
library(tidyverse)      # To wrangle data
library(sf)             # To plot spatial features
library(geomtextpath)   # For angled labels in ggplot
library(ggspatial)      # To add scale bars and north arrows
library(ggpubr)         # To extract legends
library(colorspace)     # To lighten or darken colors

# Reload data of interpatch connectivity and subset to entries relating to
# egression
egre <- "03_Data/03_Results/InterpatchConnectivityBootstrapped.rds" %>%
  read_rds() %>%
  subset(Type == "Egression") %>%
  subset(SourceArea == "Overall" & CurrentArea != "Overall") %>%
  mutate(DispersalSuccess = round(DispersalSuccess)) %>%
  mutate(CurrentArea = factor(CurrentArea, levels = 7:14)) %>%
  mutate(FloodLevel = factor(FloodLevel, levels = c("Min", "Max")))

# Load stuff that we would like to plot
water  <- read_sf("03_Data/02_CleanData/MajorWaters.shp")
river  <- read_sf("03_Data/02_CleanData/MajorRivers.shp")
areas  <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
cutli  <- read_sf("03_Data/02_CleanData/Cutlines.shp")

# I also want to add the Okavango river
kava <- raster("03_Data/01_RawData/MERIT/Rivers.tif")
kava <- crop(kava, extent(c(20.5, 22, -18.2, -17.5)))
kava <- rasterToPolygons(kava, fun = function(x) {x == 1})
kava <- aggregate(kava)
kava <- buffer(kava, width = 0.1 / 111, dissolve = T)
kava <- st_as_sf(kava)

# Generate area labels
labels_areas <- st_coordinates(st_point_on_surface(areas))
labels_areas <- cbind(labels_areas, st_drop_geometry(areas))

# Prepare the plot
p1 <- ggplot(egre, mapping = aes(
        x    = as.factor(CurrentArea)
      , y    = DispersalSuccess
      , fill = FloodLevel
      , col  = FloodLevel
    )
  ) +
  geom_bar(stat = "identity"
    , position = position_dodge(width = 1)
    , size     = 0.2
    , width    = 0.75
    , alpha    = 0.8
  ) +
  geom_errorbar(aes(ymin = DispersalSuccess - DispersalSuccessSD, ymax = DispersalSuccess + DispersalSuccessSD)
    , position = position_dodge(width = 1)
    , width    = 0.2
  ) +
  geom_textpath(aes(y = DispersalSuccess + 200, label = format(DispersalSuccess, big.mark = ",", scientific = FALSE))
    , position = position_dodge(width = 1)
    , col      = "black"
    , size     = 2
    , angle    = 90
  ) +
  coord_polar(direction = -1, start = 2.35) +
  theme_void() +
  scale_fill_manual(
      values = c("orange", "cornflowerblue")
    , name   = "Flood-Level"
  ) +
  scale_color_manual(
      values = lighten(c("orange", "cornflowerblue"), amount = 0.5)
    , name   = "Flood-Level"
  ) +
  theme(legend.position = "none")

# Overlay the plot above with a plot of the emigration zones
scale <- 1.5
p2 <- ggplot() +
  geom_sf(data = water, fill = "cornflowerblue", col = NA, alpha = 0.2) +
  geom_sf(data = river, col = "cornflowerblue", alpha = 0.2) +
  geom_sf(data = cutli, lty = 2, size = 0.3) +
  geom_sf(
      data  = subset(areas, Type == "Buffer")
    , col   = "transparent"
    , fill  = "purple"
    , alpha = 0.2
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
  theme(
      legend.position  = "none"
    , panel.border     = element_rect(colour = "gray30", fill = NA, size = 1)
    , panel.background = element_blank()
    , strip.background = element_rect(fill = "gray95", color = "transparent")
  ) +
  xlab("") +
  ylab("") +
  coord_sf(
      crs    = 4326
    , xlim   = c(21, 24.5)
    , ylim   = c(-17.8, -20.8)
    , expand = F
  ) +
  annotation_custom(
      ggplotGrob(p1)
    , xmin = 22.8 - scale
    , xmax = 22.8 + scale
    , ymin = -19.26 - scale
    , ymax = -19.26 + scale
  )

# Create a separate legend
l <- ggplot(egre, aes(x = as.factor(CurrentArea), y = DispersalSuccess, fill = FloodLevel, col = FloodLevel)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), size = 0.3, width = 0.75, alpha = 0.8) +
  scale_fill_manual(values = c("orange", "cornflowerblue"), name = "Flood-Level") +
  scale_color_manual(values = c("orange", "cornflowerblue"), name = "Flood-Level") +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  guides(color = guide_legend(
          title.position = "top"
        , title.hjust    = 0.5
        , label.position = "bottom"
        , label.hjust    = 0.5
        , keywidth       = 3
      )
    ) +
  guides(fill = guide_legend(
          title.position = "top"
        , title.hjust    = 0.5
        , label.position = "bottom"
        , label.hjust    = 0.5
        , keywidth       = 3
      )
    )
l <- get_legend(l)

# Put them together
p <- ggarrange(p2, l, heights = c(0.9, 0.15), nrow = 2)

# Store the plot to file
ggsave("04_Manuscript/Figures/Egression.png"
  , plot   = p
  , bg     = "white"
  , height = 7
  , width  = 7
  , scale  = 0.8
  , device = png
)
