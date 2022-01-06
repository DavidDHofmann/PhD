################################################################################
#### Plots of the Study Area
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(sf)           # To handle spatial data
library(ggplot2)      # For nice plots
library(ggspatial)    # To add scale bars and north arrows

# Load stuff that we would like to plot
africa <- read_sf("03_Data/02_CleanData/00_General_Africa.shp")
dogs   <- read_sf("03_Data/02_CleanData/00_General_WildDogs.shp")
area   <- read_sf("03_Data/02_CleanData/00_General_Shapefile.shp")
prot   <- read_sf("03_Data/02_CleanData/02_LandUse_ProtectedAreas.shp")
water  <- read_sf("03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters.shp")

# Reorder the levels of national parks
prot$Desig <- factor(prot$Desig, levels = c("National Park", "Forest Reserve", "Protected"))

# Prepare plot of africa
p1 <- ggplot() +
  geom_sf(data = africa, fill = "gray90", col = "white", lwd = 0.1) +
  geom_sf(data = dogs, fill = "gray80", col = NA) +
  geom_sf(data = area, fill = "cornflowerblue", col = "cornflowerblue", lwd = 0.5, alpha = 0.5, lty = 2) +
  theme_minimal() +
  theme(
      panel.grid = element_blank()
    , axis.ticks = element_blank()
    , axis.text  = element_blank()
  )

# Prepare plot of study area
p2 <- ggplot() +
  geom_sf(data = prot, aes(fill = Desig), col = NA, alpha = 0.7) +
  geom_sf(data = water, fill = "cornflowerblue", col = NA) +
  geom_sf(data = area, fill = NA) +
  scale_fill_brewer(palette = "Greens", direction = -1, name = "Protection Status") +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
    , pad_x      = unit(0.7, "cm")
    , pad_y      = unit(0.6, "cm")
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
  theme_minimal()

# Show the plots
p1
p2
