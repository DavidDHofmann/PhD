################################################################################
#### Covariates
################################################################################
# Description: Plot of the covariates

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(tidyverse)     # For data wrangling
library(sf)            # For plotting
library(ggspatial)     # For north arrow and scale bar

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load hyena gps data and covariates
gps <- read_rds("03_Data/Hyena_GPS.rds")
cov <- stack("03_Data/Hyena_Covariates.tif")

################################################################################
#### Plot of Covariates and Observed Trajectories
################################################################################
# Crop to extent
ext <- extent(c(range(gps$x), range(gps$y))) + 20000

# Prepare the data
dat_cov <- cov %>%
  crop(., ext) %>%
  as.data.frame(xy = T) %>%
  gather(key = covariate, value = value, 3:5) %>%
  group_by(covariate) %>%
  mutate(value = (value - min(value, na.rm = T)) / (max(value, na.rm = T) - min(value, na.rm = T))) %>%
  mutate(covariate = factor(covariate, c("Water", "DistanceToWater", "Trees")))
dat_mov <- gps %>%
  subset(DOP < 10)

# Plot
p <- ggplot(dat_cov, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_path(data = dat_mov, inherit.aes = F, aes(x = x, y = y), linewidth = 0.1) +
    scale_fill_viridis_c(
        option = "viridis"
      , name   = "Value"
      , breaks = c(0, 0.5, 1)
      , labels = c("Low", "Medium", "High")
    ) +
    coord_sf() +
    theme_minimal() +
    facet_wrap(~ covariate) +
    annotation_scale(
        location   = "bl"
      , width_hint = 0.3
      , line_width = 0.5
      , pad_y      = unit(0.40, "cm")
      , pad_x      = unit(0.40, "cm")
      , height     = unit(0.15, "cm")
      , bar_cols   = "white"
      , text_col   = "white"
    ) +
    annotation_north_arrow(
        location = "br"
      , height   = unit(1.0, "cm"),
      , width    = unit(0.8, "cm"),
      , style    = north_arrow_fancy_orienteering(
            fill      = c("white", "white")
          , line_col  = NA
          , text_col  = "white"
          , text_size = 8
        )
    ) +
    xlab("Easting") +
    ylab("Northing") +
    theme(
        axis.title.y     = element_text(angle = 90, vjust = 0.5)
      , axis.text.x      = element_text(angle = 45, vjust = 0.5)
      , strip.background = element_rect(fill = "gray95", color = "white")
    )

# Store to file
ggsave("04_Manuscript/99_CaseStudyCovariates.png"
  , plot   = p
  , width  = 6
  , height = 2.5
  , bg     = "white"
  , scale  = 1.3
  , device = png
)
