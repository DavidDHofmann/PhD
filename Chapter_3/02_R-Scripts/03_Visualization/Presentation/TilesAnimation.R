################################################################################
#### Animations of Tiles
################################################################################
# Description: Visualizing the downloaded Sentinel tiles

# Clear R's brain
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(terra)               # To handle rasters
library(lubridate)           # To handle dates
library(ggspatial)           # To add scale bar and north arrow to plots
library(tidyverse)           # For data wrangling
library(sf)                  # For plotting
library(animation)           # To store animation
library(pbmcapply)           # For multicore abilities with progress bar
library(magick)              # To load images
library(ggdark)              # To use dark themes

# Load data to plot in the background
r       <- rast("03_Data/02_CleanData/00_General_Raster.tif")
r       <- as.data.frame(r, xy = T)
area    <- read_sf("03_Data/02_CleanData/00_General_Shapefile.shp")
water   <- read_sf("03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters.shp")
africa  <- read_sf("03_Data/02_CleanData/00_General_Africa.shp")
windows <- read_rds("/media/david/Elements/Windows.rds")

# Loop through the data and plot the tiles that we downloaded
plotnames <- sapply(1:nrow(windows), function(i) {

  # Prepare the plot
  p <- ggplot() +
    geom_sf(data = africa, color = "gray50", fill = NA) +
    geom_sf(data = water, color = NA, fill = "cornflowerblue") +
    geom_point(data = windows$GPS[[i]], aes(x = x, y = y)) +
    geom_sf(data = st_as_sf(windows$Window[[i]])
      , fill  = NA
      , color = "white"
      , alpha = 0.25
      , lwd   = 0.5
    ) +
    geom_sf(data = st_as_sf(windows$Tiles[[i]])
      , fill  = NA
      , color = "orange"
      , alpha = 0.25
      , lwd   = 0.5
    ) +
    coord_sf(
        crs    = 4326
      , xlim   = c(min(r$x), max(r$x))
      , ylim   = c(min(r$y), max(r$y))
      , expand = F
    ) +
    ggtitle(paste(month.abb[windows$Month[[i]]], windows$Year[[i]])) +
    dark_theme_void() +
    theme(
        plot.title       = element_text(vjust = -10)
      , panel.background = element_blank()
      , panel.border     = element_blank()
      , plot.background  = element_blank()
    ) +
    annotation_scale(
        location   = "bl"
      , width_hint = 0.2
      , line_width = 0.5
      , text_cex   = 0.5
      , height     = unit(0.15, "cm")
      , bar_cols   = c("white", "black")
      , text_col   = "white"
    ) +
    annotation_north_arrow(
        location = "br"
      , height   = unit(1.5, "cm"),
      , width    = unit(1.2, "cm"),
      , style    = north_arrow_fancy_orienteering(
            fill      = c("white", "white")
          , line_col  = NA
          , text_col  = "white"
          , text_size = 12
        )
    )

  # Store it to a file and return the filename
  filename <- tempfile(fileext = ".png")
  suppressWarnings(
    ggsave(
        plot     = p
      , filename = filename
      , width    = 2000
      , height   = 1800
      , units    = "px"
      , bg       = "transparent"
    )
  )
  return(filename)

})

# Create animation
ani.options(interval = 0.5, ani.width = 2000, ani.height = 1800)
im.convert(plotnames, "05_Presentation/Tiles.gif")
