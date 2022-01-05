################################################################################
#### Animations of Dynamic Covariates
################################################################################
# Description: Preparation of animations of all dynamic covariates

# Clear R's brain
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(terra)               # To handle rasters
library(RColorBrewer)        # To get color scales
library(lubridate)           # To handle dates
library(ggspatial)           # To add scale bar and north arrow to plots
library(tidyverse)           # For data wrangling
library(sf)                  # For plotting
library(animation)           # To store animation
library(pbmcapply)           # For multicore abilities with progress bar
library(magick)              # To load images

# Load data to plot in the background
r      <- rast("03_Data/02_CleanData/00_General_Raster.tif")
africa <- vect("03_Data/02_CleanData/00_General_Africa.shp") %>% crop(r)

# Coerce data to formats that are plottable with ggplot
r <- as.data.frame(r, xy = T)
africa <- st_as_sf(africa)

################################################################################
#### Temperature
################################################################################
# Load data
temp <- dir(path = "03_Data/02_CleanData/00_Tempmaps", pattern = ".grd$", full.names = T) %>%
  rast()

# Identify dates
dates <- names(temp) %>%
  substr(start = 1, stop = 10) %>%
  ymd()

# Keep only maps from 2019
temp <- temp[[year(dates) == 2019]]

# Prepare frames and store them to file
plot_names <- paste0(tempfile(), "_", 1:nlyr(temp), ".png")
invisible(pbmclapply(1:nlyr(temp), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {

  # Get the date and time and prepare a title from those
  datetime <- names(temp[[x]])
  date <- substr(datetime, start = 1, stop = 10)
  time <- substr(datetime, start = 12, stop = 13)
  datetime <- make_datetime(year(date), month(date), day(date), time, 0, 0)
  datetime <- format(datetime, "%Y-%m-%d %H:%M:%S")

  # Convert the raster to a dataframe
  map <- disagg(temp[[x]], fact = 4, method = "bilinear")
  map <- as.data.frame(map, xy = T)
  names(map)[3] <- "value"

  # Make the plot
  p <- ggplot(map, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = africa, inherit.aes = F, fill = NA, col = "gray30") +
    geom_sf_text(data = africa, inherit.aes = F, aes(label = COUNTRY), nudge_y = 0.1, col = "gray30") +
    coord_sf(
        crs    = 4326
      , xlim   = c(min(r$x), max(r$x))
      , ylim   = c(min(r$y), max(r$y))
      , expand = F
    ) +
    scale_fill_gradientn(
        colors = brewer.pal(11, name = "RdYlBu")[11:1]
      , labels = function(x){paste0(x, "°")}
      , limits = c(0, 40)
      , guide  = guide_colorbar(
        , title          = expression("Temperature (C)")
        , show.limits    = T
        , title.position = "top"
        , title.hjust    = 0.5
        , ticks          = T
        , barheight      = unit(0.6, "cm")
        , barwidth       = unit(10.0, "cm")
      )
    ) +
    labs(
        x        = NULL
      , y        = NULL
      , fill     = NULL
      , title    = "Temperature"
      , subtitle = datetime
    ) +
    theme(
        legend.position  = "bottom"
      , legend.box       = "vertical"
      , panel.background = element_blank()
      , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    annotation_scale(
        location   = "bl"
      , width_hint = 0.2
      , line_width = 0.5
      , height     = unit(0.15, "cm")
    ) +
    annotation_north_arrow(
        location = "br"
      , height   = unit(1.5, "cm"),
      , width    = unit(1.2, "cm"),
      , style    = north_arrow_fancy_orienteering(
            fill      = c("black", "black")
          , line_col  = NA
          , text_col  = "black"
          , text_size = 12
        )
    )

  # Store the plot to file
  suppressWarnings(
    ggsave(
        plot     = p
      , filename = plot_names[x]
      , width    = 2000
      , height   = 1800
      , units    = "px"
    )
  )

}))

# Animate the plots
ani.options(interval = 1 / 20, ani.width = 2000, ani.height = 1800)
saveVideo({
  for (i in 1:length(plot_names)){
    img <- image_read(plot_names[[i]])
    plot(img)
  }
}, video.name = "Temperature.mp4")

################################################################################
#### Continue Here
################################################################################

################################################################################
#### Precipitation
################################################################################
# Load data
prec <- dir(path = "03_Data/02_CleanData/00_Rainmaps", pattern = ".grd$", full.names = T) %>%
  rast()

# Identify dates
dates <- names(prec) %>%
  substr(start = 1, stop = 10) %>%
  ymd()

# Keep only maps from 2019
prec <- prec[[year(dates) == 2019]]

# Prepare frames
plots_prec <- pbmclapply(1:nlyr(prec), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {

  # Get the date and time and prepare a title from those
  datetime <- names(prec[[x]])
  date <- substr(datetime, start = 1, stop = 10)
  time <- substr(datetime, start = 12, stop = 13)
  datetime <- make_datetime(year(date), month(date), day(date), time, 0, 0)
  datetime <- format(datetime, "%Y-%m-%d %H:%M:%S")

  # Convert the raster to a dataframe
  map <- disagg(prec[[x]], fact = 4, method = "bilinear")
  map <- as.data.frame(map, xy = T)
  names(map)[3] <- "value"

  # Make the plot
  p <- ggplot(map, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = africa, inherit.aes = F, fill = NA, col = "gray70") +
    geom_sf_text(data = africa, inherit.aes = F, aes(label = COUNTRY), nudge_y = 0.1, col = "gray70") +
    coord_sf(
        crs    = 4326
      , xlim   = c(min(r$x), max(r$x))
      , ylim   = c(min(r$y), max(r$y))
      , expand = F
    ) +
    scale_fill_gradientn(
        # colors = c("white", brewer.pal(11, name = "RdBu")[11:1])
        colors = c("white", rev(brewer.pal(11, name = "RdYlBu")))
      , limits = c(0, 10)
      , oob    = squish
      , guide  = guide_colorbar(
        , title          = expression("Precipitation (mm)")
        , show.limits    = T
        , title.position = "top"
        , title.hjust    = 0.5
        , ticks          = T
        , barheight      = unit(0.6, "cm")
        , barwidth       = unit(10.0, "cm")
      )
    ) +
    labs(
        x        = NULL
      , y        = NULL
      , fill     = NULL
      , title    = "Precipitation"
      , subtitle = datetime
    ) +
    theme(
        legend.position  = "bottom"
      , legend.box       = "vertical"
      , panel.background = element_blank()
      , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    annotation_scale(
        location   = "bl"
      , width_hint = 0.2
      , line_width = 0.5
      , height     = unit(0.15, "cm")
    ) +
    annotation_north_arrow(
        location = "br"
      , height   = unit(1.5, "cm"),
      , width    = unit(1.2, "cm"),
      , style    = north_arrow_fancy_orienteering(
            fill      = c("black", "black")
          , line_col  = NA
          , text_col  = "black"
          , text_size = 12
        )
    )

  # Return the plot
  return(p)
})

# Show one of the plots
plots_prec[[10]]

# Store the plots to temporary files
plot_names <- paste0(tempfile(), "_", 1:length(plots_prec), ".png")
pbmclapply(1:length(plots_prec), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  suppressWarnings(
    ggsave(
        plot     = plots_prec[[x]]
      , filename = plot_names[x]
      , width    = 2000
      , height   = 1800
      , units    = "px"
    )
  )
})

# Animate the plots
ani.options(interval = 1 / 20, ani.width = 2000, ani.height = 1800)
saveVideo({
  for (i in 1:length(plot_names)){
    img <- image_read(plot_names[[i]])
    plot(img)
  }
}, video.name = "Precipitation.mp4")

################################################################################
#### NDVI
################################################################################
# Load data
ndvi <- rast("03_Data/02_CleanData/01_LandCover_NDVI.grd")

# Prepare frames
plots_ndvi <- pbmclapply(1:nlyr(ndvi), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {

  # Get the map date
  date <- ymd(names(ndvi[[x]]))

  # Make resolution a bit coarser for quicker plotting
  map <- aggregate(ndvi[[x]], fact = 2000 / (res(ndvi)[1] * 111000), fun = "mean")

  # Convert raster to a regular dataframe
  map <- as.data.frame(map, xy = T)
  names(map)[3] <- "value"

  # Make the plot
  p <- ggplot(map, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = africa, inherit.aes = F, fill = NA, col = "gray30") +
    geom_sf_text(data = africa, inherit.aes = F, aes(label = COUNTRY), nudge_y = 0.1, col = "gray30") +
    coord_sf(
        crs    = 4326
      , xlim   = c(min(r$x), max(r$x))
      , ylim   = c(min(r$y), max(r$y))
      , expand = F
    ) +
    scale_fill_gradientn(
        colors = brewer.pal(11, name = "RdYlGn")[4: 11]
      , limits = c(-0.5, 1)
      , guide  = guide_colorbar(
        , title          = expression("NDVI-Value")
        , show.limits    = T
        , title.position = "top"
        , title.hjust    = 0.5
        , ticks          = T
        , barheight      = unit(0.6, "cm")
        , barwidth       = unit(10.0, "cm")
      )
    ) +
    labs(
        x        = NULL
      , y        = NULL
      , fill     = NULL
      , title    = "NDVI"
      , subtitle = date
    ) +
    theme(
        legend.position  = "bottom"
      , legend.box       = "vertical"
      , panel.background = element_blank()
      , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    annotation_scale(
        location   = "bl"
      , width_hint = 0.2
      , line_width = 0.5
      , height     = unit(0.15, "cm")
    ) +
    annotation_north_arrow(
        location = "br"
      , height   = unit(1.5, "cm"),
      , width    = unit(1.2, "cm"),
      , style    = north_arrow_fancy_orienteering(
            fill      = c("black", "black")
          , line_col  = NA
          , text_col  = "black"
          , text_size = 12
        )
    )

  # Return the plot
  return(p)
})

# Store the plots to temporary files
plot_names <- paste0(tempfile(), "_", 1:length(plots_ndvi), ".png")
pbmclapply(1:length(plots_ndvi), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  suppressWarnings(
    ggsave(
        plot     = plots_ndvi[[x]]
      , filename = plot_names[x]
      , width    = 2000
      , height   = 1800
      , units    = "px"
    )
  )
})

# Animate the plots
ani.options(interval = 1 / 10, ani.width = 2000, ani.height = 1800)
saveVideo({
  for (i in 1:length(plot_names)){
    img <- image_read(plot_names[[i]])
    plot(img)
  }
}, video.name = "NDVI.mp4")
