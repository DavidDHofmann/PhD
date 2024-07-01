################################################################################
#### Animations of Dynamic Covariates
################################################################################
# Description: Preparation of animations of all dynamic covariates

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
library(scales)              # To handle out of bounds (oob) values
library(ggpubr)              # To grab legends

# Load data to plot in the background
r      <- rast("03_Data/02_CleanData/00_General_Raster.tif")
africa <- vect("03_Data/02_CleanData/00_General_Africa.shp") %>% crop(r)

# Coerce data to formats that are plottable with ggplot
r <- as.data.frame(r, xy = T)
africa <- st_as_sf(africa)

################################################################################
#### Precipitation, Temperature, NDVI
################################################################################
# Loop through the different covariates and prepare the plots
covs <- c("Rainmaps", "Tempmaps", "NDVI")
for (i in covs) {
  # i <- covs[1]
  filename <- paste0("05_Presentation/", i, ".gif")
  if (file.exists(filename)) {
    next
  }

  # Load the data
  dat <- dir(
      path    = paste0("03_Data/02_CleanData/00_", i)
    , pattern = ".grd$", full.names = T
  ) %>% rast()

  # Retrieve dates from layernames
  dates <- names(dat) %>%
    substr(start = 1, stop = 10) %>%
    ymd()

  # Subset to a selected date range
  if (i == "NDVI") {
    subdat <- dat
  } else {
    subdat <- dat[[dates > "2019-10-24" & dates < "2019-11-10"]]
  }

  # Specify colors
  cols <- case_when(
      i == "Rainmaps" ~ c(hcl.colors(n = 50, palette = "RdYlBu"), "gray20")
    , i == "Tempmaps" ~ c(hcl.colors(n = 50, palette = "Inferno", rev = T), "transparent")
    , i == "NDVI" ~ c(hcl.colors(n = 50, palette = "Greens"), "transparent")
  )

  # Define a minimum and maximum values
  minmax <- case_when(
      i == "Rainmaps" ~ c(0, 5)
    , i == "Tempmaps" ~ c(0, 35)
    , i == "NDVI" ~ c(-1, 1)
  )

  # Plot labels
  labs <- case_when(
      i == "Rainmaps" ~ "Precipitation (mm)"
    , i == "Tempmaps" ~ "Temperature (Degrees C)"
    , i == "NDVI" ~ "NDVI"
  )

  # Prepare the frames
  cat("Preparing plots for", tolower(i), "\n")
  plot_names <- paste0(tempfile(), "_", 1:nlyr(subdat), ".png")
  invisible(pbmclapply(1:nlyr(subdat), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
    # x <- 1

    # Get the date and time and prepare a title from those
    datetime <- names(subdat[[x]])
    if (i != "NDVI") {
      date <- substr(datetime, start = 1, stop = 10)
      time <- substr(datetime, start = 12, stop = 13)
      datetime <- make_datetime(year(date), month(date), day(date), as.numeric(time), 0, 0)
      datetime <- format(datetime, "%Y-%m-%d %H:%M:%S")
    }

    # Convert the raster to a dataframe
    if (i != "NDVI") {
      map <- disagg(subdat[[x]], fact = 4, method = "bilinear")
    } else {
      map <- aggregate(subdat[[x]], fact = 4, method = "mean")
    }
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
          # colors = c("white", rev(brewer.pal(11, name = "RdYlBu")))
          colors = rev(cols)
        , limits = c(minmax[1], minmax[2])
        , oob    = squish
        , guide  = guide_colorbar(
          , title          = labs
          , show.limits    = T
          , title.position = "top"
          , title.hjust    = 0.5
          , ticks          = T
          , barheight      = unit(0.3, "cm")
          , barwidth       = unit(10.0, "cm")
          , direction      = "horizontal"
        )
      ) +
      labs(
          x        = NULL
        , y        = NULL
        , fill     = NULL
        , title    = labs
        , subtitle = datetime
      ) +
      theme_void() +
      theme(
          legend.position   = c(0.7, 1.05)
        , panel.background  = element_rect(color = "gray20")
        , legend.text       = element_text(color = "white")
        , legend.title      = element_blank()
        , legend.key        = element_blank()
        , legend.background = element_rect(fill = "black")
        , panel.grid.minor  = element_blank()
        , panel.grid.major  = element_blank()
        , plot.background   = element_rect(fill = "black")
        , plot.title        = element_text(color = "white")
        , plot.subtitle     = element_text(color = "white")
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

    # Store the plot to file
    suppressWarnings(
      ggsave(
          plot     = p
        , filename = plot_names[x]
        , width    = 1000
        , height   = 900
        , units    = "px"
        , bg       = "transparent"
        , dpi      = 150
      )
    )

  }))

  # Randomly visualize a few of them
  # plot(image_read(sample(plot_names, 1)))

  # Animate the plots
  cat("Preparing animation for", tolower(i), "\n")
  ani.options(interval = 1 / 10)
  im.convert(plot_names, filename)

#   # Prepare an animation
#   plot_names %>%
#     image_read() %>%
#     image_morph(10) %>%
#     image_animate(fps = 20, optimize = T) %>%
#     image_write(path = "Test.gif")


}

################################################################################
#### Flood
################################################################################
# Load the data
dat <- dir(
    path    = paste0("03_Data/02_CleanData/00_Floodmaps/02_Resampled")
  , pattern = ".tif$", full.names = T
) %>% rast()

# Get an extent
ext <- ext(dat)

# We also want to add other waterbodies in the background
water1 <- rast("/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/02_CleanData/01_LandCover_WaterCoverStatic.tif")
water2 <- rast("/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/02_CleanData/03_LandscapeFeatures_Rivers.tif")
water1[ext] <- 0
water <- max(water1, water2)

# Retrieve dates from layernames
dates <- names(dat) %>%
  substr(start = 1, stop = 10) %>%
  ymd()

# Let's determine the flood extent in each image
flood_summary <- dat %>%
  expanse(unit = "km", byValue = T) %>%
  as.data.frame() %>%
  pivot_wider(
    , id_cols     = layer
    , names_from  = value
    , values_from = area
    , values_fill = 0
  ) %>%
  rename(Flood = "0", Dryland = "1", Cloud = "2") %>%
  mutate(Total = Flood + Dryland + Cloud) %>%
  mutate(Cloud = Cloud / Total)

# Subset to a selected date range
subdat <- dat[[which(flood_summary$Cloud < 0.10)]]
subdat <- subdat[[1:300]]
subdat <- subst(subdat, 2, 0)
subnam <- names(subdat)

# Put layers together
subdat <- extend(subdat, water)
subdat <- mask(water, subdat, maskvalue = 1, updatevalue = 1)
names(subdat) <- subnam

# Prepare the frames
filename <- paste0("05_Presentation/", "Flood", ".gif")
if (!file.exists(filename)) {
  cat("Preparing plots for floodmaps\n")
  plot_names <- paste0(tempfile(), "_", 1:nlyr(subdat), ".png")
  invisible(pbmclapply(1:nlyr(subdat), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {

    # Get the date and prepare a title from it
    datetime <- names(subdat[[x]])
    date <- substr(datetime, start = 1, stop = 10)

    # Convert the raster to a dataframe
    map <- aggregate(subdat[[x]], fact = 4, fun = "max")
    map <- as.data.frame(map, xy = T)
    names(map)[3] <- "value"

    # Make the plot
    p <- ggplot(map, aes(x = x, y = y, fill = as.factor(value))) +
      geom_raster() +
      geom_sf(data = africa, inherit.aes = F, fill = NA, col = "gray70") +
      geom_sf_text(data = africa, inherit.aes = F, aes(label = COUNTRY), nudge_y = 0.1, col = "gray70") +
      coord_sf(
          crs    = 4326
        , xlim   = c(min(r$x), max(r$x))
        , ylim   = c(min(r$y), max(r$y))
        , expand = F
      ) +
      scale_fill_manual(
          values = c("gray20", "cornflowerblue")
        , guide  = guide_colorbar(
          , title          = "Water"
          , show.limits    = T
          , title.position = "top"
          , title.hjust    = 0.5
          , ticks          = T
          , barheight      = unit(0.3, "cm")
          , barwidth       = unit(10.0, "cm")
        )
      ) +
      labs(
          x        = NULL
        , y        = NULL
        , fill     = NULL
        , title    = "Water"
        , subtitle = datetime
      ) +
      theme_void() +
      theme(
          legend.position   = "bottom"
        , legend.box        = "vertical"
        , panel.background  = element_rect(color = "gray20")
        , legend.text       = element_text(color = "white")
        , legend.title      = element_blank()
        , legend.key        = element_blank()
        , legend.background = element_rect(fill = "black")
        , panel.grid.minor  = element_blank()
        , panel.grid.major  = element_blank()
        , plot.background   = element_rect(fill = "black")
        , plot.title        = element_text(color = "white")
        , plot.subtitle     = element_text(color = "white")
      ) +
      annotation_scale(
          location   = "bl"
        , width_hint = 0.2
        , line_width = 0.5
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

    # Store the plot to file
    suppressWarnings(
      ggsave(
          plot     = p
        , filename = plot_names[x]
        , width    = 1000
        , height   = 900
        , units    = "px"
        , bg       = "transparent"
        , dpi      = 150
      )
    )

  }))

  # Randomly visualize a few of them
  # plot(image_read(sample(plot_names, 1)))

  # Animate the plots
  cat("Preparing animation for the flood\n")
  ani.options(interval = 1 / 10, ani.width = 1000, ani.height = 900)
  im.convert(plot_names, filename)
}
