################################################################################
#### Time Kernel of Dispersal
################################################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(raster)         # To handle raster data
library(rgdal)          # To load shapefiles
library(tidyverse)      # To wrangle data
library(davidoff)       # Custom functions
library(RColorBrewer)   # For colors
library(viridis)        # For colors
library(sf)             # To plot spatial objects with ggplot
library(ggspatial)      # For north arrow and scale bar
library(cowplot)        # To grab legends
library(ggpubr)         # To arrange ggplots

################################################################################
#### Load Required Data
################################################################################
# Load simulations
sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Load betweenness maps
betweenness <- stack("03_Data/03_Results/99_Betweenness.grd")

# Apply focal filter to buffer/smooth maps
betweenness <- lapply(1:nlayers(betweenness), function(x){
  focal(betweenness[[x]], w = matrix(1, 3, 3), fun = mean)
}) %>% stack()

# Load shapefiles for background
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
africa  <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
prot    <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Load reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Convert to sf for plotting with ggplot
kaza    <- st_as_sf(kaza)
africa  <- st_as_sf(africa)

# Convert maps to dataframes
betweenness <- lapply(1:nlayers(betweenness), function(x){
  df <- as.data.frame(betweenness[[x]], xy = T)
  names(df) <- c("x", "y", "layer")
  return(df)
})

# Set names
names(betweenness) <- c("125", "500", "2000")

# Also convert the reference raster
r <- as.data.frame(r, xy = T)

################################################################################
#### Plot
################################################################################
# Prepare color palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

# Main Plot
p1 <- ggplot() +
  geom_raster(
      data    = betweenness[[3]]
    , mapping = aes(x = x, y = y, fill = sqrt(layer))
  ) +
  scale_fill_gradientn(
      colours = magma(100)
    , guide   = guide_colorbar(
      , title          = expression("Betweenness Score" ^ "0.5")
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  geom_sf(
      data        = kaza
    , col         = "white"
    , fill        = NA
    , lty         = 1
    , lwd         = 0.5
    , show.legend = F
  ) +
  geom_sf(
      data        = africa
    , col         = "white"
    , fill        = NA
    , lty         = 2
    , lwd         = 0.2
    , show.legend = F
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    , title    = "Betweenness"
    , subtitle = "After 2000 Steps"
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
    , line_width = 1
    , height     = unit(0.15, "cm")
    , bar_cols   = c("white", "white")
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

# Plot for separate legend
p2 <- ggplot() +
  geom_raster(
      data        = betweenness[[3]]
    , mapping     = aes(x = x, y = y, fill = layer)
    , show.legend = F
  ) +
  geom_sf(
      data        = kaza
    , mapping     = aes(color = "KAZA-TFCA Borders")
    , fill        = NA
    , lty         = 1
    , lwd         = 0.5
    , show.legend = "line"
  ) +
  geom_sf(
      data        = africa
    , mapping     = aes(color = "Country Borders")
    , fill        = NA
    , lty         = 2
    , lwd         = 0.2
    , show.legend = "line"
  ) +
  scale_color_manual(
      values = c("Country Borders" = "white", "KAZA-TFCA Borders" = "white")
    , guide  = guide_legend(
        override.aes = list(
          linetype = c(2, 1)
        , shape    = c(NA, NA)
        , lwd      = c(0.2, 0.5)
      )
    )
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  theme(
      legend.position       = c(0.20, 0.90)
    , legend.title          = element_blank()
    , legend.box            = "vertical"
    , legend.background     = element_blank()
    , legend.box.background = element_rect(fill  = "black", color = "white")
    , legend.margin         = margin(0, 8, 2, 6)
    , legend.text           = element_text(color = "white")
    , legend.key            = element_blank()
    , legend.key.size       = unit(0.8, "lines")
    , legend.key.width      = unit(1.2, "lines")
    , panel.background      = element_blank()
  )

# Extract legend
legend <- get_legend(p2)

# Put into main plot
p3 <- p1 + annotation_custom(
      grob = legend
    , xmin = 18.75
    , xmax = 21
    , ymin = -13
    , ymax = -13.5
  )

# Show it
p3

# Store it to file
ggsave("04_Manuscript/99_Betweenness.png", plot = p3)

################################################################################
#### Plots over Time
################################################################################
# Write a function to plot a heatmap
plotBetweenness <- function(x, subtitle = NULL, legend = T, barwidth = 10){

  # Prepare dataframe
  x <- as.data.frame(x, xy = T)
  names(x) <- c("x", "y", "layer")

  # Main Plot
  p1 <- ggplot() +
    geom_raster(
        data        = x
      , mapping     = aes(x = x, y = y, fill = sqrt(layer))
    ) +
    scale_fill_gradientn(
        colours = magma(100)
      , guide   = guide_colorbar(
        , title          = expression("Betweenness Score" ^ "0.5")
        , show.limits    = T
        , title.position = "top"
        , title.hjust    = 0.5
        , ticks          = T
        , barheight      = unit(0.6, "cm")
        , barwidth       = unit(10.0, "cm")
      )
    ) +
    geom_sf(
        data        = kaza
      , col         = "white"
      , fill        = NA
      , lty         = 1
      , lwd         = 0.5
      , show.legend = F
    ) +
    geom_sf(
        data        = africa
      , col         = "white"
      , fill        = NA
      , lty         = 2
      , lwd         = 0.2
      , show.legend = F
    ) +
    coord_sf(
        crs    = 4326
      , xlim   = c(min(r$x), max(r$x))
      , ylim   = c(min(r$y), max(r$y))
      , expand = F
    ) +
    labs(
        x        = NULL
      , y        = NULL
      , fill     = NULL
      , title    = "Betweenness"
      , subtitle = subtitle
    ) +
    theme(
        legend.position      = "top"
      , legend.justification = "right"
      , legend.box           = "vertical"
      , legend.box.margin    = margin(-10, 0, -10, 0)
      , panel.background     = element_blank()
      , panel.border         = element_rect(colour = "black", fill = NA, size = 1)
      , plot.subtitle        = element_text(margin = margin(t = 0, b = -30))
    ) +
    annotation_scale(
        location   = "bl"
      , width_hint = 0.2
      , line_width = 1
      , height     = unit(0.15, "cm")
      , bar_cols   = c("white", "white")
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

  if (legend){

    # Plot for separate legend
    p2 <- ggplot() +
      geom_raster(
          data        = x
        , mapping     = aes(x = x, y = y, fill = layer)
        , show.legend = F
      ) +
      geom_sf(
          data        = kaza
        , mapping     = aes(color = "KAZA-TFCA Borders")
        , fill        = NA
        , lty         = 1
        , lwd         = 0.5
        , show.legend = "line"
      ) +
      geom_sf(
          data        = africa
        , mapping     = aes(color = "Country Borders")
        , fill        = NA
        , lty         = 2
        , lwd         = 0.2
        , show.legend = "line"
      ) +
      scale_color_manual(
          values = c("white", "white")
        , guide  = guide_legend(
            override.aes = list(
              linetype = c(2, 1)
            , shape    = c(NA, NA)
            , lwd      = c(0.2, 0.5)
          )
        )
      ) +
      coord_sf(
          crs    = 4326
        , xlim   = c(min(x$x), max(x$x))
        , ylim   = c(min(x$y), max(x$y))
        , expand = F
      ) +
      xlim(min(x$x), max(x$x)) +
      ylim(min(x$y), max(x$y)) +
      theme(
          legend.position       = c(0.20, 0.90)
        , legend.title          = element_blank()
        , legend.box            = "vertical"
        , legend.background     = element_blank()
        , legend.box.background = element_rect(fill  = "black", color = "white")
        , legend.margin         = margin(0, 8, 2, 6)
        , legend.text           = element_text(color = "white")
        , legend.key            = element_blank()
        , legend.key.size       = unit(0.8, "lines")
        , legend.key.width      = unit(1.2, "lines")
        , panel.background      = element_blank()
      )

    # Extract legend
    legend <- get_legend(p2)

    # Put into main plot
    p3 <- p1 + annotation_custom(
          grob = legend
        , xmin = 18.75
        , xmax = 21
        , ymin = -13
        , ymax = -13.5
      )
  } else {
    p3 <- p1
  }

  # Return it
  return(p3)

}

# Let's apply the function to get all desired plots
maps <- lapply(1:length(betweenness), function(x){
  subtitle <- paste0("After ", names(betweenness)[[x]], " Steps")
  map <- plotBetweenness(
      x        = betweenness[[x]]
    , subtitle = subtitle
    , legend   = F
    , barwidth = 7
  )
  return(map)
})

# Arrange plots nicely
p <- ggarrange(maps[[1]], maps[[2]], maps[[3]], nrow = 1)

# Store the arranged plot
ggsave("04_Manuscript/99_BetweennessIndividual.png"
  , plot   = p
  , scale  = 2
  , height = 3
  , width  = 9
)
