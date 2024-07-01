################################################################################
#### Connectivity
################################################################################
# Description: Visualization of Connectivity

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(tidyverse)    # To wrangle data
library(terra)        # To handle spatial data
library(sf)           # To visualize spatial data
library(RColorBrewer) # For color scales
library(ggspatial)    # For scale bar and north arrow
library(scales)       # To squish colorscales
library(igraph)       # For network plots
library(ggnetwork)    # For network plots
library(colorspace)   # To darken and lighten colors

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load shapefiles that we want to plot
water <- read_sf("03_Data/02_CleanData/MajorWaters.gpkg")
areas <- read_sf("03_Data/02_CleanData/Sources.gpkg") %>% mutate(ID = 1:n())
roads <- read_sf("03_Data/02_CleanData/Roads.gpkg")
afric <- read_sf("03_Data/02_CleanData/Africa.gpkg")
vills <- read_sf("03_Data/02_CleanData/Villages.gpkg")
vills <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y)

# Prepare centroids
lay <- st_point_on_surface(areas)
lay <- st_coordinates(lay)

# Create country labels
labels_countries <- data.frame(
    x     = c(25.5, 26, 25.7, 21.5, 23.5)
  , y     = c(-19.3, -18.2, -17.6, -17.6, -17.8)
  , Label = c("BOTSWANA", "ZIMBABWE", "ZAMBIA", "ANGOLA", "NAMIBIA")
)

# areas$Name <- gsub(areas$Name, pattern = "\\s", replacement = " \\\n")

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.1, 25.6)
  , y     = c(-19.1, -18.2, -17.5, -20.7)
  , Label = c("Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans")
)

# Generate area labels
labels_areas <- st_coordinates(st_point_on_surface(areas))
labels_areas <- cbind(labels_areas, st_drop_geometry(areas))

# Reload connectivity results
dat <- "03_Data/03_Results/Connectivity.rds" %>%
  read_rds() %>%
  mutate(Data = map(Filename, function(x) {
    dat             <- read_rds(x)
    dat$Heatmap     <- unwrap(dat$Heatmap)
    dat$Betweenness <- unwrap(dat$Betweenness)
    return(dat)
  })) %>%
  mutate(ModelCode = substr(ModelCode, start = 1, stop = 3)) %>%
  mutate(ModelCode = factor(ModelCode, levels = c("SSS", "DMD"))) %>%
  mutate(Formula = factor(Formula, levels = c("Simple", "Full"), labels = c("Simple Formula", "Complex Formula")))

# Extract interpatch connectivity
ipc <- dat %>% mutate(Data = map(Data, function(x) {
  x$Interpatch
  })) %>%
  unnest(Data) %>%
  rename(
      SourceArea        = source
    , CurrentArea       = area
    , DispersalDuration = duration
    , DispersalSuccess  = success
  )

# ggplot(ipc, aes(x = SourceArea, y = CurrentArea, fill = DispersalDuration)) +
#   geom_tile() +
#   facet_grid(Formula ~ ModelCode) +
#   coord_equal()

# Sum heatmaps under the same configurations
heat <- dat %>%
  mutate(Data = map(Data, function(x) {
    x$Heatmap
  })) %>%
  nest(Data = -c(ModelCode, Formula)) %>%
  mutate(Data = map(Data, function(x) {
    result <- do.call(c, x$Data)
    result <- sum(result)
    return(result)
  }))

# # Compute map correlations
# heat %>%
#   subset(ModelCode == "SSS") %>%
#   pull(Data) %>%
#   do.call(c, .) %>%
#   layerCor(fun = "cor")
# heat %>%
#   subset(ModelCode == "DMD") %>%
#   pull(Data) %>%
#   do.call(c, .) %>%
#   layerCor(fun = "cor")

# Some smoothing
heat <- heat %>%
  mutate(Data = map(Data, function(x) {
    result <- disagg(x, fact = 4, method = "bilinear")
    result <- focal(result, w = matrix(1, 5, 5), fun = mean)
    result <- aggregate(result, fact = 2, fun = "mean")
    return(result)
  }))

# Quick plot
plot(do.call(c, heat$Data), main = paste(heat$Formula, heat$ModelCode))

# Create a single dataframe
heat <- heat %>%
  mutate(Data = map(Data, function(x) {
    result <- as.data.frame(x, xy = T)
    names(result)[3] <- "Value"
    return(result)
  })) %>% unnest(Data)
#
# library(raster)
# library(gstat)
# test <- heat %>%
#   nest(Data = -c(Formula, ModelCode)) %>%
#   mutate(Variogram = map(Data, function(x) {
#     pts   <- vect(x, geom =c("x", "y"))
#     pts   <- as(pts, "Spatial")
#     vario <- variogram(object = Value ~ 1, data = pts)
#     return(vario)
#   }))
#
# # Fit models
# test <- test %>%
#   mutate(Model = map(Variogram, function(x) {
#     var <- fit.variogram(x, vgm("Sph"))
#     result <- data.frame(
#         Param = c("Sill", "Range", "Nugget")
#       , Value = c(var$psill[2], var$range[2], var$psill[1])
#     )
#     return(result)
#   }))
#
# test %>%
#   dplyr::select(Formula, ModelCode, Model) %>%
#   unnest(Model)
#
# test %>%
#   dplyr::select(Formula, ModelCode, Variogram) %>%
#   unnest(Variogram) %>%
#   ggplot(aes(x = dist, y = gamma, color = ModelCode, linetype = Formula)) +
#     geom_line() +
#     theme_awesome()

################################################################################
#### Functions to Plot
################################################################################
# Function to plot the different metrics (except for the inter-patch
# connectivity)
spectral <- colorRampPalette(rev(brewer.pal(11, name = "Spectral")))
plotMet <- function(
    data
  , formula     = ~ ModelCode
  , area        = NULL
  , barwidth    = unit(16, "cm")
  , colorscheme = c("dark", "light")
  , palette     = magma(100)
  , col_area    = c("black", "white")
  , col_country = c("black", "white")
  , col_water   = NA
  , fill_water  = NA
  , name        = NULL
  , trans       = "identity"
  , limits      = NULL
  , ext         = NULL
  , lomehi      = F
  , cutlines    = F
  , diffmap     = F
  ) {

    # Testing
    # data        <- heat
    # formula     <- ~ ModelCode
    # barwidth    <- 16
    # labels      <- F
    # colorscheme <- "dark"
    # col_area    <- c("black")
    # col_country <- c("black")
    # palette     <- spectral(100)
    # name        <- "Test"
    # area        <- 1
    # trans       <- "identity"
    # limits      <- range(data$Value)
    # col_water   <- "gray50"
    # fill_water  <- NA
    # lomehi      <- T
    # cutlines    <- F
    # ext <- ext(dat$Data[[1]]$Heatmap)
    # diffmap <- F

  # Some checks
  if (diffmap) {
      legend_labels <- c("< 0", "0", "> 0")
    } else{
      legend_labels <- c("Low", "Medium", "High")
  }
  breaks          <- c(
      min(data$Value, na.rm = T)
    , (min(data$Value, na.rm = T) + max(data$Value, na.rm = T)) / 2
    , max(data$Value, na.rm = T)
  )
  if (!is.null(limits)) {
    breaks <- c(min(limits), mean(limits), max(limits))
  }
  if (is.null(ext)) {
    ext <- extent(r)
  }

  # The plot
  ggplot() +
    geom_raster(
        data    = data
      , mapping = aes(x = x, y = y, fill = Value)
    ) +
    geom_sf(
        data  = water
      , fill  = fill_water
      , col   = col_water
      , alpha = 0.25
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
        data        = areas
      , color       = col_area
      , fill        = col_area
      , lty         = 1
      , linewidth   = 0.2
      , show.legend = F
      , alpha       = 0.15
    ) +
    geom_sf(
        data      = afric
      , linewidth = 0.4
      , col       = col_country
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
      , col      = col_country
      , size     = 2
      , fontface = 2
    ) +
    geom_text(
        data     = labels_areas
      , mapping  = aes(x = X, y = Y, label = Name)
      , col      = col_area
      , fontface = 3
      , size     = 2
    ) +
    geom_text(
        data     = labels_areas
      , mapping  = aes(x = X, y = Y, label = Name)
      , col      = col_area
      , fontface = 3
      , size     = 2
    ) +
    geom_text(
        data     = subset(vills, place == "City")
      , mapping  = aes(x = x, y = y, label = name)
      , col      = ifelse(colorscheme == "dark", "gray90", "gray30")
      , fontface = 3
      , size     = 2
      , nudge_y  = c(0.1, -0.1, 0.1)
    ) +
    scale_size_manual(values = c(1.0, 0.25)) +
    scale_color_manual(values = c(col_area, "red")) +
    {
      if (lomehi) {
        scale_fill_gradientn(
            colors  = palette
          , breaks  = breaks
          , labels  = legend_labels
          , trans   = trans
          , limits  = limits
          , oob     = squish
          , guide   = guide_colorbar(
            , title          = name
            , show.limits    = T
            , title.position = "bottom"
            , title.hjust    = 0.5
            , ticks          = F
            , barheight      = unit(0.2, "cm")
            , barwidth       = barwidth
          )
        )
      } else {
        scale_fill_gradientn(
            colors  = palette
          # , breaks  = c(0, 125, 250)
          # , labels  = c("Low", "Medium", "High")
          , labels  = function(x){format(x, big.mark = "'")}
          , trans   = trans
          , limits  = limits
          , oob     = squish
          , guide   = guide_colorbar(
            , title          = name
            , show.limits    = T
            , title.position = "bottom"
            , title.hjust    = 0.5
            , ticks          = F
            , barheight      = unit(0.2, "cm")
            , barwidth       = barwidth
          )
        )
      }
    } +
    scale_x_continuous(breaks = seq(21.5, 26.5, by = 1)) +
    scale_y_continuous(breaks = seq(-20.5, -17.5, by = 1)) +
    coord_sf(
        crs    = 4326
      , xlim   = ext[1:2]
      , ylim   = ext[3:4]
      , expand = F
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
    ) +
    facet_grid(formula)
}

# Generate Heatmaps
ext <- ext(dat$Data[[1]]$Heatmap)
p <- plotMet(subset(heat, Formula == "Complex Formula")
  , formula     = ~ ModelCode
  , barwidth    = unit(8, "cm")
  , colorscheme = "dark"
  , col_area    = "black"
  , col_country = "black"
  , col_water   = NA
  , fill_water  = NA
  , palette     = spectral(100)
  , name        = "# Traversing Trajectories"
  , trans       = "identity"
  , ext         = ext - c(0, 0, 0, 1)
  , lomehi      = T
)

# Recolor text
p <- p + theme(
    axis.text         = element_text(color = "white")
  , axis.ticks        = element_line(color = "white")
  , legend.text       = element_text(color = "white")
  , legend.title      = element_text(color = "white")
  , panel.background  = element_blank()
  , plot.background   = element_blank()
  , legend.background = element_blank()
)

# Store plot to file
ggsave("05_Presentation/HeatmapComparison.png"
  , plot   = p
  , bg     = "transparent"
  , height = 4
  , width  = 6
  , scale  = 1.6
  , device = png
)

# ################################################################################
# #### Interpatch Connectivity
# ################################################################################
# # Function to plot the inter-patch connectivity maps
# plotNet <- function(
#     sources     = NULL
#   , targets     = NULL
#   , pairs       = NULL
#   , reverse     = F
#   , area        = NULL
#   , col_area    = "grey"
#   , curvature   = 0.2
#   , col_water   = NA
#   , fill_water  = NA
#   , alpha_water = 1
#   , ext         = NULL
#   , max_size    = 1
#   , trans       = "identity"
#   , barwidth    = unit(8, "cm")
#   , days        = T
#   ) {
#
#   # Testing
#   formula <- ~ ModelCode
#   sources <- NULL
#   targets <- NULL
#   pairs   <- rbind(c(1, 2), c(2, 3), c(1, 3))
#   barwidth    = unit(10, "cm")
#   reverse     = T
#   fill_water  = "cornflowerblue"
#   alpha_water = 1
#   ext         = ext
#   max_size    = 1.5
#   col_water   = NA
#   col_area    = "black"
#   area        <- 1
#   curvature   = 0.2
#   trans    <- "identity"
#
#   # Subset to desired data
#   sub <- select(ipc, SourceArea, CurrentArea, DispersalSuccess, DispersalDuration, ModelCode)
#   if (!is.null(sources) & !is.null(targets)) {
#       subi <- subset(sub, SourceArea %in% sources & CurrentArea %in% targets)
#       if (reverse) {
#         subi2 <- subset(sub, SourceArea %in% targets & CurrentArea %in% sources)
#         subi <- rbind(subi, subi2)
#       }
#       sub <- subi
#       sub <- subset(sub, SourceArea != CurrentArea)
#     } else if (!is.null(pairs)) {
#       all <- list()
#       for (i in 1:nrow(pairs)) {
#         subi <- subset(sub, SourceArea == pairs[i, 1] & CurrentArea == pairs[i, 2])
#         if (reverse) {
#           subi2 <- subset(sub, SourceArea == pairs[i, 2] & CurrentArea == pairs[i, 1])
#           subi <- rbind(subi, subi2)
#         }
#         all[[i]] <- subi
#       }
#       sub <- do.call(rbind, all)
#       sub <- subset(sub, SourceArea != CurrentArea)
#   }
#
#   # Generate a network
#   net <- graph_from_data_frame(
#       d        = sub
#     , vertices = unique(areas$ID)
#     , directed = T
#   )
#
#   # Prepare networks for ggplotting with ggplot
#   net_p <- suppressWarnings(ggnetwork(net, layout = lay, arrow.gap = 0.1, scale = F))
#   net_p <- subset(net_p, !is.na(ModelCode))
#
#   # Plot
#   ggplot() +
#     geom_sf(
#         data  = water
#       , fill  = fill_water
#       , col   = col_water
#       , alpha = alpha_water
#     ) +
#     # geom_sf(
#     #     data  = subset(areas, Type == "Main")
#     #   , col   = "orange"
#     #   , fill  = "orange"
#     #   , alpha = 0.2
#     # ) +
#     # geom_sf(
#     #     data  = subset(areas, Type == "Buffer")
#     #   , col   = "purple"
#     #   , fill  = "purple"
#     #   , alpha = 0.2
#     # ) +
#     geom_sf(
#       #         data        = subset(areas, Type == "Main")
#         data        = areas
#       , color       = col_area
#       , fill        = col_area
#       , lty         = 1
#       , linewidth   = 0.2
#       , show.legend = F
#       , alpha       = 0.15
#     ) +
#     {
#       if (!is.null(area)) {
#         geom_sf(
#             data      = subset(areas, ID %in% area)
#           , col       = "red"
#           , fill      = NA
#           , linewidth = 0.5
#         )
#       }
#     } +
#     geom_text(
#         data     = labels_waters
#       , mapping  = aes(x = x, y = y, label = Label)
#       , col      = darken(fill_water, 1.5)
#       , fontface = 3
#       , size     = 1.3
#     ) +
#     geom_edges(
#         data      = net_p
#       , mapping   = aes(
#           x    = x
#         , y    = y
#         , xend = xend
#         , yend = yend
#         , size = DispersalSuccess
#         , col  = DispersalDuration
#       )
#       , curvature = curvature
#       , arrow     = arrow(length = unit(6, "pt"), type = "closed", angle = 10)
#     ) +
#     geom_text(
#         data     = labels_areas
#       , mapping  = aes(x = X, y = Y, label = Name)
#       , col      = "black"
#       , fontface = 3
#       , size     = 2
#     ) +
#     scale_fill_manual(
#         values = c("transparent", "cornflowerblue")
#       , name   = ""
#     ) +
#     facet_wrap(formula) +
#     coord_sf(
#         crs    = 4326
#       , xlim   = ext[1:2]
#       , ylim   = ext[3:4]
#       , expand = F
#     ) +
#     scale_size_area(
#         name     = "Dispersal Frequency\n(number of trajectories)"
#       , max_size = max_size
#       , breaks   = c(200, 400, 600)
#       , labels   = c(200, 400, 600)
#       , guide    = guide_legend(
#           order          = 2
#         , title.hjust    = 0.5
#         , title.position = "bottom"
#         # , label.position = "bottom"
#       )
#     ) +
#     scale_color_gradientn(
#         colors = rev(brewer.pal(11, "RdYlGn"))
#       , trans  = trans
#       , guide  = guide_colorbar(
#         , title          = "Duration"
#         , show.limits    = T
#         , title.position = "bottom"
#         , title.hjust    = 0.5
#         , ticks          = F
#         , barheight      = unit(0.2, "cm")
#         , barwidth       = barwidth
#       )
#     ) +
#     scale_x_continuous(breaks = seq(21.5, 25.5, by = 1)) +
#     scale_y_continuous(breaks = seq(-20.5, -18.5, by = 1)) +
#     theme(
#         panel.background = element_blank()
#       , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
#       , legend.position  = "bottom"
#       , legend.box       = "horizontal"
#       , legend.key       = element_blank()
#       , strip.background = element_rect(fill = "gray95", color = "transparent")
#     ) +
#     annotation_scale(
#         location   = "bl"
#       , width_hint = 0.2
#       , line_width = 0.5
#       , text_cex   = 0.5
#       , height     = unit(0.1, "cm")
#     ) +
#     annotation_north_arrow(
#         location = "br"
#       , height   = unit(0.7, "cm"),
#       , width    = unit(0.6, "cm"),
#       , style    = north_arrow_fancy_orienteering(
#             fill      = "black"
#           , line_col  = NA
#           , text_col  = "black"
#           , text_size = 4
#         )
#     ) +
#     xlab("") +
#     ylab("")
# }
#
#
# # Visualize heatmaps
# p1 <- ggplot(as.data.frame(heat[[2]], xy = T), aes(x = x, y = y, fill = DMD_F)) +
#   geom_raster() +
#   scale_fill_distiller(palette = "Spectral")
#
# p2 <- as.data.frame(heat[[2]], xy = T) %>%
#   group_by(x) %>%
#   summarize(value = mean(DMD_F)) %>%
#   ggplot(aes(x = x, y = value)) +
#     geom_line() +
#     scale_y_reverse()
# p3 <- as.data.frame(heat[[2]], xy = T) %>%
#   group_by(y) %>%
#   summarize(value = mean(DMD_F)) %>%
#   ggplot(aes(x = value, y = y)) +
#     geom_path()
#
# library(cowplot)
# plot_grid(
#   p1 + theme(legend.position = "none"), p3, p2
#   , nrow = 2
#   , align = "hv"
#   , axis = "t"
#   , rel_widths = c(1, 0.3)
#   , rel_heights = c(1, 0.3)
# )
