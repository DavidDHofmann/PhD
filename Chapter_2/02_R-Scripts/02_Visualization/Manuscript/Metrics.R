################################################################################
#### Metrics
################################################################################
# Description: Visualization of all spatial metrics

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(terra)          # To handle spatial data
library(raster)         # To handle spatial data
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(ggspatial)      # For scale bar and north arrow
library(igraph)         # For network analysis
library(ggnetwork)      # To plot network using ggplot
library(rgdal)          # To handle spatial data
library(sf)             # To handle spatial data
library(scales)         # To squish oob values
library(latex2exp)      # For easy latex code
library(RColorBrewer)   # For custom colors
library(ggpubr)         # To arrange multiple plots
library(viridis)        # For color palettes
library(smoothr)        # To smooth geometries

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Load Data
################################################################################
# Load shapefiles that we want to plot
areas <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
intre <- read_sf("03_Data/02_CleanData/AreasOfInterest.shp")
roads <- read_sf("03_Data/02_CleanData/Roads.shp")
afric <- read_sf("03_Data/02_CleanData/Africa.shp")
vills <- read_sf("03_Data/02_CleanData/Villages.shp")
cutli <- read_sf("03_Data/02_CleanData/Cutlines.shp")[c(4, 8), ]
cutli <- st_crop(cutli, st_bbox(c(xmin = 21.5, xmax = 23.6, ymin = -20.1, ymax = -18.2)))
vills <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y)

# Load results on interpatch connectivity
ipc <- read_rds("03_Data/03_Results/InterpatchConnectivityBootstrapped.rds")
ipc <- subset(ipc, SourceArea != "Overall" & CurrentArea != "Overall" & SourceArea != CurrentArea)

# Prepare centroids
lay <- st_point_on_surface(areas)
lay <- st_coordinates(lay)

# Also prepare water layers
water <- rast("03_Data/02_CleanData/WaterCover.tif")[[c("min", "max")]]
water_min <- st_as_sf(as.polygons(subst(water[[1]], 0, NA)))
water_max <- st_as_sf(as.polygons(subst(water[[2]], 0, NA)))
names(water_min)[1] <- "Water"
names(water_max)[1] <- "Water"
water <- rbind(water_min, water_max)
water$FloodLevel <- factor(c("Min", "Max"), levels = c("Min", "Max"))

# Prepare a custom color ramp (I don't like ggplots version of spectral)
spectral <- colorRampPalette(rev(brewer.pal(11, name = "Spectral")))

# Load the reference raster
r <- raster("03_Data/02_CleanData/ReferenceRaster.tif")

# Create country labels
labels_countries <- data.frame(
    x     = c(25.5, 26, 25.7, 21.5, 23.5)
  , y     = c(-19.3, -18.2, -17.6, -17.6, -17.8)
  , Label = c("BOTSWANA", "ZIMBABWE", "ZAMBIA", "ANGOLA", "NAMIBIA")
)

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.1, 25.6)
  , y     = c(-19.1, -18.2, -17.5, -20.7)
  , Label = c("Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans")
)

# Generate area labels
labels_areas <- st_coordinates(st_point_on_surface(areas))
labels_areas <- cbind(labels_areas, st_drop_geometry(areas))

################################################################################
#### Preparing Heatmaps
################################################################################
# Load the heatmaps and keep only desired columns
maps1 <- "03_Data/03_Results/HeatmapsBetweennessGlobal.rds" %>%
  read_rds() %>%
  mutate(SourceArea = NA, Level = "Global") %>%
  select(Steps, SourceArea, FloodLevel, Level, Data = Heatmap) %>%
  arrange(Steps, SourceArea, FloodLevel, Level)
maps2 <- "03_Data/03_Results/HeatmapsBetweennessLocal.rds" %>%
  read_rds() %>%
  mutate(Level = "Local") %>%
  select(Steps, SourceArea, FloodLevel, Level, Data = Heatmap) %>%
  arrange(Steps, SourceArea, FloodLevel, Level)
heat <- rbind(maps1, maps2)

# I also want to plot a difference map
diff <- subset(heat, Level == "Global" & Steps == 2000)
diff <- diff$Data[[1]] - diff$Data[[2]]
diff <- tibble(
    Steps      = 2000
  , SourceArea = NA
  , FloodLevel = "Difference"
  , Level      = "Difference"
  , Data       = list(diff)
)
heat <- rbind(heat, diff)

# Convert the maps to dataframes
heat <- heat %>%
  mutate(Data = map(Data, function(x) {
    result <- as.data.frame(x, xy = T)
    names(result) <- c("x", "y", "Value")
    return(result)
  })) %>% unnest(Data)

# Make sure the levels are correctly ordered
heat$FloodLevel <- factor(heat$FloodLevel, levels = c("Min", "Max", "Difference"))
heat$Metric     <- "Heatmap"

################################################################################
#### Preparing Betweenness Maps
################################################################################
# Load the heatmaps and keep only desired columns
maps1 <- "03_Data/03_Results/HeatmapsBetweennessGlobal.rds" %>%
  read_rds() %>%
  mutate(SourceArea = NA, Level = "Global") %>%
  select(Steps, SourceArea, FloodLevel, Level, Data = Betweenness) %>%
  arrange(Steps, SourceArea, FloodLevel, Level)
maps2 <- "03_Data/03_Results/HeatmapsBetweennessLocal.rds" %>%
  read_rds() %>%
  mutate(Level = "Local") %>%
  select(Steps, SourceArea, FloodLevel, Level, Data = Betweenness) %>%
  arrange(Steps, SourceArea, FloodLevel, Level)
betw <- rbind(maps1, maps2)

# I also want to plot a difference map
diff <- subset(betw, Level == "Global" & Steps == 2000)
diff <- diff$Data[[1]] - diff$Data[[2]]
diff <- tibble(
    Steps      = 2000
  , SourceArea = NA
  , FloodLevel = "Difference"
  , Level      = "Difference"
  , Data       = list(diff)
)
betw <- rbind(betw, diff)

# Convert maps to dataframes
betw <- betw %>%
  mutate(Data = map(Data, function(x) {
    result <- as.data.frame(x, xy = T)
    names(result) <- c("x", "y", "Value")
    return(result)
  })) %>% unnest(Data)

# Make sure the levels are correctly ordered
betw$FloodLevel <- factor(betw$FloodLevel, levels = c("Min", "Max", "Difference"))
betw$Metric     <- "Betweenness"

################################################################################
#### Preparing Human Wildlife Interaction Maps
################################################################################
# Load the distance maps
hwcm <- read_rds("03_Data/03_Results/Distance.rds")

# I also want to plot a difference map
diff <- subset(hwcm, Level == "Global")
diff <- diff$Distance[[2]] - diff$Distance[[1]]
diff <- tibble(
    Filename   = NA
  , SourceArea = NA
  , FloodLevel = "Difference"
  , Level      = "Difference"
  , Distance   = list(diff)
)
hwcm <- rbind(hwcm, diff)

# Convert rasterlayers to dataframe
hwcm$Distance <- lapply(hwcm$Distance, function(x) {
  x <- as.data.frame(x, xy = T)
  names(x) <- c("x", "y", "Value")
  return(x)
})
hwcm <- unnest(hwcm, Distance)

# Make sure the levels are correctly ordered
hwcm$Steps      <- 2000
hwcm$FloodLevel <- factor(hwcm$FloodLevel, levels = c("Min", "Max", "Difference"))
hwcm$Filename   <- NULL
hwcm$Metric     <- "HWC"

################################################################################
#### Preparing Alternative Human Wildlife Interaction Maps
################################################################################
# For the alternative method, we'll combine human influence (continuous) with
# The heatmap. Thus, we need to load the heatmaps and human influence map
hwcm_alternative <- "03_Data/03_Results/HeatmapsBetweennessGlobal.rds" %>%
  read_rds() %>%
  mutate(SourceArea = NA, Level = "Global") %>%
  select(Steps, SourceArea, FloodLevel, Level, Data = Heatmap) %>%
  arrange(Steps, SourceArea, FloodLevel, Level)
humans <- raster("03_Data/02_CleanData/HumanInfluence.tif")

# Make sure the human influence layer has the same resolution, extent, and
# origin
humans <- aggregate(humans, fact = res(hwcm_alternative$Data[[1]]) / res(humans)[1])
humans <- resample(humans, hwcm_alternative$Data[[1]], method = "ngb")

# Go through the heatmaps and compute alternative human wildlife conflict maps
hwcm_alternative$Data <- lapply(hwcm_alternative$Data, function(x) {
  hwcm_alt <- x * humans
  return(hwcm_alt)
})

# I also want to plot a difference map
diff <- subset(hwcm_alternative, Level == "Global" & Steps == 2000)
diff <- diff$Data[[1]] - diff$Data[[2]]
diff <- tibble(
    Steps      = 2000
  , SourceArea = NA
  , FloodLevel = "Difference"
  , Level      = "Difference"
  , Data       = list(diff)
)
hwcm_alternative <- rbind(hwcm_alternative, diff)

# Convert rasterlayers to dataframe
hwcm_alternative$Data <- lapply(hwcm_alternative$Data, function(x) {
  x <- as.data.frame(x, xy = T)
  names(x) <- c("x", "y", "Value")
  return(x)
})
hwcm_alternative <- unnest(hwcm_alternative, Data)

# Make sure the levels are correctly ordered
hwcm_alternative$Steps      <- 2000
hwcm_alternative$FloodLevel <- factor(hwcm_alternative$FloodLevel, levels = c("Min", "Max", "Difference"))
hwcm_alternative$Filename   <- NULL
hwcm_alternative$Metric     <- "HWC_Alternative"

################################################################################
#### Combine Metrics
################################################################################
maps <- rbind(heat, betw, hwcm, hwcm_alternative)

################################################################################
#### Functions to Plot
################################################################################
# Function to plot the different metrics (except for the inter-patch
# connectivity)
plotMet <- function(
    data
  , formula     = ~FloodLevel
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
  , aois        = F
  , diffmap     = F
  , which_aois  = c("Maun", "Lake Ngami", "Linyanti Swamp", "Panhandle")
  ) {

    # Testing
    # data        <- subset(maps, Steps == 2000 & Level == "Global" & Metric == "HWC")
    # formula     <- ~FloodLevel
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
    # aois        <- T
    # lomehi      <- T
    # cutlines    <- F
    # ext <- ext(r)

  # Some checks
  colorscheme     <- match.arg(colorscheme)
  col_area        <- match.arg(col_area)
  col_country     <- match.arg(col_country)
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
    {
      if (cutlines) {
        geom_sf(
            data      = cutli
          , col       = "black"
          , lty       = 2
          , linewidth = 0.4
        )
      }
    } +
    {
      if (cutlines) {
        geom_text(
            data     = data.frame(x = 22.3, y = -18.7, Label = "Line of Separation")
          , mapping  = aes(x = x, y = y, label = Label)
          , col      = "black"
          , fontface = 3
          , size     = 1.3
          , angle    = -45
        )
      }
    } +
    {
      if (aois) {
        geom_sf(
            data      = subset(intre, Name %in% which_aois)
          , col       = "black"
          , fill      = NA
          , lty       = 2
          , linewidth = 0.4
        )
      }
    } +
    geom_point(
        data        = subset(vills, place == "City")
      , mapping     = aes(x = x, y = y, size = place)
      , col         = ifelse(colorscheme == "dark", "gray90", "gray30")
      , shape       = 15
      , show.legend = F
      , size        = 1
    ) +
    geom_sf(
        data        = subset(areas, Type == "Main")
      , color       = col_area
      , fill        = col_area
      , lty         = 1
      , linewidth   = 0.2
      , show.legend = F
      , alpha       = 0.15
    ) +
    {
      if (!is.null(area)) {
        geom_sf(
            data      = subset(areas, ID %in% area)
          , col       = "red"
          , fill      = NA
          , linewidth = 0.5
        )
      }
    } +
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
        data     = subset(labels_areas, Type == "Main")
      , mapping  = aes(x = X, y = Y, label = ID)
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
    {
      if (aois) {
        geom_sf_text(
            data        = subset(intre, Name %in% which_aois & !(Name %in% c("Maun", "Linyanti Swamp")))
          , mapping     = aes(label = Name)
          , col         = ifelse(colorscheme == "dark", "gray90", "gray30")
          , fontface    = 3
          , size        = 2
        )
      }
    } +
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
    scale_x_continuous(breaks = seq(21.5, 25.5, by = 1)) +
    scale_y_continuous(breaks = seq(-20.5, -18.5, by = 1)) +
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

# Function to plot the inter-patch connectivity maps
plotNet <- function(
    sources     = NULL
  , targets     = NULL
  , pairs       = NULL
  , reverse     = F
  , area        = NULL
  , col_area    = "grey"
  , curvature   = 0.2
  , col_water   = NA
  , fill_water  = NA
  , alpha_water = 1
  , ext         = NULL
  , max_size    = 1
  , trans       = "identity"
  , zones       = c("Main", "Buffer")
  , barwidth    = unit(8, "cm")
  , days        = T
  ) {

  # Set extent if none is given
  if (is.null(ext)) {
    ext <- extent(r)
  }

  # Convert to days if desired
  if (days) {
    ipc$DispersalDuration   <- ipc$DispersalDuration / 5
    ipc$DispersalDurationSD <- ipc$DispersalDurationSD / 5
  }

  # Subset to desired data
  sub <- select(ipc, SourceArea, CurrentArea, DispersalSuccess, DispersalDuration, FloodLevel)
  if (!is.null(sources) & !is.null(targets)) {
    subi <- subset(sub, SourceArea %in% sources & CurrentArea %in% targets)
    if (reverse) {
      subi2 <- subset(sub, SourceArea %in% targets & CurrentArea %in% sources)
      subi <- rbind(subi, subi2)
    }
    sub <- subi
    sub <- subset(sub, SourceArea != CurrentArea)
  } else if (!is.null(pairs)) {
    all <- list()
    for (i in 1:nrow(pairs)) {
      subi <- subset(sub, SourceArea == pairs[i, 1] & CurrentArea == pairs[i, 2])
      if (reverse) {
        subi2 <- subset(sub, SourceArea == pairs[i, 2] & CurrentArea == pairs[i, 1])
        subi <- rbind(subi, subi2)
      }
      all[[i]] <- subi
    }
    sub <- do.call(rbind, all)
    sub <- subset(sub, SourceArea != CurrentArea)
  }

  # Generate a network
  net <- graph_from_data_frame(
      d        = sub
    , vertices = unique(areas$ID)
    , directed = T
  )

  # Prepare networks for ggplotting with ggplot
  net_p <- suppressWarnings(ggnetwork(net, layout = lay, arrow.gap = 0.1, scale = F))
  net_p <- subset(net_p, !is.na(FloodLevel))
  net_p$FloodLevel <- factor(net_p$FloodLevel, levels = c("Min", "Max"))

  # Plot
  ggplot() +
    geom_sf(
        data  = water
      , fill  = fill_water
      , col   = col_water
      , alpha = alpha_water
    ) +
    # geom_sf(
    #     data  = subset(areas, Type == "Main")
    #   , col   = "orange"
    #   , fill  = "orange"
    #   , alpha = 0.2
    # ) +
    # geom_sf(
    #     data  = subset(areas, Type == "Buffer")
    #   , col   = "purple"
    #   , fill  = "purple"
    #   , alpha = 0.2
    # ) +
    geom_sf(
      #         data        = subset(areas, Type == "Main")
        data        = subset(areas, Type %in% zones)
      , color       = col_area
      , fill        = col_area
      , lty         = 1
      , linewidth   = 0.2
      , show.legend = F
      , alpha       = 0.15
    ) +
    {
      if (!is.null(area)) {
        geom_sf(
            data      = subset(areas, ID %in% area)
          , col       = "red"
          , fill      = NA
          , linewidth = 0.5
        )
      }
    } +
    geom_text(
        data     = labels_waters
      , mapping  = aes(x = x, y = y, label = Label)
      , col      = darken(fill_water, 1.5)
      , fontface = 3
      , size     = 1.3
    ) +
    {
      if ("Buffer" %in% zones) {
        geom_point(
            data     = subset(labels_areas, Type == "Buffer")
          , mapping  = aes(x = X, y = Y)
          , col      = "gray"
          , size     = 4
        )
      }
    } +
    geom_edges(
        data      = net_p
      , mapping   = aes(
          x    = x
        , y    = y
        , xend = xend
        , yend = yend
        , size = DispersalSuccess
        , col  = DispersalDuration
      )
      , curvature = curvature
      , arrow     = arrow(length = unit(6, "pt"), type = "closed", angle = 10)
    ) +
    geom_text(
        data     = subset(labels_areas, Type %in% zones)
      , mapping  = aes(x = X, y = Y, label = ID)
      , col      = "black"
      , fontface = 3
      , size     = 2
    ) +
    scale_fill_manual(
        values = c("transparent", "cornflowerblue")
      , name   = ""
    ) +
    facet_wrap(~ FloodLevel) +
    coord_sf(
        crs    = 4326
      , xlim   = ext[1:2]
      , ylim   = ext[3:4]
      , expand = F
    ) +
    scale_size_area(
        name     = "Dispersal Frequency\n(number of trajectories)"
      , max_size = max_size
      , breaks   = c(200, 400, 600)
      , labels   = c(200, 400, 600)
      , guide    = guide_legend(
          order          = 2
        , title.hjust    = 0.5
        , title.position = "bottom"
        # , label.position = "bottom"
      )
    ) +
    scale_color_gradientn(
        colors = rev(brewer.pal(11, "RdYlGn"))
      , trans  = trans
      , guide  = guide_colorbar(
        , title          = ifelse(days, "Dispersal Duration\n(days)", "Dispersal Duration\n(steps)")
        , show.limits    = T
        , title.position = "bottom"
        , title.hjust    = 0.5
        , ticks          = F
        , barheight      = unit(0.2, "cm")
        , barwidth       = barwidth
      )
    ) +
    scale_x_continuous(breaks = seq(21.5, 25.5, by = 1)) +
    scale_y_continuous(breaks = seq(-20.5, -18.5, by = 1)) +
    theme(
        panel.background = element_blank()
      , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
      , legend.position  = "bottom"
      , legend.box       = "horizontal"
      , legend.key       = element_blank()
      , strip.background = element_rect(fill = "gray95", color = "transparent")
    ) +
    annotation_scale(
        location   = "bl"
      , width_hint = 0.2
      , line_width = 0.5
      , text_cex   = 0.5
      , height     = unit(0.1, "cm")
    ) +
    annotation_north_arrow(
        location = "br"
      , height   = unit(0.7, "cm"),
      , width    = unit(0.6, "cm"),
      , style    = north_arrow_fancy_orienteering(
            fill      = "black"
          , line_col  = NA
          , text_col  = "black"
          , text_size = 4
        )
    ) +
    xlab("") +
    ylab("")
}

################################################################################
#### Connectivity Plots
################################################################################
# Define plot extent
ext <- c(21, 25.1, -20.7, -17.9)

# Create plot of interpatch connectivity between selected areas (neighboring
# ones)
p <- plotNet(
    pairs = rbind(
        c(1, 2)
      , c(1, 5)
      , c(1, 6)
      , c(2, 3)
      , c(2, 6)
      , c(3, 4)
      , c(3, 6)
      , c(4, 5)
      , c(4, 6)
      , c(5, 6)
    )
  , barwidth    = unit(10, "cm")
  , reverse     = T
  , fill_water  = "cornflowerblue"
  , alpha_water = 1
  , ext         = ext
  , max_size    = 1.5
  , zones       = c("Main", "Buffer")
)

# Store
ggsave("04_Manuscript/Figures/InterpatchConnectivityMap.png"
  , plot   = p
  , width  = 5
  , height = 3
  , scale  = 1.5
  , bg     = "white"
  , device = png
)

# Generate Heatmaps
p1 <- plotMet(subset(maps, Level == "Global" & Steps == 2000 & Metric == "Heatmap")
  , barwidth    = unit(8, "cm")
  , colorscheme = "dark"
  , col_area    = "black"
  , col_country = "black"
  , col_water   = NA
  , fill_water  = NA
  , palette     = spectral(100)
  , name        = "# Traversing Trajectories"
  , trans       = "identity"
  , ext         = ext
  , lomehi      = F
  , aois        = F
)

# Generate Betweenness maps
p2 <- plotMet(subset(maps, Level == "Global" & Steps == 2000 & Metric == "Betweenness")
  , barwidth    = unit(8, "cm")
  , colorscheme = "dark"
  , col_area    = "white"
  , col_country = "white"
  , col_water   = NA
  , fill_water  = "gray80"
  , palette     = magma(100)
  , name        = "Betweenness"
  , trans       = "sqrt"
  , ext         = ext
  , lomehi      = F
  , aois        = F
)

# Generate human wildlife conflict maps
p3 <- plotMet(subset(maps, Level == "Global" & Steps == 2000 & Metric == "HWC")
  , barwidth    = unit(8, "cm")
  , colorscheme = "light"
  , col_area    = "black"
  , col_country = "black"
  , col_water   = NA
  , fill_water  = "gray60"
  , palette     = terrain.colors(100, rev = T)
  , name        = "Human-Wildlife Conflict Potential"
  , trans       = "identity"
  , ext         = ext
  , lomehi      = T
  , aois        = T
  , which_aois  = c("Panhandle", "Maun")
)

# Generate alternative human wildlife conflict maps
p3_appendix <- plotMet(subset(maps, Level == "Global" & Steps == 2000 & Metric == "HWC_Alternative")
  , barwidth    = unit(8, "cm")
  , colorscheme = "light"
  , col_area    = "black"
  , col_country = "black"
  , col_water   = NA
  , fill_water  = "gray60"
  , palette     = terrain.colors(100, rev = T)
  , name        = "Human-Wildlife Conflict Potential"
  , trans       = "identity"
  , ext         = ext
  , lomehi      = T
  , aois        = T
  , which_aois  = c("Panhandle", "Maun")
)

################################################################################
#### Difference Maps
################################################################################
# We'll use the maximum water extent for now
water_backup <- water
water        <- subset(water, FloodLevel == "Max")

# Create difference map for heatmaps
p4 <- plotMet(subset(maps, Level == "Difference" & Metric == "Heatmap")
  , barwidth    = unit(6, "cm")
  , formula     = ~"Max - Min"
  , colorscheme = "light"
  , col_area    = "black"
  , col_country = "black"
  , col_water   = NA
  , fill_water  = NA
  , palette     = colorRampPalette(c("orange", "white", "cornflowerblue"))(100)
  , name        = "Difference"
  , trans       = "identity"
  , ext         = ext
  , lomehi      = T
  , diffmap     = T
  , aoi         = F
  , limits      = c(-350, 350)
) + rremove("y.text") + rremove("y.ticks")

# Create difference map for betweenness maps
p5 <- plotMet(subset(maps, Level == "Difference" & Metric == "Betweenness")
  , barwidth    = unit(6, "cm")
  , formula     = ~"Max - Min"
  , colorscheme = "light"
  , col_area    = "black"
  , col_country = "black"
  , col_water   = NA
  , fill_water  = "gray60"
  , palette     = colorRampPalette(c("orange", "white", "cornflowerblue"))(100)
  , name        = "Difference"
  , trans       = "identity"
  , ext         = ext
  , lomehi      = T
  , diffmap     = T
  , aoi         = F
  , limits      = c(-5e8, 5e8)
) + rremove("y.text") + rremove("y.ticks")

# Create difference map for hwc maps
p6 <- plotMet(subset(maps, Level == "Difference" & Metric == "HWC")
  , barwidth    = unit(6, "cm")
  , formula     = ~"Max - Min"
  , colorscheme = "light"
  , col_area    = "black"
  , col_country = "black"
  , col_water   = NA
  , fill_water  = "gray60"
  , palette     = colorRampPalette(c("orange", "white", "cornflowerblue"))(100)
  , name        = "Difference"
  , trans       = "identity"
  , ext         = ext
  , lomehi      = T
  , diffmap     = T
  , aoi         = T
  , which_aois  = c("Panhandle", "Maun")
  , limits      = c(-100, 100)
) + rremove("y.text") + rremove("y.ticks")

# Create difference map for alternative hwc maps
p6_appendix <- plotMet(subset(maps, Level == "Difference" & Metric == "HWC_Alternative")
  , barwidth    = unit(6, "cm")
  , formula     = ~"Max - Min"
  , colorscheme = "light"
  , col_area    = "black"
  , col_country = "black"
  , col_water   = NA
  , fill_water  = "gray60"
  , palette     = colorRampPalette(c("orange", "white", "cornflowerblue"))(100)
  , name        = "Difference"
  , trans       = "identity"
  , ext         = ext
  , lomehi      = T
  , diffmap     = T
  , aoi         = T
  , which_aois  = c("Panhandle", "Maun")
  , limits      = c(-750, 750)
) + rremove("y.text") + rremove("y.ticks")

# Also create a map of human influence only
humans <- "03_Data/02_CleanData/HumanInfluence.tif" %>%
  rast() %>%
  aggregate(fact = 4, fun = "sum") %>%
  as.data.frame(xy = T) %>%
  setNames(c("x", "y", "Value"))
p_humans <- plotMet(test
  , barwidth    = unit(6, "cm")
  , formula     = ~"Human Influence"
  , colorscheme = "dark"
  , col_area    = "white"
  , col_country = "white"
  , col_water   = NA
  , fill_water  = "gray80"
  , palette     = hcl.colors(100)
  , name        = "Human Influence"
  , trans       = "identity"
  , ext         = ext
  , lomehi      = T
  , aois        = F
  , which_aois  = c("Panhandle", "Maun")
) + rremove("y.text") + rremove("y.ticks")

# Arrange all
p7 <- ggarrange(p1, p4, widths = c(2, 1), labels = c("a1", "a2"), label.x = c(0, -0.05))
p8 <- ggarrange(p2, p5, widths = c(2, 1), labels = c("b1", "b2"), label.x = c(0, -0.05))
p9 <- ggarrange(p3, p6, widths = c(2, 1), labels = c("c1", "c2"), label.x = c(0, -0.05))
p  <- ggarrange(p7, p8, p9, nrow = 3)

# Also prepare a HWC plot for the appendixe
p1_appendix <- ggarrange(p1, p_humans, widths = c(2, 1), labels = c("a", "b"), label.x = c(0, -0.05))
p2_appendix <- ggarrange(p3_appendix, p6_appendix, widths = c(2, 1), labels = c("c1", "c2"), label.x = c(0, -0.05))
p_appendix  <- ggarrange(p1_appendix, p2_appendix, nrow = 2)

# Store
ggsave("04_Manuscript/Figures/Metrics.png"
  , plot   = p
  , width  = 5
  , height = 5
  , scale  = 2
  , bg     = "white"
  , device = png
)

# Store
ggsave("04_Manuscript/Figures/HumanWildlifeConflictAlternative.png"
  , plot   = p_appendix
  , width  = 5
  , height = 5 / 3 * 2
  , scale  = 2
  , bg     = "white"
  , device = png
)

################################################################################
#### Source Specific Maps: Heatmaps
################################################################################
# Use min-max, watermaps again
water <- water_backup

# I also want to generate source-area specific plots
ids <- na.omit(unique(maps$SourceArea))
p   <- list()
for (i in ids) {
  p[[i]] <- plotMet(subset(maps, Level == "Local" & Steps == 2000 & Metric == "Heatmap" & SourceArea == i)
    , barwidth    = unit(10, "cm")
    , colorscheme = "dark"
    , col_area    = "black"
    , col_country = "black"
    , col_water   = NA
    , fill_water  = NA
    , palette     = spectral(100)
    , name        = "# Traversing Trajectories"
    , trans       = "identity"
    , ext         = ext
    , lomehi      = F
    , area        = i
  )
}

# Put the plots together
p1 <- ggarrange(p[[1]], p[[2]], p[[3]], ncol = 1, labels = c("a", "b", "c"))
p2 <- ggarrange(p[[4]], p[[5]], p[[6]], ncol = 1, labels = c("d", "e", "f"))
p3 <- ggarrange(p1, p2, ncol = 2)

# Store the plot
ggsave("04_Manuscript/Figures/HeatmapsIndividual.png"
  , plot   = p3
  , bg     = "white"
  , height = 7
  , width  = 8
  , scale  = 1.4
  , device = png
)

################################################################################
#### Source Specific Maps: Betweenness
################################################################################
# I also want to generate source-area specific plots
ids <- na.omit(unique(maps$SourceArea))
p   <- list()
for (i in ids) {
  p[[i]] <- plotMet(subset(maps, Level == "Local" & Steps == 2000 & Metric == "Betweenness" & SourceArea == i)
    , barwidth    = unit(10, "cm")
    , colorscheme = "dark"
    , col_area    = "white"
    , col_country = "white"
    , col_water   = NA
    , fill_water  = "gray80"
    , palette     = magma(100)
    , name        = "Betweenness"
    , trans       = "sqrt"
    , ext         = ext
    , lomehi      = F
    , area        = i
  )
}

# Put the plots together
p1 <- ggarrange(p[[1]], p[[2]], p[[3]], ncol = 1, labels = c("a", "b", "c"))
p2 <- ggarrange(p[[4]], p[[5]], p[[6]], ncol = 1, labels = c("d", "e", "f"))
p3 <- ggarrange(p1, p2, ncol = 2)

# Store the plot
ggsave("04_Manuscript/Figures/BetweennessIndividual.png"
  , plot   = p3
  , bg     = "white"
  , height = 7
  , width  = 8
  , scale  = 1.4
  , device = png
)

################################################################################
#### Source Specific Maps: Human Wildlife Interaction
################################################################################
# I also want to generate source-area specific plots
ids <- na.omit(unique(maps$SourceArea))
p   <- list()
for (i in ids) {
  p[[i]] <- plotMet(subset(maps, Level == "Local" & Steps == 2000 & Metric == "HWC" & SourceArea == i)
    , barwidth    = unit(10, "cm")
    , colorscheme = "light"
    , col_area    = "black"
    , col_country = "black"
    , col_water   = NA
    , fill_water  = "gray60"
    , palette     = terrain.colors(100, rev = T)
    , name        = "Human-Wildlife Conflict Potential"
    , trans       = "identity"
    , ext         = ext
    , lomehi      = T
    , area        = i
    , aois        = T
  )
}

# Put the plots together
p1 <- ggarrange(p[[1]], p[[2]], p[[3]], ncol = 1, labels = c("a", "b", "c"))
p2 <- ggarrange(p[[4]], p[[5]], p[[6]], ncol = 1, labels = c("d", "e", "f"))
p3 <- ggarrange(p1, p2, ncol = 2)

# Store the plot
ggsave("04_Manuscript/Figures/HumanWildlifeConflictIndividual.png"
  , plot   = p3
  , bg     = "white"
  , height = 7
  , width  = 8
  , scale  = 1.4
  , device = png
)

################################################################################
#### Source Specific Maps: Inter-Patch Connectivity
################################################################################
# I also want to generate source-area specific plots
ids <- na.omit(unique(maps$SourceArea))

# Let's generate the same plot for any other source area
p <- list()
for (i in 1:6) {
  p[[i]] <- plotNet(i, base::setdiff(1:6, i)
    , barwidth    = unit(5, "cm")
    , reverse     = T
    , fill_water  = "cornflowerblue"
    , alpha_water = 1
    , ext         = ext
    , zones       = c("Main", "Buffer")
    , area        = i
  )
}

# Put the plots together
p1 <- ggarrange(p[[1]], p[[2]], p[[3]], ncol = 1, labels = c("a", "b", "c"))
p2 <- ggarrange(p[[4]], p[[5]], p[[6]], ncol = 1, labels = c("d", "e", "f"))
p3 <- ggarrange(p1, p2, ncol = 2)

# And repeat the same for the borders
p <- list()
for (i in 1:6) {
  p[[i]] <- plotNet(i, 7:14
    , barwidth    = unit(5, "cm")
    , reverse     = T
    , fill_water  = "cornflowerblue"
    , alpha_water = 1
    , ext         = ext
    , curvature   = 0
    , area        = i
  )
}

# Put the plots together
p4 <- ggarrange(p[[1]], p[[2]], p[[3]], ncol = 1, labels = c("a", "b", "c"))
p5 <- ggarrange(p[[4]], p[[5]], p[[6]], ncol = 1, labels = c("d", "e", "f"))
p6 <- ggarrange(p4, p5, ncol = 2)

# Store the plot
ggsave("04_Manuscript/Figures/IPCMain.png"
  , plot   = p3
  , bg     = "white"
  , height = 6
  , width  = 8
  , scale  = 1.8
  , device = png
)
ggsave("04_Manuscript/Figures/IPCBuffer.png"
  , plot   = p6
  , bg     = "white"
  , height = 6
  , width  = 8
  , scale  = 1.8
  , device = png
)
