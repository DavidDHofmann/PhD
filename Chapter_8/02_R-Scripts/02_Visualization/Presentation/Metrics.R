################################################################################
#### Metrics
################################################################################
# Description: Visualization of all spatial metrics

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
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
#### Some General Metrics
################################################################################
# Load data needed for this
hum <- rast("03_Data/02_CleanData/HumanDensity.tif")
agr <- rast("03_Data/02_CleanData/Agriculture.tif")
cat <- rast("03_Data/02_CleanData/CattleDensity.tif")
aoi <- vect("03_Data/02_CleanData/AreasOfInterest.shp")[c(1, 2, 4), ]
aoi <- smooth(aoi)

# Let's compute the human density per km2 and the cattle density within the
# different polygons as well as the cover by agricultural fields
aoi$HumanDensity        <- terra::extract(hum, aoi, fun = "sum")[, 2]
aoi$CattleDensityPerKM2 <- terra::extract(cat, aoi, fun = "mean")[, 2]
aoi$Agriculture         <- terra::extract(agr, aoi, fun = "mean")[, 2] * 100
aoi$Area                <- expanse(aoi, unit = "km")
aoi$HumanDensityPerKM2  <- aoi$HumanDensity / aoi$Area
aoi <- values(aoi)

# Store the relevant metrics
writeLines(as.character(round(aoi$HumanDensityPerKM2[2])), "04_Manuscript/99_HumanDensityMaun.tex")
writeLines(as.character(round(aoi$HumanDensityPerKM2[3])), "04_Manuscript/99_HumanDensityPanhandle.tex")
writeLines(as.character(round(aoi$CattleDensityPerKM2[2])), "04_Manuscript/99_CattleDensityMaun.tex")
writeLines(as.character(round(aoi$CattleDensityPerKM2[3])), "04_Manuscript/99_CattleDensityPanhandle.tex")
writeLines(as.character(round(aoi$Agriculture[2])), "04_Manuscript/99_AgricultureMaun.tex")
writeLines(as.character(round(aoi$Agriculture[3])), "04_Manuscript/99_AgriculturePanhandle.tex")

################################################################################
#### Load Data
################################################################################
# Load shapefiles that we want to plot
areas <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
intre <- read_sf("03_Data/02_CleanData/AreasOfInterest.shp")[c(1, 2, 4), ]
intre <- smooth(intre)
roads <- read_sf("03_Data/02_CleanData/Roads.shp")
afric <- read_sf("03_Data/02_CleanData/Africa.shp")
vills <- read_sf("03_Data/02_CleanData/Villages.shp")
cutli <- read_sf("03_Data/02_CleanData/Cutlines.shp")[c(4, 8), ]
cutli <- st_crop(cutli, st_bbox(c(xmin = 21.5, xmax = 23.6, ymin = -20.1, ymax = -18.2)))
vills <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y)

# Load results on interpatch connectivity
ipc <- read_rds("03_Data/03_Results/BootstrappedInterpatchConnectivity.rds")

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

# test1 <- diff$Data[[1]]
# test2 <- diff$Data[[2]]
# var(test1[])
# var(test2[])
# Moran(test1, direc)
# par(mfrow = c(1, 2))
# install.packages("usdm")
# library(usdm)
# test1_vario <- Variogram(test1, lag = 0.1)
# test2_vario <- Variogram(test2, lag = 0.1)
# plot(NA, xlim = c(0, 2.5), ylim = c(0, 30000))
# lines(gamma ~ distance, test1_vario@variogram, type = "o")
# lines(gamma ~ distance, test2_vario@variogram, type = "o")
# library(pgirmess)
# correlo1 <- correlog(coords = coordinates(test1), z = test1[], nbclass = 12)
# correlo2 <- correlog(coords = coordinates(test1), z = test1[], nbclass = 10)

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
#### Preparing Human Wildlife Conflict Maps
################################################################################
# Load the distance maps
maps1 <- read_rds("03_Data/03_Results/Distance.rds")
maps1$Level = "Local"

# Create "Global Metrics"
maps2 <- lapply(c("Min", "Max"), function(x) {
  global <- maps1 %>%
    subset(FloodLevel == x) %>%
    pull(Distance) %>%
    stack() %>%
    rast() %>%
    sum()
  global <- tibble(
      FloodLevel = x
    , SourceArea = NA
    , Filename   = NA
    , Distance   = list(raster(global))
    , Level      = "Global"
  )
  return(global)
}) %>% do.call(rbind, .)

# Put all maps together
hwcm <- rbind(maps1, maps2)

# Smoothen them
hwcm$Distance <- lapply(hwcm$Distance, function(x) {
  x <- rast(x)
  x <- aggregate(x, fact = 10, fun = "sum")
  x <- focal(x, c(3, 3), fun = mean)
  return(x)
})

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
#### Combine Metrics
################################################################################
maps <- rbind(heat, betw, hwcm)

################################################################################
#### Functions to Plot
################################################################################
# Function to plot the different metrics (except for the inter-patch
# connectivity)
plotMet <- function(
    data
  , formula     = ~FloodLevel
  , area        = unique(areas$ID)
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
  ) {

  # Testing
  # data        <- subset(maps, Steps == 2000 & Level == "Global" & Metric == "Heatmap")
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
  # lomehi      <- F
  # cutlines    <- F
  # ext <- ext(r)

  # Some checks
  areas$Highlight <- areas$ID %in% area
  colorscheme     <- match.arg(colorscheme)
  col_area        <- match.arg(col_area)
  col_country     <- match.arg(col_country)
  breaks          <- c(
      min(data$Value)
    , (min(data$Value) + max(data$Value)) / 2
    , max(data$Value)
  )
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
            data      = intre
          , col       = "orange"
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
      , mapping     = aes(col = Highlight)
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
            data        = intre[-2, ]
          , mapping     = aes(label = Name)
          , col         = ifelse(colorscheme == "dark", "gray90", "gray30")
          , fontface = 3
          , size     = 2
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
          , labels  = c("Low", "Medium", "High")
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
        legend.position   = "bottom"
      , legend.box        = "vertical"
      , legend.background = element_blank()
      , legend.text       = element_text(color = "white")
      , legend.title      = element_text(color = "white")
      , legend.key        = element_blank()
      , panel.background  = element_blank()
      , panel.border      = element_blank()
      , panel.grid.minor  = element_blank()
      , panel.grid.major  = element_blank()
      , plot.background   = element_blank()
      , axis.text         = element_text(color = "white")
      , axis.ticks        = element_line(color = "white")
      , strip.background  = element_blank()
      , strip.text.x      = element_blank()
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
  , curvature   = 0.2
  , col_water   = NA
  , fill_water  = NA
  , alpha_water = 1
  , ext         = NULL
  , max_size    = 1
  , trans       = "identity"
  , zones       = c("Main", "Buffer")
  , barwidth    = unit(8, "cm")
  ) {

  # Set extent if none is given
  if (is.null(ext)) {
    ext <- extent(r)
  }

  # Subset to desired data
  sub <- select(ipc, SourceArea, CurrentArea, Freq, StepNumber, FloodLevel)
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
        data  = subset(areas, Type %in% zones)
      , col   = "gray"
      , fill  = "gray"
      , alpha = 0.2
    ) +
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
        , size = Freq
        , col  = StepNumber
      )
      , curvature = curvature
      , arrow     = arrow(length = unit(6, "pt"), type = "closed", angle = 10)
    ) +
    geom_text(
        data     = subset(labels_areas, Type %in% zones)
      , mapping  = aes(x = X, y = Y, label = ID)
      , col      = "white"
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
        name     = "Frequency"
      , max_size = max_size
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
        , title          = "Duration (in steps)"
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
        legend.position   = "bottom"
      , legend.box        = "vertical"
      , legend.background = element_blank()
      , legend.text       = element_text(color = "white")
      , legend.title      = element_text(color = "white")
      , legend.key        = element_blank()
      , panel.background  = element_blank()
      , panel.border      = element_blank()
      , panel.grid.minor  = element_blank()
      , panel.grid.major  = element_blank()
      , plot.background   = element_blank()
      , axis.text         = element_blank()
      , axis.ticks        = element_blank()
      , strip.background  = element_blank()
      , strip.text.x      = element_blank()
    ) +
    xlab("") +
    ylab("")
}

################################################################################
#### Global Plots
################################################################################
# Define plot extent
ext <- c(21, 25.1, -20.7, -17.9)

# Generate Heatmaps
p1 <- plotMet(subset(maps, Level == "Global" & Steps == 2000 & Metric == "Heatmap")
  , barwidth    = unit(16, "cm")
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
  , barwidth    = unit(16, "cm")
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
  , barwidth    = unit(16, "cm")
  , colorscheme = "dark"
  , col_area    = "white"
  , col_country = "white"
  , col_water   = NA
  , fill_water  = "gray80"
  , palette     = viridis::plasma(100)
  , name        = "Human Wildlife Conflict Potential"
  , trans       = "identity"
  , ext         = ext
  , lomehi      = T
  , aois        = T
)

# Create plot of interpatch connectivity between selected source areas
p4 <- plotNet(
    pairs = rbind(
        c(1, 2)
      , c(1, 5)
      , c(5, 4)
      , c(4, 3)
      , c(3, 2)
      , c(1, 6)
      , c(2, 6)
      , c(3, 6)
      , c(4, 6)
      , c(5, 6)
      # , c(1, 7)
      # , c(1, 8)
      # , c(2, 9)
      # , c(2, 10)
      # , c(3, 10)
      # , c(3, 11)
      # , c(4, 11)
      # , c(4, 12)
      # , c(5, 13)
      # , c(5, 14)
    )
  , barwidth    = unit(10, "cm")
  , reverse     = T
  , fill_water  = "cornflowerblue"
  , alpha_water = 1
  , ext         = ext
  , max_size    = 1.5
  , zones       = "Main"
)

# Store them
ggsave("05_Presentation/99_Metrics1.png"
  , plot   = p1
  , bg     = "transparent"
  , width  = 8
  , height = 4
  , scale  = 1.1
)
ggsave("05_Presentation/99_Metrics2.png"
  , plot   = p2
  , bg     = "transparent"
  , width  = 8
  , height = 4
  , scale  = 1.1
)
ggsave("05_Presentation/99_Metrics3.png"
  , plot   = p3
  , bg     = "transparent"
  , width  = 8
  , height = 4
  , scale  = 1.1
)
ggsave("05_Presentation/99_Metrics4.png"
  , plot   = p4
  , bg     = "transparent"
  , width  = 8
  , height = 4
  , scale  = 1.1
)

################################################################################
#### Source Specific Maps: Heatmaps
################################################################################
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

# Store them
for (i in 1:length(p)) {
  plotname <- paste0("05_Presentation/99_HeatmapsIndividual_", i, ".png")
  ggsave(plotname
    , plot   = p[[i]]
    , bg     = "transparent"
    , width  = 8
    , height = 4
    , scale  = 1.1
  )
}

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

# Store them
for (i in 1:length(p)) {
  plotname <- paste0("05_Presentation/99_BetweennessIndividual_", i, ".png")
  ggsave(plotname
    , plot   = p[[i]]
    , bg     = "transparent"
    , width  = 8
    , height = 4
    , scale  = 1.1
  )
}

################################################################################
#### Source Specific Maps: Human Wildlife Conflict
################################################################################
# I also want to generate source-area specific plots
ids <- na.omit(unique(maps$SourceArea))
p   <- list()
for (i in ids) {
  p[[i]] <- plotMet(subset(maps, Level == "Local" & Steps == 2000 & Metric == "HWC" & SourceArea == i)
    , barwidth    = unit(10, "cm")
    , colorscheme = "dark"
    , col_area    = "white"
    , col_country = "white"
    , col_water   = NA
    , fill_water  = "gray80"
    , palette     = viridis::plasma(100)
    , name        = "Human Wildlife Conflict Potential"
    , trans       = "identity"
    , ext         = ext
    , lomehi      = T
    , aois        = T
    , area        = i
  )
}

# Store them
for (i in 1:length(p)) {
  plotname <- paste0("05_Presentation/99_HumanWildlifeConflictIndividual_", i, ".png")
  ggsave(plotname
    , plot   = p[[i]]
    , bg     = "transparent"
    , width  = 8
    , height = 4
    , scale  = 1.1
  )
}

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
  )
}

# Store them
for (i in 1:length(p)) {
  plotname <- paste0("05_Presentation/99_IPCMainIndividual_", i, ".png")
  ggsave(plotname
    , plot   = p[[i]]
    , bg     = "transparent"
    , width  = 8
    , height = 4
    , scale  = 1.1
  )
}

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
  )
}

# Store them
for (i in 1:length(p)) {
  plotname <- paste0("05_Presentation/99_IPCMainBuffer_", i, ".png")
  ggsave(plotname
    , plot   = p[[i]]
    , bg     = "transparent"
    , width  = 8
    , height = 4
    , scale  = 1.1
  )
}
