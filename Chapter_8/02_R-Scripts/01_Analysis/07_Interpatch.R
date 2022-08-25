################################################################################
#### Interpatch Connectivity
################################################################################
# Description: Computing interpatch connectivity from simulated dispersal paths

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(pbmcapply)      # To run stuff in parallel
library(igraph)         # For network analysis
library(rgeos)          # To manipulate spatial objects
library(ggnetwork)      # To plot network using ggplot
library(sf)             # To plot spatial features
library(ggspatial)      # To add scale bars etc to plots

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load reference raster
r <- rast("03_Data/02_CleanData/ReferenceRaster.tif")

# Load source areas
area <- vect("03_Data/02_CleanData/SourceAreas.shp")

# Rasterize them to the reference raster
area_r <- terra::rasterize(area, y = r, field = "ID")

# Visualize them
plot(area_r, main = "Source Areas")
plot(area, add = T)
text(area, "ID", cex = 0.5, halo = T)

################################################################################
#### Prepare Simulations
################################################################################
# Load dispersal simulations
sims <- read_rds("03_Data/03_Results/DispersalSimulation.rds")
sims <- subset(sims, FloodLevel != "Mean")

# Keep only desired columns
sims <- sims[, c("x", "y", "TrackID", "StepNumber", "SourceArea", "FloodLevel")]

# # Subset for now
# sims <- subset(sims, TrackID %in% sample(sims$TrackID, 1000))

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# Identify through which national parks the dispersers moved
visits <- data.frame(
    TrackID    = sims$TrackID
  , StepNumber = sims$StepNumber
  , FloodLevel = sims$FloodLevel
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , Area       = raster::extract(raster(area_r), sims)
)

# Calculate for each step the distance to the first coordinate. We'll use this
# to determine how far an individual had to disperse before reaching another
# area
visits <- visits %>%
  nest(data = -TrackID) %>%
  mutate(data = pbmclapply(data
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x) {

      # Project coordinates
      coords <- reprojCoords(
          xy   = x[, c("x", "y")]
        , from = CRS("+init=epsg:4326")
        , to   = CRS("+init=epsg:32734")
      )

      # Compute distance to first coordinate
      first <- coords[1, ]
      distance <- sqrt((coords[, 1] - first[1]) ** 2 + (coords[, 2] - first[2]) ** 2)
      x$DistanceFromFirst <- distance

      # Return the resulting object
      return(x)
  })) %>%
  unnest(data)

# Add the information on the reached parks and distance traveled to the
# simulations
sims$CurrentArea       <- visits$Area
sims$DistanceFromFirst <- visits$DistanceFromFirst

# Convert simulations to regular dataframe
sims <- as.data.frame(sims, xy = T)
sims$xy <- NULL

# Identify how long it takes to reach the different areas
visits <- sims %>%
  rename(From = SourceArea, To = CurrentArea) %>%
  group_by(TrackID, FloodLevel, From, To) %>%
  summarize(
      StepNumber        = min(StepNumber)
    , DistanceFromFirst = min(DistanceFromFirst)
    , .groups           = "drop"
  ) %>%
  subset(!is.na(From) & !is.nan(To) & !is.na(To)) %>%
  arrange(TrackID, StepNumber)

# # Replace national park IDs with proper park names. Also identify the country
# # for each of the national parks.
# visits$FromArea <- as.character(area$ID[match(visits$From, area$ID)])
# visits$ToArea   <- as.character(area$ID[match(visits$To, area$ID)])

# Summarize connections by area. Compute average number of steps required to
# make connection, it's standard deviation, as well as how often the connection
# has been made.
visits_area <- visits %>%
  group_by(FloodLevel, From, To) %>%
  summarize(
      MeanStepNumber = mean(StepNumber)
    , SDStepNumber   = sd(StepNumber)
    , Frequency      = n()
    , .groups        = "drop"
  )

# Store visits to file
write_rds(visits, "03_Data/03_Results/InterpatchConnectivity.rds")

# # Remove self-loops and zero-links
# visits <- subset(visits_area, From != To & Frequency != 0)
#
# # Create a network
# lay <- coordinates(gCentroidWithin(as(area, "Spatial")))
# net_p <- lapply(unique(visits$FloodLevel), function(x) {
#
#   # Make network
#   net <- graph_from_data_frame(
#       d        = dplyr::select(subset(visits, FloodLevel == x), From, To, Frequency, everything())
#     , vertices = unique(area$ID)
#     , directed = T
#   )
#
#   # Prepare networks for ggplotting with ggplot
#   net <- ggnetwork(net, layout = lay, arrow.gap = 0.1, scale = F)
#
#   # Remove NA entries
#   net <- subset(net, !is.na(FloodLevel))
#
#   # Return it
#   return(net)
#
# }) %>% do.call(rbind, .)
#
# # Main Plot
# p1 <- ggplot() +
#   geom_sf(
#       data        = st_as_sf(area)
#     , col         = "#6BA36B"
#     , fill        = "darkgreen"
#     , alpha       = 0.5
#     , lwd         = 0
#     , show.legend = F
#   ) +
#   geom_edges(
#       data      = net_p
#     , mapping   = aes(
#         x    = x
#       , y    = y
#       , xend = xend
#       , yend = yend
#       , size = Frequency
#       , col  = MeanStepNumber
#     )
#     , curvature = 0.2
#     , arrow     = arrow(length = unit(6, "pt"), type = "closed", angle = 10)
#   ) +
#   geom_sf_text(
#       data     = st_as_sf(area)
#     , mapping  = aes(label = ID)
#     , col      = "black"
#     , fontface = 3
#     , size     = 3
#   ) +
#   scale_size_area(
#       name     = "Relative Frequency"
#     , max_size = 1
#   ) +
#   scale_color_gradientn(
#       colors  = viridis::magma(50)
#     , guide   = guide_colorbar(
#         title          = "Duration (Steps)"
#       , show.limits    = T
#       , title.position = "top"
#       , title.hjust    = 0.5
#       , ticks          = T
#       , barheight      = unit(0.6, "cm")
#       , barwidth       = unit(3.0, "cm")
#       , order = 1
#     )
#   ) +
#   # scale_fill_manual(
#   #   values = c("#70ab70", "#d9f0d3")
#   # ) +
#   # coord_sf(
#   #     crs    = 4326
#   #   , xlim   = c(min(r$x), max(r$x))
#   #   , ylim   = c(min(r$y), max(r$y))
#   #   , expand = F
#   # ) +
#   labs(
#       x        = NULL
#     , y        = NULL
#     , fill     = NULL
#     # , title    = "Interpatch Connectivity"
#     # , subtitle = "In Relation to Dispersal Duration"
#   ) +
#   guides(
#       size  = guide_legend(title.position = "top", order = 2)
#   ) +
#   theme(
#       legend.position      = "bottom"
#     , legend.box           = "horizontal"
#     , legend.title.align   = 0.5
#     , panel.background     = element_blank()
#     , panel.border         = element_rect(colour = "black", fill = NA, size = 1)
#     , legend.title         = element_text(size = 10),
#     , legend.text          = element_text(size = 8)
#     , legend.margin        = margin(c(0, 0, 0, 0))
#     , legend.key           = element_blank()
#   ) +
#   annotation_scale(
#       location   = "tl"
#     , width_hint = 0.2
#     , line_width = 1
#     , height     = unit(0.15, "cm")
#     , bar_cols   = c("black", "white")
#     , text_col   = "black"
#   ) +
#   annotation_north_arrow(
#       location = "tr"
#     , height   = unit(1.5, "cm"),
#     , width    = unit(1.2, "cm"),
#     , style    = north_arrow_fancy_orienteering(
#           fill      = c("black", "black")
#         , line_col  = NA
#         , text_col  = "black"
#         , text_size = 12
#       )
#   ) +
#   facet_wrap(~ FloodLevel, ncol = 1)
#
# # Store the plot
# ggsave("04_Manuscript/Interpatch.png", plot = p, width = 8, height = 6)
