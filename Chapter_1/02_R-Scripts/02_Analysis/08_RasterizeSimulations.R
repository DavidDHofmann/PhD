################################################################################
#### Rasterization of Simulated Dispersal Trajectories
################################################################################
# Description: In this script, we rasterize the simulated dispersal
# trajectories and create "heatmaps".

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(terra)        # For quick raster manipulation
library(raster)       # For general raster manipulation
library(rgdal)        # To read spatial data
library(tidyverse)    # For data wrangling
library(rgeos)        # For manipulating vector data
library(lubridate)    # For working with dates
library(viridis)      # For nicer colors
library(pbmcapply)    # To show progress bar in mclapply calls
library(tmap)         # For nice spatial plots
library(davidoff)     # Custom functions
library(spatstat)     # For quick rasterization
library(maptools)     # For quick rasterization

################################################################################
#### Load and Prepare Data
################################################################################
# Load the reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Load the simulated dispersal trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Ungroup them
sims <- ungroup(sims)

# Remove undesired columns
sims <- sims[, c("x", "y", "TrackID", "StepNumber", "Area")]

# Reproject coordinates to utm (required for spatstat)
sims[, c("x", "y")] <- reprojCoords(
    xy   = sims[, c("x", "y")]
  , from = CRS("+init=epsg:4326")
  , to   = CRS("+init=epsg:32734")
)

# Prepare extent that encompassess all coordinates + some buffer
ext <- extent(min(sims$x), max(sims$x), min(sims$y), max(sims$y)) +
  c(-1000, +1000, -1000, +1000)

# Span a raster with desired resolution
r <- raster(ext, res = 1000)
values(r) <- runif(ncell(r))
crs(r) <- CRS("+init=epsg:32734")

# Collect garbage
gc()

# Check out the number of rows
nrow(sims) / 1e6

################################################################################
#### Function to Rasterize Tracks
################################################################################
# Function to rasterize trajectories after desired number of steps and from
# desired source area
rasterizeSims <- function(
      simulations = NULL      # Simulated trajectories
    , raster      = NULL      # Raster onto which we rasterize
    , steps       = 500       # How many steps should be considered
    , area        = "Main"    # Simulations from which areas?
    , messages    = T         # Print update messages?
    , mc.cores    = detectCores() - 1
  ){

  # Subset to corresponding data
  sub <- simulations[which(
      simulations$StepNumber <= steps
    & simulations$Area %in% area
  ), ]

  # Make sure raster values are all 0
  values(raster) <- 0

  # Create spatial lines
  sub_traj <- sims2tracks(
      simulations = sub
    , id          = "TrackID"
    , messages    = messages
    , mc.cores    = mc.cores
  )

  # Remove data by converting into spatial lines
  sub_traj <- as(sub_traj, "SpatialLines")
  crs(sub_traj) <- CRS("+init=epsg:32734")

  # Rasterize lines onto the cropped raster
  if (messages){
    cat("Rasterizing spatial lines...\n")
  }
  heatmap <- rasterizeSpatstat(
      l        = sub_traj
    , r        = raster
    , mc.cores = 1
  )

  # Store heatmap to temporary file
  heatmap <- writeRaster(heatmap, tempfile())
  crs(heatmap) <- CRS("+init=epsg:32734")

  # Return the resulting heatmap
  return(heatmap)
}

################################################################################
#### Rasterize Trajectories
################################################################################
# Create a dataframe with all source points and points in time at which we want
# to rasterize trajectories
rasterized <- as_tibble(
  expand.grid(
      steps            = c(68, 125, 250, 500, 1000, 2000)
    , area             = unique(sims$Area)
    , stringsAsFactors = F
  )
)

# Add a column for temporary but unique filename. Make sure the tempdir has
# plenty of storage.
rasterized$filename <- tempfile(
    pattern = paste0(
        "steps_", rasterized$steps
      , "_area_", rasterized$area
      , "_"
    )
  , fileext = ".tif"
)

# Loop through the study design and reasterize trajectories
heatmaps <- list()
for (i in 1:nrow(rasterized)){

  # Create heatmap
  heatmaps[[i]] <- rasterizeSims(
      simulations = sims
    , raster      = r
    , steps       = rasterized$steps[i]
    , area        = rasterized$area[i]
    , messages    = T
    , mc.cores    = detectCores() - 1
  )

  # Clean garbage
  gc()

  # Print update
  cat(i, "/", nrow(rasterized), "done...\n")

}

# Combine maps
combined <- stack(heatmaps)

# Reproject them
combined <- rast(combined)
combined <- terra::project(combined, CRS("+init=epsg:4326"), method = "bilinear")
combined <- stack(combined)

# Crop them to the buffer
buffer <- readOGR("03_Data/03_Results/99_BufferArea.shp")
combined <- crop(combined, buffer)

# Store to file
writeRaster(combined, "03_Data/03_Results/99_Heatmaps.grd", overwrite = T)

# Add maps to the tibble
rasterized <- mutate(rasterized, heatmap = lapply(1:nlayers(combined), function(x){
  combined[[x]]
}))

# Store to file
write_rds(rasterized, "03_Data/03_Results/99_Heatmaps.rds")
rasterized <- read_rds("03_Data/03_Results/99_Heatmaps.rds")

library(rasterVis)
plot(rasterized$heatmap[[1]], col = cols)
library(RColorBrewer)
cols <- rev(colorRampPalette(brewer.pal(11, 'Spectral'))(256))
myTheme <- rasterTheme(region = cols)
levelplot(rasterized$heatmap[[1]], margin = F, par.settings = myTheme, cuts = 20)
levelplot(rasterized$heatmap[[1]])

tm_shape(rasterized$heatmap[[1]]) + tm_raster(palette = "-Spectral", style = "cont")
p1 <- grid.grab()
tm_shape(rasterized$heatmap[[2]]) + tm_raster(palette = "-Spectral", style = "cont")
p2 <- grid.grab()
tm_shape(rasterized$heatmap[[6]] + rasterized$heatmap[[12]]) + tm_raster(palette = "-Spectral", style = "cont")
p3 <- grid.grab()
library(cowplot)
plot_grid(p1, p2, p3, )
test
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))
print(p1, vp = viewport(layout.pos.col = 1, layout.pos.row = 1))
print(p2, vp = viewport(layout.pos.col = 2, layout.pos.row = 1))
print(p3, vp = viewport(layout.pos.col = 2, layout.pos.row = 2))

p1 <- levelplot(combined[[1]], margin = F)
p2 <- levelplot(combined[[2]], margin = F)
p3 <- levelplot(combined[[3]], margin = F)
p4 <- levelplot(combined[[4]], margin = F)
p5 <- levelplot(combined[[5]], margin = F)
p6 <- levelplot(combined[[6]], margin = F)
p7 <- levelplot(combined[[6]] + combined[[12]], margin = F)

library(gridExtra)
lattice.options(
  layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0))
)
lattice.options(
  layout.heights=list(bottom.padding=list(x=-1), top.padding=list(x=-1)),
  layout.widths=list(left.padding=list(x=0), right.padding=list(x=0))
)
p1 <- levelplot(combined[[1]]
  , margin = F
  , xlab = NULL
  , ylab = NULL
  , colorkey = F
  , scales = list(
      x = list(draw = FALSE)
    , y = list(draw = FALSE)
  )
)
p1 <- levelplot(combined[[1]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p2 <- levelplot(combined[[2]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p3 <- levelplot(combined[[3]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p4 <- levelplot(combined[[4]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p5 <- levelplot(combined[[5]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p6 <- levelplot(combined[[6]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p7 <- levelplot(combined[[7]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p8 <- levelplot(combined[[8]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p9 <- levelplot(combined[[9]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p10 <- levelplot(combined[[10]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p11 <- levelplot(combined[[11]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p12 <- levelplot(combined[[12]], margin = F, xlab = NULL, ylab = NULL, colorkey = F, scales = list(x = list(draw = FALSE), y = list(draw = FALSE)), par.settings = myTheme, cuts = 99)
p13 <- levelplot(combined[[6]] + combined[[12]], margin = F, xlab = NULL, ylab = NULL, par.settings = myTheme, cuts = 99)
lay <- rbind(
    c(1, 2, 3, 4, 5, 6)
  , c(13, 13, 13, 13, 13, 13)
  , c(13, 13, 13, 13, 13, 13)
  , c(13, 13, 13, 13, 13, 13)
  , c(13, 13, 13, 13, 13, 13)
  , c(7, 8, 9, 10, 11, 12)
)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, layout_matrix = lay)
levelplot(combined[[6]] + combined[[12]], par.settings = myTheme, cuts = 90, margin = F)
plot(combined[[6]] + combined[[12]], col = rev(cols(100)), horizontal = T, box = F, axes = F)
p13 <- tm_shape(combined[[6]] + combined[[12]]) + tm_raster(palette = "-Spectral", style = "cont")

library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(n = 11, "Spectral"))

################################################################################
#### Visualizations
################################################################################
# Required Data
# rasterized <- read_rds("03_Data/03_Results/99_RasterizedSimulations.rds")
# heatmaps <- stack("03_Data/03_Results/99_RasterizedSimulations.tif")
points1 <- shapefile("03_Data/03_Results/99_SourcePoints.shp")
points2 <- shapefile("03_Data/03_Results/99_SourcePoints2.shp")
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>%
  readOGR() %>%
  as("SpatialLines")
africa <- "03_Data/02_CleanData/00_General_Africa.shp" %>%
  readOGR() %>%
  as("SpatialLines")
africa_crop <- "03_Data/02_CleanData/00_General_Africa.shp" %>%
  readOGR() %>%
  crop(kaza)
nati        <- shapefile("03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(3Classes)")

# Subset to national parks
nati <- subset(nati, Desig == "National Park")

# Subset to national parks that we want to plot
nati <- subset(nati, Name %in% c("Mavinga", "Luengue-Luiana", "Kafue"
  , "Hwange", "Central Kalahari", "Chobe", "Moremi", "Matusadona", "Khaudum"))

# There is a double entry for Kafue, get rid of the erronous one
nati$Area <- gArea(nati, byid = TRUE)
nati <- subset(nati, Area != min(Area))

# Create a separate shapefile for the text. We have to change some of the
# coordinates to make sure that they don't overlap
nati_text <- nati
nati_text$x <- coordinates(nati_text)[, 1]
nati_text$y <- coordinates(nati_text)[, 2]
nati_text <- nati_text@data
nati_text$y[nati_text$Name == "Kafue"] <-
  nati_text$y[nati_text$Name == "Kafue"] + 0.5
nati_text$y[nati_text$Name == "Chobe"] <-
  nati_text$y[nati_text$Name == "Chobe"] - 0.1
nati_text$y[nati_text$Name == "Matusadona"] <-
  nati_text$y[nati_text$Name == "Matusadona"] - 0.1
coordinates(nati_text) <- c("x", "y")
crs(nati_text) <- CRS("+init=epsg:4326")

# Check how they align
plot(nati)
points(nati_text)

# Add "NP" to the text (on a new line)
head(nati_text)
nati_text$Name <- paste0(nati_text$Name, "\nNP")

# We only keep the countries of interest in the cropped africa file
africa_crop <- subset(africa_crop, COUNTRY %in% c(
    "Angola"
  , "Namibia"
  , "Botswana"
  , "Zimbabwe"
  , "Zambia")
)

# Normalize heatmaps
for (i in 1:nlayers(heatmaps)){
  heatmaps[[i]] <- normalizeMap(heatmaps[[i]])
}

# Prepare plot of each map
p <- list()
for (i in 1:nlayers(heatmaps)){
  p[[i]] <- tm_shape(heatmaps[[i]]) +
      tm_raster(
          palette         = "-Spectral"
        , style           = "cont"
        , title           = "Traversal Frequency"
        , labels          = c("Low-Frequency", "", "High-Frequency")
        , legend.reverse  = T
      ) +
    tm_grid(
        n.x                 = 5
      , n.y                 = 5
      , labels.inside.frame = FALSE
      , lines               = FALSE
      , ticks               = TRUE
    ) +
    tm_shape(nati) +
      tm_borders(
          col   = "black"
        , alpha = 0.6
      ) +
    tm_shape(kaza) +
      tm_lines(
          col = "black"
        , lwd = 2
      ) +
    tm_shape(nati_text) +
      tm_text("Name"
        , col       = "black"
        , alpha     = 0.6
        , fontface  = 3
        , size      = 0.5
        , shadow    = F
      ) +
    tm_shape(africa) +
      tm_lines(
          col = "black"
        , lwd = 1
        , lty = 2
      ) +
    tm_shape(africa_crop) +
      tm_text("COUNTRY"
        , col       = "black"
        , just      = "bottom"
        , fontface  = 2
      ) +
    tm_layout(
        legend.text.color   = "white"
      , legend.title.color  = "white"
    ) +
    tm_scale_bar(
        position    = "left"
      , text.size   = 0.5
      , text.color  = "white"
      , width       = 0.125
    ) +
    tm_compass(
        text.color   = "white"
      , color.dark   = "white"
      , color.light  = "white"
    ) +
    tm_credits(paste0("Points: ", rasterized$sampling[[i]], "\nSteps: ", rasterized$steps[[i]])
    , position  = c("right", "top")
    , size      = 1.5
    , col       = "white"
  )
}

# Store the plots
for (i in 1:length(p)){
  name <- paste0(
      "04_Manuscript/99_RasterizedSims_Points"
    , rasterized$sampling[[i]]
    , "_Steps"
    , rasterized$steps[[i]]
    , ".png"
  )
  tmap_save(tm = p[[i]], name)
}

# ################################################################################
# #### Combine Heatmaps
# ################################################################################
# # Extend all heatmaps to that they match the extent of the reference raster
# rasterized <- mutate(rasterized, heatmaps_full = map(heatmaps, function(x){
#
#   # Match extent with reference raster
#   x <- terra::expand(x, r, value = 0)
#
#   # Reclassify NaNs to 0s (I don't know really why "value = 0" doesn't do the
#   # job)
#   x <- classify(x, rcl = data.frame(from = NaN, to = 0))
#
#   # Store to temporary file on the external drive
#   x <- writeRaster(x, tempfile(
#       fileext = ".tif"
#     , tmpdir  = "/media/david/SharedSpace/RStuff/Temporary"
#   )
#
#   # Collect garbage
#   gc()
#
#   # Return the expanded raster
#   return(x)
# }))
#
# # Combine all heatmaps belonging to the same number of steps and the same
# # sampling type
# rasterized <- rasterized %>%
#
#   # Make sure we group by steps and sampling type
#   group_by(steps, sampling) %>%
#
#   # Nest data and create tibble
#   nest() %>%
#
#   # Create a new column that contains the summed raster
#   mutate(heatmap = map(data, function(x){
#
#     # Sum the respective rasters
#     summed <- sum(rast(x$heatmaps_full))
#
#     # Create write to file
#     summed <- writeRaster(summed, tempfile(fileext = ".tif"))
#     return(summed)
#   }))
#
# # Put the summed heatmaps into a single stack and assign sensible names
# static <- subset(rasterized, sampling == "Static")
# random <- subset(rasterized, sampling == "Random")
# summed_static <- rast(static$heatmap)
# summed_random <- rast(random$heatmap)
#
# # Assign sensible layer names
# names(summed_static) <- paste0("steps_", static$steps)
# names(summed_random) <- paste0("steps_", random$steps)
#
# # Write stack to file
# writeRaster(summed_static
#   , "03_Data/03_Results/99_Simulations_Static.tif"
#   , overwrite = T
# )
# writeRaster(summed_random
#   , "03_Data/03_Results/99_Simulations.Random.tif"
#   , overwrite = T
# )
#
# # Visualize
# plot(log(summed_static + 1), col = magma(50))
# plot(log(summed_random + 1), col = magma(50))
#
# ################################################################################
# #### Compare Heatmaps
# ################################################################################
# # To compare heatmaps we need to normalize them so that they cover the same
# # range
# normalizeMap <- function(raster){
#   xmin <- minmax(raster)[1, ]
#   xmax <- minmax(raster)[2, ]
#   raster_norm <- (raster - xmin)/(xmax - xmin)
#   return(raster_norm)
# }
# normalizeMap(summed_static)
#
# # Compare heatmaps based on correlation
# pairs(summed_static)
# pairs(summed_random)





# # Load protected areas
# prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
#
# # Create SpatialPoints for first location of each trajectory
# first <- sims %>%
#   dplyr::select("x", "y", "TrackID") %>%
#   group_by(TrackID) %>%
#   slice(1) %>%
#   SpatialPointsDataFrame(
#       coords      = cbind(.[["x"]], .[["y"]])
#     , proj4string = CRS("+init=epsg:4326")
#   )
#
# # Assess from which protected area each trajectory left
# first$Origin <- as.character(over(first, prot)$Name)
# first <- first@data[, c("TrackID", "Origin")]
#
# # Join information to simulated tracks
# sims <- left_join(sims, first, by = "TrackID")
# sims$Origin[sims$Area == "Buffer"] <- "Buffer"
#
# # Make sure there are no NAs
# sum(is.na(sims$Origin))
