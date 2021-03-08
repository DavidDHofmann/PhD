################################################################################
#### Rasterization of Simulated Dispersal Trajectories
################################################################################
# Description: In this script, we rasterize the simulated dispersal
# trajectories. To achieve this, we'll first rasterize trajectories by their
# source points. Afterwards, we can combine them them into a single "heatmap".

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(terra)            # For quick raster manipulation
library(raster)           # For general raster manipulation
library(tidyverse)        # For data wrangling
library(rgeos)            # For manipulating vector data
library(lubridate)        # For working with dates
library(viridis)          # For nicer colors
library(tictoc)           # To keep track of processing time
library(pbmcapply)        # To show progress bar in mclapply calls
library(tmap)             # For nice spatial plots
library(Cairo)            # To store plots
library(davidoff)         # Custom functions

################################################################################
#### Function to Rasterize Simulated Trajectory
################################################################################
# Function to rasterize trajectories of a desired source point after a desired
# amount of time (i.e. after a desired number of steps)
rasterizeSims <- function(
      simulations = NULL      # Simulated trajectories
    , steps       = 400       # How many steps should be considered
    , sampling    = "Static"  # For which sampling method we want to rasterize
  ){

  # Create tracks that fulfill the above requirements
  sub_traj <- sims2tracks(
      simulations = simulations
    , steps       = steps
    , sampling    = sampling
  )

  # Remove data by creating spatial lines
  sub_traj <- as(sub_traj, "SpatialLines")

  # Assign correct CRS
  crs(sub_traj) <- CRS("+init=epsg:4326")

  # Coerce to "vect" for faster rasterization
  sub_traj <- vect(sub_traj)

  # Speed up rasterization by cropping the raster
  r_crop <- crop(r, ext(sub_traj))

  # Rasterize lines onto the cropped raster
  heatmap <- rasterizeTerra(sub_traj, r_crop)

  # Return the resulting heatmap
  return(heatmap)
}

################################################################################
#### Load and Prepare Data
################################################################################
# Load the reference raster
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")

# For this part we can reduce the resolution of the raster drastically
r <- aggregate(r, fact = 10)

# Load the simulated dispersal trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Check out the number of rows
nrow(sims) / 1e6

################################################################################
#### Rasterize Trajectories Once
################################################################################
# Create a dataframe with all source points and points in time at which we want
# to rasterize trajectories
rasterized <- as_tibble(
  expand.grid(
      steps     = c(68, 125, 250, 500, 1000, 2000)
    , sampling  = c("Static", "Random")
  )
)

# Add a column for temporary but unique filename. Make sure the tempdir has
# plenty of storage.
rasterized$filename <- tempfile(
    pattern = paste0(
        "steps_", rasterized$steps
      , "_sampling_", rasterized$sampling
      , "_"
    )
  , fileext = ".tif"
)

# Rasterize simulated trajectories. Note that the number of cores you can use
# depends a bit on the amount of ram that is available. In some cases memory may
# overflow.
heatmaps <- pbmclapply(1:nrow(rasterized)
  , mc.cores            = 1
  , ignore.interactive  = T
  , FUN                 = function(i){

  # Rasterize trajectories
  heatmap <- rasterizeSims(
      simulations = sims
    , steps       = rasterized$steps[i]
    , sampling    = rasterized$sampling[i]
  )

  # Make sure the map is not stored in memory but on disk
  heatmap <- terra::writeRaster(heatmap, rasterized$filename[i], overwrite = T)

  # Clear cache
  gc()

  # Return the final raster
  return(heatmap)

})

# For whatever reason the final list does not correctly link to the rasters and
# we need to reload the rasters using their temporary filenames
heatmaps <- lapply(1:nrow(rasterized), function(x){
  rast(rasterized$filename[x])
})

# Prepare nice layernames
names <- paste0("Steps_", rasterized$steps, "_Sampling_", rasterized$sampling)

# Put heatmaps into a stack
heatmaps <- do.call(c, heatmaps)
names(heatmaps) <- names

# Store the final stack
heatmaps <- writeRaster(heatmaps
  , "03_Data/03_Results/99_RasterizedSimulations.tif"
  , overwrite = T
)

# Put them into the tibble
rasterized$heatmap <- vector(mode = "list", length = nrow(rasterized))
for (i in 1:nrow(rasterized)){
  rasterized$heatmap[[i]] <- heatmaps[[i]]
}

# Write the tibble to file too
write_rds(rasterized, "03_Data/03_Results/99_RasterizedSimulations.rds")

################################################################################
#### Rasterize Trajectories Repeatedly
################################################################################
# To get a sense of variability, we will need to create multiple heatmaps for
# each study design. Let's do so and create 50 heatmaps for each study "facette"
rasterized2 <- as_tibble(
  expand.grid(
      steps     = c(68, 125, 250, 500, 1000, 2000)
    , sampling  = c("Static", "Random")
    , bootstrap = 1:10
  )
)

# Add a column for temporary but unique filename. Make sure the tempdir has
# plenty of storage.
rasterized2$filename <- tempfile(
    pattern = paste0(
        "steps_", rasterized2$steps
      , "_sampling_", rasterized2$sampling
      , "_repetition_", rasterized2$bootstrap
      , "_"
    )
  , fileext = ".tif"
)

# Rasterize simulated trajectories. Note that the number of cores you can use
# depends a bit on the amount of ram that is available. In some cases memory may
# overflow.
heatmaps2 <- pbmclapply(1:nrow(rasterized2)
  , mc.cores            = 1
  , ignore.interactive  = T
  , FUN                 = function(i){

  # Only keep 50 simulations per source point and point sampling method
  sims_sub <- sims %>%
    subset(.
      , StepNumber <= rasterized2$steps[i]
      & PointSampling == rasterized2$sampling[i]
    ) %>%
    group_by(StartPoint, ID) %>%
    nest() %>%
    group_by(StartPoint) %>%
    sample_n(50) %>%
    unnest(data)

  # Rasterize trajectories
  heatmap <- rasterizeSims(
      simulations = sims_sub
    , steps       = rasterized2$steps[i]
    , sampling    = rasterized2$sampling[i]
  )

  # Make sure the map is not stored in memory but on disk
  heatmap <- terra::writeRaster(heatmap, rasterized2$filename[i], overwrite = T)

  # Clear cache
  gc()

  # Return the final raster
  return(heatmap)

})

# For whatever reason the final list does not correctly link to the rasters and
# we need to reload the rasters using their temporary filenames
heatmaps2 <- lapply(1:nrow(rasterized2), function(x){
  rast(rasterized2$filename[x])
})

# Prepare nice layernames
names <- paste0(
    "Steps_", rasterized2$steps
  , "_Sampling_", rasterized2$sampling
  , "_Bootstrap_", rasterized2$bootstrap
)

# Put heatmaps into a stack
heatmaps2 <- do.call(c, heatmaps)
names(heatmaps2) <- names

# Store the final stack
heatmaps2 <- writeRaster(heatmaps2
  , "03_Data/03_Results/99_RasterizedSimulationsBootstrap.tif"
  , overwrite = T
)

# Put them into the tibble
rasterized2$heatmap2 <- vector(mode = "list", length = nrow(rasterized2))
for (i in 1:nrow(rasterized2)){
  rasterized2$heatmap2[[i]] <- heatmaps2[[i]]
}

# Write the tibble to file too
write_rds(rasterized2, "03_Data/03_Results/99_RasterizedSimulationsBootstrap.rds")

################################################################################
#### Visualizations
################################################################################
# Required Data
rasterized <- read_rds("03_Data/03_Results/99_RasterizedSimulations.rds")
heatmaps <- stack("03_Data/03_Results/99_RasterizedSimulations.tif")
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
