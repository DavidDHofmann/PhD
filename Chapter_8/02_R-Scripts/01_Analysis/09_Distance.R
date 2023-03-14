################################################################################
#### Distance to Humans
################################################################################
# Description: Computing the distance to humans for each simulated step

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(terra)      # To handle spatial data
library(raster)     # To handle spatial data
library(tidyverse)  # To wrangle data
library(maptools)   # To compute point density
library(spatstat)   # To compute point density
library(pbmcapply) # For multicore progress bar

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load distance to humans layer
humans <- rast("03_Data/02_CleanData/DistanceToHumans.tif")

# Load dispersal simulations
sims <- read_rds("03_Data/03_Results/DispersalSimulation.rds")

# Keep only desired columns
sims <- sims[, c("x", "y", "sl_", "ta_", "TrackID", "StepNumber", "SourceArea", "FloodLevel")]

# For later, extract the unique combinations of Track_IDs and FloodLevels
allid <- sims %>% select(TrackID, FloodLevel, SourceArea) %>% distinct()

################################################################################
#### Distance to Humans
################################################################################
# Create directories into which we will store the output rasters
dir.create("03_Data/03_Results/99_Distances", showWarnings = F)

# Calculate the distance to humans for each step
sims$DistanceToHumans <- terra::extract(humans, cbind(sims$x, sims$y))[, 1]

# We only care about locations that are closer than 500 meters to humans (this
# will strictly reduce the number of datapoints we have to keep track off, but
# also implies that we have to keep track of the original TrackIDs (allid) when
# we bootstrap.
sims <- subset(sims, DistanceToHumans < 500)

# Design through which to loop
design <- expand_grid(
    FloodLevel = unique(allid$FloodLevel)
  , SourceArea = unique(allid$SourceArea)
)
design$Filename <- with(design, paste0("03_Data/03_Results/99_Distances/Distance_SourceArea", SourceArea, "_FloodLevel", FloodLevel, ".tif"))

# Loop throgh the design and generate density maps of the locations that are
# within a threshold distance to the nearest "human-influenced" pixel
if (!file.exists("03_Data/03_Results/Distance.tif")) {
    distmaps <- list()
    for (i in 1:nrow(design)) {
      cat("Computing distance map", i, "out of", nrow(design), "\n")
      if (file.exists(design$Filename[i])) {
        distmaps[[i]] <- rast(design$Filename[[i]])
      } else {
        closeto <- subset(sims,
            FloodLevel == design$FloodLevel[i] &
            SourceArea == design$SourceArea[i]
        )
        closeto <- vect(x = cbind(closeto$x, closeto$y)
          , crs  = "epsg:4326"
          , atts = closeto
        )
        distmaps[[i]] <- terra::rasterize(closeto
          , humans
          , fun        = function(i) { length(i) }
          , background = 0
        )
        writeRaster(distmaps[[i]], design$Filename[i], overwrite = T)
      }
    }

    # Combine maps, reproject them, crop, and store them
    combined <- do.call(c, combined)
    writeRaster(combined, "03_Data/03_Results/Distance.tif", overwrite = T)
  } else {
    combined <- stack("03_Data/03_Results/Distance.tif")
}

# Add maps to the tibble
design <- mutate(design, Distance = lapply(1:nlayers(combined), function(x) {
  combined[[x]]
}))

# Create "Global Metrics"
maps1 <- design
maps1$Level = "Local"
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
  x <- raster(x)
  return(x)
})

# Store the tibble
write_rds(hwcm, "03_Data/03_Results/Distance.rds")

################################################################################
#### Quantify Differences with Bootstrapping
################################################################################
# To quantify differences in different focal regions, let's load the shapefile
# of the areas of interest
aoi <- vect("03_Data/02_CleanData/AreasOfInterest.shp")

# Compute the area of the different areas
aoi$Area <- expanse(aoi, unit = "km")

# Create a reference map that has a resolution of 1 km2
ref <- rast(ext(humans), resolution = 1 / 111)

# We will need to subsample the simulated trajectories. Thus, nest them
# accordingly.
sims_nested <- nest(sims, Data = -c(FloodLevel, TrackID))

# Size of the samples to take
n_sample <- allid %>%
  count(FloodLevel) %>%
  pull(n) %>%
  first()

# Design through which we will loop
design <- expand_grid(
    Replicate  = 1:1000
  , FloodLevel = unique(sims_nested$FloodLevel)
)

# Repeatedly rasterize coordinates within a certain distance to humans, then
# extract the point density below our areas of interest
design$HWC <- pbmclapply(1:nrow(design), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {

  # Sample IDs (with replacement)
  keep <- allid %>%
    subset(FloodLevel == design$FloodLevel[x]) %>%
    slice_sample(replace = T, n = n_sample)

  # Only keep simulations that fall into those IDs
  closeto <- sims_nested %>%
    subset(TrackID %in% keep$TrackID) %>%
    unnest(Data)

  # Count the number of locations in the different areas of interest
  n_points <- terra::extract(aoi, cbind(closeto$x, closeto$y))
  n_points <- count(n_points, Name)

  # Remove those not falling into an aoi
  n_points <- na.omit(n_points)

  # Calculate density
  n_points <- left_join(n_points, values(aoi), by = "Name")
  n_points <- mutate(n_points, Density = n / Area)

  # Return those
  return(n_points)
})

# Compute summary statistics
metrics <- design %>%
  unnest(HWC) %>%
  group_by(FloodLevel, Name) %>%
  summarize(Number = mean(Density), SE = sd(Density)) %>%
  ungroup() %>%
  arrange(Name, FloodLevel) %>%
  group_by(Name) %>%
  mutate(NumberRev = rev(Number)) %>%
  mutate(Percent = round(Number / NumberRev * 100 - 100, 0)) %>%
  mutate(FloodLevel = factor(FloodLevel, levels = c("Min", "Max"))) %>%
  select(-NumberRev) %>%
  mutate(combined = paste0(format(round(Number), big.mark = "'"), "$ \\pm $", sprintf("%.2f", round(SE, 2))))

# Store them to file
write_rds(metrics, "03_Data/03_Results/DistanceAOI.rds")

# Store some specific metrics to file as well
metrics %>% subset(Name == "Maun" & FloodLevel == "Max") %>% pull(Percent) %>%
  as.character() %>%
  writeLines("04_Manuscript/99_DensityHWCMaunPercentage.tex")
metrics %>% subset(Name == "Panhandle" & FloodLevel == "Max") %>% pull(Percent) %>%
  as.character() %>%
  writeLines("04_Manuscript/99_DensityHWCMaunPercentage.tex")

# ggplot(metrics, aes(x = FloodLevel, y = Number, ymin = Number - SE, ymax = Number + SE, col = FloodLevel, fill = FloodLevel)) +
#   geom_col(alpha = 0.5) +
#   geom_errorbar(width = 0.2) +
#   facet_grid(~ Name) +
#   theme_minimal() +
#   scale_color_manual(values = c("orange", "cornflowerblue")) +
#   scale_fill_manual(values = c("orange", "cornflowerblue"))
