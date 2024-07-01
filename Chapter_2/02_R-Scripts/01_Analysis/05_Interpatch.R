################################################################################
#### Interpatch Connectivity
################################################################################
# Description: Computing interpatch connectivity from simulated dispersal paths

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

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

# Specify if you'd like durations to be returned in days or steps
days <- T

# Load reference raster
r <- rast("03_Data/02_CleanData/ReferenceRaster.tif")

# Load source areas
area <- vect("03_Data/02_CleanData/SourceAreas.shp")

# Rasterize them to the reference raster
area_r <- terra::rasterize(area, y = r, field = "ID")

# Prepare a tibble containing the different metrics we'd like to compute
metrics <- expand_grid(
    Metric = c("Dispersal Success", "Dispersal Duration", "Egression", "Immigration")
  , Level  = c("Overall", "BySource")
  , Data   = list(NA)
)

################################################################################
#### Identifying Visits
################################################################################
# Before we can compute metrics, we need to identify through which area each of
# the simulated trajectories moved. Specifically, we are interested in
# identifying through which source areas they moved, as well as into which
# egression zones they moved.

# Load dispersal simulations
sims <- read_rds("03_Data/03_Results/DispersalSimulation.rds")

# Keep only desired columns
sims <- sims[, c("x", "y", "TrackID", "StepNumber", "SourceArea", "FloodLevel")]

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# Identify through which designated areas the dispersers moved
visits <- data.frame(
    TrackID    = sims$TrackID
  , StepNumber = sims$StepNumber
  , FloodLevel = sims$FloodLevel
  , SourceArea = sims$SourceArea
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , Area       = raster::extract(raster(area_r), sims)
)

# Ignore any step in an "na" or "nan" area
visits <- subset(visits, !is.na(Area) & !is.nan(Area))

# Calculate, for each step, the distance to the first coordinate
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

# Identify how long (in steps or days) it takes to reach the different areas
visits <- visits %>%
  rename(SourceArea = SourceArea, CurrentArea = Area) %>%
  group_by(TrackID, FloodLevel, SourceArea, CurrentArea) %>%
  summarize(
      StepNumber        = min(StepNumber) * ifelse(days, 1 / 5, 1)
    , DistanceFromFirst = min(DistanceFromFirst)
    , .groups           = "drop"
  ) %>%
  arrange(TrackID, StepNumber)

# Bootstrap inter-patch connectivity metrics of interest
visits_nested <- nest(visits, Data = -c(TrackID, FloodLevel, SourceArea))
n_sample <- 1000
bootstrapped <- pbmclapply(1:1000, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {

  # Resample data and count connections between areas and the steps needed to
  # make those
  dat_sampled <- visits_nested %>%
    group_by(FloodLevel, SourceArea) %>%
    slice_sample(replace = T, n = n_sample) %>%
    ungroup() %>%
    mutate(TrackID = 1:n()) %>%
    unnest(Data) %>%
    group_by(TrackID, FloodLevel, SourceArea, CurrentArea) %>%
    summarize(DispersalDuration = min(StepNumber), .groups = "drop") %>%
    arrange(TrackID, DispersalDuration)

  # Compute dispersal from source to source
  ipc_source_source_dispersal <- dat_sampled %>%
    subset(SourceArea != CurrentArea & CurrentArea <= 6) %>%
    group_by(FloodLevel, SourceArea, CurrentArea) %>%
    summarize(
        DispersalDuration = mean(DispersalDuration)
      , DispersalSuccess  = length(unique(TrackID))
      , .groups           = "drop"
    ) %>%
    mutate(Type = "Dispersal")

  # Compute egression from source to buffer (egression) zone
  ipc_source_buffer_egression <- dat_sampled %>%
    subset(SourceArea != CurrentArea & CurrentArea > 6) %>%
    group_by(TrackID) %>%
    slice_head(n = 1) %>% # This ensures that only the first area through which the individual leaves the study area is considered
    group_by(FloodLevel, SourceArea, CurrentArea) %>%
    summarize(
        DispersalDuration = mean(DispersalDuration)
      , DispersalSuccess  = length(unique(TrackID))
      , .groups           = "drop"
    ) %>%
    mutate(Type = "Egression")

  # Compute dispersal from from source to overall
  ipc_source_overall_dispersal <- dat_sampled %>%
    subset(SourceArea != CurrentArea & CurrentArea <= 6) %>%
    group_by(FloodLevel, SourceArea) %>%
    summarize(
        DispersalDuration = mean(DispersalDuration)
      , DispersalSuccess  = length(unique(TrackID))
      , .groups           = "drop"
    ) %>%
    mutate(CurrentArea = "Overall", Type = "Dispersal")

  # Compute egression from source to overall
  ipc_source_overall_egression <- dat_sampled %>%
    subset(SourceArea != CurrentArea & CurrentArea > 6) %>%
    group_by(TrackID) %>%
    slice_head(n = 1) %>% # This ensures that only the first area through which the individual leaves the study area is considered
    group_by(FloodLevel, SourceArea) %>%
    summarize(
        DispersalDuration = mean(DispersalDuration)
      , DispersalSuccess  = length(unique(TrackID))
      , .groups           = "drop"
    ) %>%
    mutate(CurrentArea = "Overall", Type = "Egression")

  # Compute dispersal from overall to source
  ipc_overall_source_dispersal <- dat_sampled %>%
    subset(SourceArea != CurrentArea & CurrentArea <= 6) %>%
    group_by(FloodLevel, CurrentArea) %>%
    summarize(
        DispersalDuration = mean(DispersalDuration)
      , DispersalSuccess  = length(unique(TrackID))
      , .groups           = "drop"
    ) %>%
    mutate(SourceArea = "Overall", Type = "Dispersal")

  # Compute connections into egression
  ipc_overall_source_egression <- dat_sampled %>%
    subset(SourceArea != CurrentArea & CurrentArea > 6) %>%
    group_by(TrackID) %>%
    slice_head(n = 1) %>% # This ensures that only the first area through which the individual leaves the study area is considered
    group_by(FloodLevel, CurrentArea) %>%
    summarize(
        DispersalDuration = mean(DispersalDuration)
      , DispersalSuccess  = length(unique(TrackID))
      , .groups           = "drop"
    ) %>%
    mutate(SourceArea = "Overall", Type = "Egression")

  # Compute overall number of individuals moving from one area to another
  ipc_overall_overall_dispersal <- dat_sampled %>%
    subset(SourceArea != CurrentArea & CurrentArea <= 6) %>%
    group_by(FloodLevel, TrackID) %>%
    summarize(DispersalDuration = min(DispersalDuration), .groups = "drop") %>%
    group_by(FloodLevel) %>%
    summarize(
        DispersalDuration = mean(DispersalDuration)
      , DispersalSuccess  = length(unique(TrackID))
      , .groups           = "drop"
    ) %>%
    mutate(SourceArea = "Overall", CurrentArea = "Overall", Type = "Dispersal")

  # Compute overall egression
  ipc_overall_overall_egression <- dat_sampled %>%
    subset(SourceArea != CurrentArea & CurrentArea > 6) %>%
    group_by(TrackID) %>%
    slice_head(n = 1) %>% # This ensures that only the first area through which the individual leaves the study area is considered
    group_by(FloodLevel) %>%
    summarize(
        DispersalDuration = mean(DispersalDuration)
      , DispersalSuccess  = length(unique(TrackID))
      , .groups           = "drop"
    ) %>%
    mutate(SourceArea = "Overall", CurrentArea = "Overall", Type = "Egression")

  # Put them together and assign bootstrap indicator
  boot <- rbind(
      ipc_source_source_dispersal
    , ipc_source_buffer_egression
    , ipc_source_overall_dispersal
    , ipc_source_overall_egression
    , ipc_overall_source_dispersal
    , ipc_overall_source_egression
    , ipc_overall_overall_dispersal
    , ipc_overall_overall_egression
  )
  boot$Bootstrap <- x

  # Return all
  return(boot)
})

# Compute summaries across replicates
ipc_bootstrapped <- bootstrapped %>%
  do.call(rbind, .) %>%
  group_by(FloodLevel, SourceArea, CurrentArea, Type) %>%
  summarize(
      DispersalSuccessSD  = sd(DispersalSuccess)
    , DispersalSuccess    = mean(DispersalSuccess)
    , DispersalDurationSD = sd(DispersalDuration)
    , DispersalDuration   = mean(DispersalDuration)
    , .groups             = "drop"
  )

# Self-loops, i.e. movements from a source area to itself are currently absent
# from the data, as we did not consider those. In addition, we may also be
# lacking connections between some areas. Hence, let's create a "full" dataset
# for reference. We will fill the missing values with NAs.
full_dispersal <- expand_grid(
    SourceArea  = c(as.character(1:6), "Overall")
  , CurrentArea = as.character(1:6)
  , FloodLevel  = c("Min", "Max")
  ) %>%
  mutate(Type = "Dispersal")
full_egression <- expand_grid(
    SourceArea  = c(as.character(1:6), "Overall")
  , CurrentArea = c(as.character(7:14), "Overall")
  , FloodLevel  = c("Min", "Max")
  ) %>%
  mutate(Type = "Egression")
full <- rbind(full_dispersal, full_egression)
ipc_bootstrapped <- full_join(full, ipc_bootstrapped
  , by = c("SourceArea", "CurrentArea", "FloodLevel", "Type")
)

# Anything which is not a self-loop, but contains NA is virtually no movement
# between the respective areas.
indices <- with(ipc_bootstrapped, SourceArea != CurrentArea & is.na(DispersalSuccess))
ipc_bootstrapped[indices, c("DispersalSuccessSD", "DispersalSuccess")] <- 0
ipc_bootstrapped[indices, c("DispersalDurationSD", "DispersalDuration")] <- Inf

# Store visits to file and reload them
write_rds(ipc_bootstrapped, "03_Data/03_Results/InterpatchConnectivityBootstrapped.rds")
ipc_bootstrapped <- read_rds("03_Data/03_Results/InterpatchConnectivityBootstrapped.rds")

################################################################################
#### Additional Metrics: Number of Individuals Reaching any other Area
################################################################################
# Get the number of individuals that reach any other area
metrics <- list()
metrics[[1]] <- ipc_bootstrapped %>%
  subset(SourceArea == "Overall" & CurrentArea == "Overall" & Type == "Dispersal") %>%
  dplyr::select(FloodLevel, Metric = DispersalSuccess, MetricSE = DispersalSuccessSD) %>%
  mutate(MetricName = "NumberReachingOtherSourceAreas")

################################################################################
#### Additional Metrics: Dispersal Duration before Reaching any other Area
################################################################################
# Identify number of steps required before reaching another area
metrics[[2]] <- ipc_bootstrapped %>%
  subset(SourceArea == "Overall" & CurrentArea == "Overall" & Type == "Dispersal") %>%
  dplyr::select(FloodLevel, Metric = DispersalDuration, MetricSE = DispersalDurationSD) %>%
  mutate(MetricName = "DurationReachingOtherSourceAreas")

################################################################################
#### Additional Metrics: Number of Individuals Moving into area 6
################################################################################
# Get the number of individuals that reach area 6
metrics[[3]] <- ipc_bootstrapped %>%
  subset(SourceArea == "Overall" & CurrentArea == "6" & Type == "Dispersal") %>%
  dplyr::select(FloodLevel, Metric = DispersalSuccess, MetricSE = DispersalSuccessSD) %>%
  mutate(MetricName = "NumberReaching6")

# Identify number of steps required before reaching area 6
metrics[[4]] <- ipc_bootstrapped %>%
  subset(SourceArea == "Overall" & CurrentArea == "6" & Type == "Dispersal") %>%
  dplyr::select(FloodLevel, Metric = DispersalDuration, MetricSE = DispersalDurationSD) %>%
  mutate(MetricName = "DurationReaching6")

################################################################################
#### Additional Metrics: Metrics for improved Connectivity
################################################################################
# For a few selected source areas we get an increase in dispersal success
better <- ipc_bootstrapped %>%
  subset(SourceArea != "Overall" & CurrentArea %in% c("1", "2", "3", "4", "5", "6") & Type == "Dispersal") %>%
  select(SourceArea, CurrentArea, FloodLevel, DispersalSuccess) %>%
  pivot_wider(names_from = FloodLevel, values_from = DispersalSuccess) %>%
  mutate(Percent = Max / Min * 100 - 100) %>%
  mutate(Better = Percent > 0) %>%
  subset(Better) %>%
  arrange(desc(Percent))

# Let's create a combined dataframe (including the SEs)
metrics[[5]] <- ipc_bootstrapped %>%
  select(FloodLevel, SourceArea, CurrentArea, DispersalSuccess, DispersalSuccessSD) %>%
  left_join(., better, by = c("SourceArea", "CurrentArea")) %>%
  subset(!is.na(Better)) %>%
  subset(Percent == max(Percent)) %>%
  mutate(Percent = round(Percent)) %>%
  dplyr::select(FloodLevel, Metric = DispersalSuccess, MetricSE = DispersalSuccessSD) %>%
  mutate(MetricName = "Connectivity5to4")

################################################################################
#### Additional Metrics: Temporary Egression
################################################################################
# Compute temporary egression
metrics[[6]] <- ipc_bootstrapped %>%
  subset(Type == "Egression" & SourceArea == "Overall" & CurrentArea == "Overall") %>%
  dplyr::select(FloodLevel, Metric = DispersalSuccess, MetricSE = DispersalSuccessSD) %>%
  mutate(MetricName = "Egression")

################################################################################
#### Combine Metrics
################################################################################
# Compute percentages
metrics <- metrics %>%
  do.call(rbind, .) %>%
  mutate(combined = paste0(
      format(round(Metric), big.mark = ",")
    , " $\\pm$ "
    , round(MetricSE, 0)
    # , sprintf("%.2f", round(MetricSE, 2))
  )) %>%
  group_by(MetricName) %>%
  mutate(MetricRev = rev(Metric)) %>%
  mutate(Percent = round(Metric / MetricRev * 100 - 100, 0))
metrics

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/05_Interpatch.rds")

# Print to terminal
cat("Done :)\n")
