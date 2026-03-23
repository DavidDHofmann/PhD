################################################################################
#### Rasterization of Simulated Dispersal Trajectories
################################################################################
# Description: In this script, we rasterize the simulated dispersal trajectories
# and create "heatmaps".

# Clear R's brain
rm(list = ls())

# Load required packages
library(terra)        # For quick raster manipulation
library(tidyverse)    # For data wrangling
library(igraph)       # To compute betweenness
library(pbmcapply)    # For multicore use with progress bar

# Change the working directory
wd <- ifelse(Sys.info()["sysname"] == "Linux"
  , "/media/david/SharedSpace/SwitchDrive/02_Academia/02_PhD/Chapter_3"
  , "D:/SwitchDrive/02_Academia/02_PhD/Chapter_3"
)
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

# Load some background maps
source <- vect("03_Data/02_CleanData/Sources.gpkg")
water  <- vect("03_Data/02_CleanData/MajorWaters.gpkg")

# Load reference raster and aggregate to a coarser resolution
r <- rast("03_Data/02_CleanData/Raster.tif")
r <- aggregate(r, fact = 8, fun = mean)
res(r) * 111

################################################################################
#### Functions
################################################################################
# Interpolate simulations if desired
interpolateSims <- function(sim, eps = 1, mc.cores = 1) {

  # Setup apply function to be used (depending on multicore use or not)
  if (mc.cores > 1) {
      applier <- function(...) {
        pbmclapply(ignore.interactive = T, mc.cores = mc.cores, ...)
      }
    } else {
      applier <- function(...) {
        lapply(...)
      }
  }

  # Run interpolation
  sim_nested <- nest(sim, Data = -id)
  sim_nested$Data <- applier(sim_nested$Data, function(x) {
    interpolateTrack(x, by = eps)
  })
  int <- unnest(sim_nested, Data)
  return(int)
}

# Reload simulations design (drop unnecessary columns)
sims <- "03_Data/03_Results/Simulations.rds" %>%
  read_rds() %>%
  dplyr::select(-c(x, y, Timestamp, Done))

# Function to compute connectivity maps and ipc
connectivity <- function(sim, r, areas = NULL, singlecount = T, metrics = c("Heatmap", "Betweenness", "Interpatch")) {

  # Check validity
  metrics <- match.arg(metrics
    , choices    = c("Heatmap", "Betweenness", "Interpatch")
    , several.ok = T
  )

  # Prepare reference raster
  r <- r[[1]]
  r[] <- 0

  # Prepare result list
  results <- vector("list", length = length(metrics))
  names(results) <- metrics

  # Generate rasterized areas
  areas$ID <- 1:length(areas)
  areas_r  <- terra::rasterize(areas, r, field = "ID")

  # Generate visitation histories
  sim$cell <- cellFromXY(r, cbind(sim$x, sim$y))
  sim$area <- areas_r[][sim$cell]

  # Create heatmap
  if ("Heatmap" %in% metrics) {
    heat_count <- sim %>%
      dplyr::select(id, cell) %>%
      distinct() %>%
      count(cell)
    heat <- r
    heat[heat_count$cell] <- heat_count$n
    names(heat) <- "Heatmap"
    results[["Heatmap"]] <- heat
  }

  # Compute betweenness if desired
  if ("Betweenness" %in% metrics) {

    # Identify cell transitions for betweenness
    cell_transitions <- sim %>%
      dplyr::select(-area) %>%
      subset(!is.na(cell)) %>%
      group_by(id) %>%
      mutate(
          from = cell
        , to   = lead(cell)
      ) %>%
      group_by(id, from, to) %>%
      subset(from != to) %>%
      na.omit() %>%
      summarize(TotalConnections = n(), .groups = "drop")

    # If each trajectory should only count once, adjust the number of connections
    if (singlecount) {
      cell_transitions$TotalConnections <- 1
    }

    # Aggregate visitation histories across all paths
    cell_transitions <- cell_transitions %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections), .groups = "drop") %>%
      mutate(weight = mean(TotalConnections) / TotalConnections)

    # Compute betweenness
    ver  <- seq_len(ncell(r))
    net  <- graph_from_data_frame(cell_transitions, vertices = as.numeric(ver))
    betw <- setValues(r, betweenness(net))
    names(betw) <- "Betweenness"
    results[["Betweenness"]] <- betw
  }

  # Compute interpatch connectivity
  if ("Interpatch" %in% metrics) {
    ipc <- sim %>%
      subset(!is.na(area) & !is.nan(area)) %>%
      group_by(id) %>%
      mutate(
          source  = first(area)
        , duration = difftime(timestamp, min(timestamp), units = "hours")
      ) %>%
      ungroup() %>%
      subset(area != source)

    # If there are no rows left, we have to return NA
    if (nrow(ipc) == 0) {
        ipc <- NA
      } else {
        ipc <- ipc %>%
          group_by(id, source, area) %>%
          summarize(duration = min(duration), .groups = "drop") %>%
          group_by(source, area) %>%
          summarize(duration = mean(duration), success = length(id), .groups = "drop") # Could add duration_sd here
    }
    results[["Interpatch"]] <- ipc

  }

  # Return metrics
  return(results)
}

################################################################################
#### Compute Connectivity
################################################################################
# Reload simulations design (drop unnecessary columns)
sims <- "03_Data/03_Results/Simulations.rds" %>%
  read_rds() %>%
  dplyr::select(-c(x, y, Timestamp)) %>%
  subset(file.exists(Filename))

# Define design through which to loop and specify filename for connectivity
# results
design <- sims %>%
  dplyr::select(Source, Formula, ModelSeasons, FittingCovariates, PredictionCovariates, ModelCode) %>%
  distinct() %>%
  mutate(Filename = paste0(
      "03_Data/03_Results/Connectivity/"
    , ModelCode, "_"
    , Source
    , ".rds")
  ) %>%
  mutate(Done = file.exists(Filename))

# Write the design to file
write_rds(design, "03_Data/03_Results/Connectivity.rds")

# Create folder
if (!dir.exists("03_Data/03_Results/Connectivity/")) {
  dir.create("03_Data/03_Results/Connectivity/")
}

# Loop through the design and generate connectivity estimates for design
# combinations that haven't been computed yet
pb <- txtProgressBar(min = 0, max = nrow(design))
for (i in seq_len(nrow(design))) {

  # # Testing
  # i <- 1

  # Skip if files already exist
  if (file.exists(design$Filename[i])) {
    setTxtProgressBar(pb, i)
    next
  }

  # Subset to relevatn simulations and load them
  sims_sub <- sims %>%
    subset(Source == design$Source[i] & ModelCode == design$ModelCode[i] & Replicate <= 1e6) %>%
    mutate(Data = map(Filename, read_rds)) %>%
    dplyr::select(Source, Replicate, Data) %>%
    unnest(Data) %>%
    rename(id = Replicate, timestamp = Timestamp)

  # Interpolate to xx meters
  cat("Interpolating simulations\n")
  sims_sub_inter <- interpolateSims(sims_sub, eps = 1000 / 111000, mc.cores = detectCores() - 1)

  # Compute connectivity
  cat("Computing connectivity\n")
  conn <- connectivity(
      sim         = sims_sub_inter
    , r           = r
    , areas       = source
    , singlecount = T
    # , metrics     = "Heatmap"
  )

  # Make sure rasters are loaded
  conn$Heatmap     <- wrap(conn$Heatmap)
  conn$Betweenness <- wrap(conn$Betweenness)

  # Store it to file
  setTxtProgressBar(pb, i)
  write_rds(conn, design$Filename[i])
}

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/02_Analysis/04_Connectivity.rds")
cat("Done :)\n")
