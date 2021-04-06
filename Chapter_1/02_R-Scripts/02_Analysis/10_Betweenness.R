################################################################################
#### Analysis of Simulated Dispersal Events
################################################################################
# Description: In this script we will analyse the simulated trajectories using a
# network approach. Ultimately, we will generate a map illustrating the
# betweenness of each pixel in our study area.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(parallel)       # To run on multiple cores
library(pbmcapply)      # To run on multiple cores with progress bar
library(igraph)         # For network analysis
library(viridis)        # For nice colors
library(RColorBrewer)   # For more nice colors
library(davidoff)       # Custom functions

# Set a seed
set.seed(12345)

# Suppress scientific notation
options(scipen = 999)

################################################################################
#### Functions to Calculate Network Metrics
################################################################################
# We will create and use networks in order to derive network metrics that we can
# plot spatially. Let's prepare a function that we can use to derive all desired
# network metrics and depict them on a raster.
netMet <- function(
    network   = NULL    # The network based on which the metrics are calculated
  , raster    = NULL    # The raster onto which the metrics are calculated
  , tempfile  = F       # Should the resulting raster go to a temporary file?
  , metrics = c("betweenness", "closeness", "degree")
  ){

    # Calculate the desired network metrics and return them as rasters
    result <- vector(mode = "list", length = 3)
    if ("betweenness" %in% metrics){
      betweenness <- raster
      values(betweenness) <- betweenness(network)
      names(betweenness) <- "betweenness"
      result[[1]] <- betweenness
    }
    if ("closeness" %in% metrics){
      closeness <- raster
      values(closeness) <- closeness(network)
      names(closeness) <- "closeness"
      result[[2]] <- closeness
    }
    if ("degree" %in% metrics){
      degree <- raster
      values(degree) <- degree(network)
      names(degree) <- "degree"
      result[[3]] <- degree
    }

    # Remove NULLs from the list
    result <- plyr::compact(result)

    # Put all into a stack
    result <- stack(result)

    # In case the user desires to store a temporary file, do so
    if (tempfile){
      result <- writeRaster(result, tempfile())
    }

    # Return the final stack
    return(result)
}

# Function to combine metrics from multiple graphs
stackMet <- function(
    metrics   = NULL                  # The stack of metrics
  , metric    = "betweenness"         # The metric to be combined
  , tempfile  = F
  , fun       = function(x){mean(x)}  # The function to apply to stack the metric
  ){
    stacked <- stack(lapply(metrics, function(z){z[[metric]]}))
    funned  <- calc(stacked, fun)
    if (tempfile){
      funned <- writeRaster(funned, tempfile())
    }
    return(funned)
}

# We now want to create a visitation history for each trajectory. This history
# depicts all transitions from one raster cell to another. To do this
# repeatedly, we write a function to retrieve the visitation history from a
# sequence of values
visitHist <- function(x, singlecount = F){
  transitions <- data.frame(from = lag(x), to = x) %>%
    group_by(from, to) %>%
    na.omit() %>%
    summarize(TotalConnections = n(), .groups = "drop")
  if (singlecount){
    transitions$TotalConnections = 1
  }
  return(transitions)
}

################################################################################
#### Load and Clean Data
################################################################################
# Load the simulated dispersal trajectories
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")

# Prepare design through which we want to loop
design <- expand_grid(
    steps     = c(125, 500, 2000)
  # , bootstrap = c(1:100)
)

# # Let's shuffle the matrix (makes prediction of calculation time in pbmclapply
# # more accurate)
# design <- design[sample(nrow(design), replace = F), ]

# Load the reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Need to coarsen resolution
# r  <- aggregate(r, fact = 5000 / 250, fun = max)
r  <- aggregate(r, fact = 2500 / 250, fun = max)

# Note that the resolution of the raster determines the number of vertices in
# the resulting network (each raster cell will be turned into a vertex). Let's
# therefore keep track of the IDs of each raster cell.
vertices <- 1:ncell(r)

# The resolution of the raster therefore also determines the layout of the
# network (i.e. the spatial location of each vertex). Let's keep track of the
# coordinates of each raster cell (i.e. each vertex) so that we can nicely plot
# the network graph afterwards.
lay  <- as.matrix(as.data.frame(r, xy = T)[, c(1, 2)])

# Fill the rasters with unique cell values (we will use the values as cell IDs)
values(r)  <- vertices

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# At each coordinate of the simulated trajectories we now extract the cell IDs
# from the different rasters
visits <- data.frame(
    TrackID    = sims$TrackID
  , StepNumber = sims$StepNumber
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , r          = raster::extract(r, sims)
)

# Remove simulations
rm(sims)
gc()

# Loop through the design and calculate betweenness map
maps <- list()
for (i in 1:nrow(design)){

  # Subset data to desired steps
  sub <- visits[visits$StepNumber <= design$steps[i], ]

  # Nest tracks
  sub <- nest(sub, data = -TrackID)

  # Create visitation history
  cat("Getting visitation history...\n")
  history <- pbmclapply(
      X                  = sub$data
    , ignore.interactive = T
    , mc.cores           = 1
    , FUN                = function(y){
        visitHist(y$r, singlecount = T)
    }) %>%
    do.call(rbind, .) %>%
    group_by(from, to) %>%
    summarize(TotalConnections = sum(TotalConnections), .groups = "drop") %>%
    ungroup() %>%
    mutate(weight = mean(TotalConnections) / TotalConnections)

  # Create network
  cat("Creating graph...\n")
  net <- graph_from_data_frame(history, vertices = vertices)

  # Calculate Betweenness
  cat("Calculating betweenness...\n")
  met <- netMet(
      network  = net
    , raster   = r
    , metrics  = "betweenness"
    , tempfile = T
  )

  # Put into list
  maps[[i]] <- met

  # Update
  cat(i, "out of", nrow(design), "done\n")

}

# Store them
maps <- stack(maps)
plot(maps)
plot(sqrt(maps[[3]]), col = magma(100))
plot(sqrt(focal(maps[[3]], w = matrix(1, 3, 3), fun = sum)), col = magma(100))
plot(sqrt(focal(maps[[3]], w = matrix(1, 3, 3), fun = max)), col = magma(100))
writeRaster(maps, "03_Data/03_Results/99_Betweenness2500.grd")

# At each coordinate of the simulated trajectories we now extract the cell IDs
# from the different rasters
visits <- data.frame(
    TrackID = sims$TrackID
  , x       = coordinates(sims)[, 1]
  , y       = coordinates(sims)[, 2]
  , R10000  = raster::extract(r10000, sims)
  , R5000   = raster::extract(r5000, sims)
  , R2500   = raster::extract(r2500, sims)
)
# Let's check what the function does exactly
visitHist(c(1, 2, 2, 3, 4, 5, 1, 2), singlecount = F)
visitHist(c(1, 2, 2, 3, 4, 5, 1, 2), singlecount = T)
visitHist(c(1, 2, 2, NA, 4, 5, 1, 2), singlecount = T)

# We want to retrieve the visitation history to each trajectory seperately, so
# let's nest them.
visits <- visits %>% nest(data = -TrackID)

# Apply the function to each trajectory for each raster. We therefore identify
# the visitation history of each trajectory for each of the different spatial
# resolutions.
visits$History10000 <- pbmclapply(1:nrow(visits)
  , mc.cores            = detectCores() / 2
  , ignore.interactive  = T
  , function(x){
    visitHist(visits$data[[x]]$R10000, singlecount = T)
  })
visits$History5000 <- pbmclapply(1:nrow(visits)
  , mc.cores            = detectCores() / 2
  , ignore.interactive  = T
  , function(x){
    visitHist(visits$data[[x]]$R5000, singlecount = T)
  })
visits$History2500 <- pbmclapply(1:nrow(visits)
  , mc.cores            = detectCores() / 2
  , ignore.interactive  = T
  , function(x){
    visitHist(visits$data[[x]]$R2500, singlecount = T)
  })

# Check the result
print(visits)
print(visits$History10000[[1]])
print(visits$History5000[[1]])
print(visits$History2500[[1]])

# Important to note: we now derived the vistation history for each simulated
# trajectory individually. We could easily put them together and calculate a
# single visitation history throughout all dispersal trajectories. We will
# actually do this now.
visits10000 <- do.call(rbind, visits$History10000) %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections)) %>%
      ungroup() %>%
      mutate(weight = mean(TotalConnections) / TotalConnections)
visits5000 <- do.call(rbind, visits$History5000) %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections)) %>%
      ungroup() %>%
      mutate(weight = mean(TotalConnections) / TotalConnections)
visits2500 <- do.call(rbind, visits$History2500) %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections)) %>%
      ungroup() %>%
      mutate(weight = mean(TotalConnections) / TotalConnections)

# Create graphs
net10000 <- graph_from_data_frame(visits10000, vertices = vertices10000)
net5000 <- graph_from_data_frame(visits5000, vertices = vertices5000)
net2500 <- graph_from_data_frame(visits2500, vertices = vertices2500)

# Check if weighted
is_weighted(net10000)
is_weighted(net5000)
is_weighted(net2500)

# Check if connected
is_connected(net10000)
is_connected(net5000)
is_connected(net2500)

# Calculate metrics
met10000 <- netMet(network = net10000, raster = r10000, metrics = "betweenness")
met5000  <- netMet(network = net5000, raster = r5000, metrics = "betweenness")
met2500  <- netMet(network = net2500, raster = r2500, metrics = "betweenness")


################################################################################
#### Load and Clean Data
################################################################################
# Load the simulated dispersal trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")

# # Subset to simulations of interest
# sims <- subset(sims, StepNumber <= 2000)
#
# For now, only keep a few trajectories
# sims <- sims %>%
#   nest(data = -TrackID) %>%
#   slice_sample(n = 10) %>%
#   unnest(data)

# Load the reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Prepare multiple rasters with different resolutions. Use the reference raster
# for this
r10000  <- aggregate(r, fact = 10000 / 250)
r5000   <- aggregate(r, fact = 5000 / 250)
r2500   <- aggregate(r, fact = 2500 / 250)

# Store them to file to save some memory
r10000  <- writeRaster(r10000, tempfile())
r5000   <- writeRaster(r5000, tempfile())
r2500   <- writeRaster(r2500, tempfile())

# Clean remaining garbage
gc()

# Note that the resolution of the raster determines the number of vertices in
# the resulting network (each raster cell will be turned into a vertex). Let's
# therefore keep track of the IDs of each raster cell.
vertices10000 <- 1:ncell(r10000)
vertices5000  <- 1:ncell(r5000)
vertices2500  <- 1:ncell(r2500)

# The resolution of the raster therefore also determines the layout of the
# network (i.e. the spatial location of each vertex). Let's keep track of the
# coordinates of each raster cell (i.e. each vertex) so that we can nicely plot
# the network graph afterwards.
lay10000  <- as.matrix(as.data.frame(r10000, xy = T)[, c(1, 2)])
lay5000   <- as.matrix(as.data.frame(r5000, xy = T)[, c(1, 2)])
lay2500   <- as.matrix(as.data.frame(r2500, xy = T)[, c(1, 2)])

# Fill the rasters with unique cell values (we will use the values as cell IDs)
values(r10000)  <- vertices10000
values(r5000)   <- vertices5000
values(r2500)   <- vertices2500

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# At each coordinate of the simulated trajectories we now extract the cell IDs
# from the different rasters
visits <- data.frame(
    TrackID = sims$TrackID
  , x       = coordinates(sims)[, 1]
  , y       = coordinates(sims)[, 2]
  , R10000  = raster::extract(r10000, sims)
  , R5000   = raster::extract(r5000, sims)
  , R2500   = raster::extract(r2500, sims)
)

# We now want to create a visitation history for each trajectory. This history
# depicts all transitions from one raster cell to another. To do this
# repeatedly, we write a function to retrieve the visitation history from a
# sequence of values
visitHist <- function(x, singlecount = F){
  transitions <- data.frame(from = lag(x), to = x) %>%
    group_by(from, to) %>%
    na.omit() %>%
    summarize(TotalConnections = n(), .groups = "drop")
  if (singlecount){
    transitions$TotalConnections = 1
  }
  return(transitions)
}

# Let's check what the function does exactly
visitHist(c(1, 2, 2, 3, 4, 5, 1, 2), singlecount = F)
visitHist(c(1, 2, 2, 3, 4, 5, 1, 2), singlecount = T)
visitHist(c(1, 2, 2, NA, 4, 5, 1, 2), singlecount = T)

# We want to retrieve the visitation history to each trajectory seperately, so
# let's nest them.
visits <- visits %>% nest(data = -TrackID)

# Apply the function to each trajectory for each raster. We therefore identify
# the visitation history of each trajectory for each of the different spatial
# resolutions.
visits$History10000 <- pbmclapply(1:nrow(visits)
  , mc.cores            = detectCores() / 2
  , ignore.interactive  = T
  , function(x){
    visitHist(visits$data[[x]]$R10000, singlecount = T)
  })
visits$History5000 <- pbmclapply(1:nrow(visits)
  , mc.cores            = detectCores() / 2
  , ignore.interactive  = T
  , function(x){
    visitHist(visits$data[[x]]$R5000, singlecount = T)
  })
visits$History2500 <- pbmclapply(1:nrow(visits)
  , mc.cores            = detectCores() / 2
  , ignore.interactive  = T
  , function(x){
    visitHist(visits$data[[x]]$R2500, singlecount = T)
  })

# Check the result
print(visits)
print(visits$History10000[[1]])
print(visits$History5000[[1]])
print(visits$History2500[[1]])

# Important to note: we now derived the vistation history for each simulated
# trajectory individually. We could easily put them together and calculate a
# single visitation history throughout all dispersal trajectories. We will
# actually do this now.
# visits10000 <- do.call(rbind, visits$History10000) %>%
#       group_by(from, to) %>%
#       summarize(TotalConnections = sum(TotalConnections)) %>%
#       ungroup() %>%
#       mutate(weight = mean(TotalConnections) / TotalConnections)
visits5000 <- do.call(rbind, visits$History5000) %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections)) %>%
      ungroup() %>%
      mutate(weight = mean(TotalConnections) / TotalConnections)
# visits2500 <- do.call(rbind, visits$History2500) %>%
#       group_by(from, to) %>%
#       summarize(TotalConnections = sum(TotalConnections)) %>%
#       ungroup() %>%
#       mutate(weight = mean(TotalConnections) / TotalConnections)

# Create graphs
# net10000 <- graph_from_data_frame(visits10000, vertices = vertices10000)
net5000 <- graph_from_data_frame(visits5000, vertices = vertices5000)
# net2500 <- graph_from_data_frame(visits2500, vertices = vertices2500)

# # Check if weighted
# is_weighted(net10000)
# is_weighted(net5000)
# is_weighted(net2500)
#
# # Check if connected
# is_connected(net10000)
# is_connected(net5000)
# is_connected(net2500)

# Calculate metrics
met10000 <- netMet(network = net10000, raster = r10000, metrics = "betweenness")
met5000  <- netMet(network = net5000, raster = r5000, metrics = "betweenness")
met2500  <- netMet(network = net2500, raster = r2500, metrics = "betweenness")

################################################################################
#### TESTING: SINGLE TRAJECTORIES
################################################################################
# Select an index
i <- 1001

# Extract first trajectory
traj <- visits$data[[i]]
coordinates(traj) <- c("x", "y")
traj <- spLines(traj)

# Visualize it
plot(r10000)
plot(traj, add = T, col = "red")

# Create graph from visitation history
graph <- graph_from_data_frame(visits$History10000[[i]], vertices = vertices10000)

# Visualize the transitions
plot(simplify(graph, remove.loops = T), layout = lay10000
  , vertex.size        = 0
  , vertex.color       = colTrans("gray")
  , vertex.frame.color = colTrans("gray")
  , vertex.label       = NA
  , edge.size          = 0.1
  , edge.arrow.size    = 0
  , edge.curved        = 0
  , edge.color         = "orange"
)

# Calculate network metrics
mets <- netMet(network = graph, raster = r10000)

# Visualize them
plot(mets)
plot(mets[[1]])

################################################################################
#### TESTING: WEIGHTS VS NO WEIGHTS
################################################################################
# Prepare data
visits_all_ww1 <- tibble(
    History10000  = list(
      do.call(rbind, visits$History10000) %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections)) %>%
      ungroup() %>%
      mutate(weight = max(TotalConnections) - TotalConnections + 1)
  )
)
visits_all_ww2 <- tibble(
    History10000  = list(
      do.call(rbind, visits$History10000) %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections)) %>%
      ungroup() %>%
      mutate(weight = mean(TotalConnections) / TotalConnections)
  )
)
visits_all_nw <- tibble(
    History10000  = list(
      do.call(rbind, visits$History10000) %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections)) %>%
      ungroup()
  )
)

# Create graphs
net_ww1 <- graph_from_data_frame(visits_all_ww1$History10000[[1]]
  , vertices = vertices10000)
net_ww2 <- graph_from_data_frame(visits_all_ww2$History10000[[1]]
  , vertices = vertices10000)
net_nw <- graph_from_data_frame(visits_all_nw$History10000[[1]]
  , vertices = vertices10000)

# Check if weighted
is_weighted(net_ww1)
is_weighted(net_ww2)
is_weighted(net_nw)

# Check if connected
is_connected(net_ww1)
is_connected(net_ww2)
is_connected(net_nw)

# Calculate metrics
res_ww1 <- netMet(network = net_ww1, raster = r10000)
res_ww2 <- netMet(network = net_ww2, raster = r10000)
res_nw  <- netMet(network = net_nw, raster = r10000)

# # Calculate centralization betweenness
# centralization.betweenness(net_ww1)$centralization
# centralization.betweenness(net_ww2)$centralization
# centralization.betweenness(net_nw)$centralization

# Compare results
par(mfrow = c(2, 2))
plot(sqrt(res_ww1[["betweenness"]]), col = viridis(50), main = "With Weights fun 1")
plot(sqrt(res_ww2[["betweenness"]]), col = viridis(50), main = "With Weights fun 2")
plot(sqrt(res_nw[["betweenness"]]), col = viridis(50), main = "Without Weights")

par(mfrow = c(2, 2))
plot(sqrt(res_ww1[["degree"]]), col = viridis(20), main = "With Weights fun 1")
plot(sqrt(res_ww2[["degree"]]), col = viridis(20), main = "With Weights fun 2")
plot(sqrt(res_nw[["degree"]]), col = viridis(20), main = "Without Weights")

# Plot betweenness against degree
par(mfrow = c(1, 2))
plot(sqrt(res_ww2[["betweenness"]]), col = viridis(20), main = "With Weights fun 2")
plot(sqrt(res_ww2[["degree"]]), col = viridis(20), main = "With Weights fun 2")
