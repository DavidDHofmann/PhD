################################################################################
#### Analysis of Simulated Dispersal Events
################################################################################
# Description: In this script we will analyse the simulated trajectories using a
# network approach. Ultimately, we will generate a map illustrating the
# betweenness of each pixel in our study area.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/media/david/My Passport/Backups/WildDogs/15. PhD/00_WildDogs"
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

################################################################################
#### Load and Clean Data
################################################################################
# Load the simulated dispersal trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Subset to simulations of interest
sims <- subset(sims
  , StepNumber    <= 200
  & PointSampling == "Static"
)

# # For now, only keep 5 trajectories per source point
# set.seed(123)
# sims <- sims %>%
#   group_by(StartPoint, ID) %>%
#   nest() %>%
#   group_by(StartPoint) %>%
#   sample_n(5) %>%
#   unnest()

# Load the reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Prepare multiple rasters with different resolutions for scaling analysis. We
# will simply use our reference raster and coarsen its resolution.
r10000  <- aggregate(r, fact = 10000 / 250)
r5000   <- aggregate(r, fact = 5000 / 250)
r2500   <- aggregate(r, fact = 2500 / 250)

# Alternatively, we could also run the following commands
# r10000  <- raster(extent(r), resolution = metersToDegrees(10000))
# r5000   <- raster(extent(r), resolution = metersToDegrees(5000))
# r2500   <- raster(extent(r), resolution = metersToDegrees(2500))

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

# At each coordinate we now extract the cell IDs from the different rasters
visits <- data.frame(
    ID      = sims$ID
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

# We want to retrieve the visitation history to each trajectory seperately, so
# let's nest them.
visits <- visits %>% nest(data = -ID)

# Apply the function to each trajectory for each raster. We therefore identify
# the visitation history of each trajectory for each of the different spatial
# resolutions.
visits$History10000 <- pbmclapply(1:nrow(visits)
  , mc.cores            = detectCores() - 1
  , ignore.interactive  = T
  , function(x){
    visitHist(visits$data[[x]]$R10000, singlecount = T)
  })
visits$History5000 <- pbmclapply(1:nrow(visits)
  , mc.cores            = detectCores() - 1
  , ignore.interactive  = T
  , function(x){
    visitHist(visits$data[[x]]$R5000, singlecount = T)
  })
visits$History2500 <- pbmclapply(1:nrow(visits)
  , mc.cores            = detectCores() - 1
  , ignore.interactive  = T
  , function(x){
    visitHist(visits$data[[x]]$R2500, singlecount = T)
  })

# Check the result
print(visits)
print(visits$History10000[[1]])

# Important to note: we now derived the vistation history for each simulated
# trajectory individually. We could easily put them together and calculate a
# single visitation history throughout all dispersal trajectories.

################################################################################
#### Calculate Network Metrics
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

    # Remove NULLs from the lsit
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
################################################################################
#### TESTING: SINGLE TRJAJECTORIES
################################################################################
# Select an index
i <- 2000

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
plot(graph, layout = lay10000, vertex.size = 0, edge.size = 0.1, edge.arrow.size = 0, vertex.label = NA)

# Calculate network metrics
mets <- netMet(network = graph, raster = r10000)

# Visualize them
plot(mets)

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
      mutate(weight = max(TotalConnections) - TotalConnections) + 1
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
net_ww1 <- graph_from_data_frame(visits_all_ww1$History10000[[1]], vertices = vertices10000)
net_ww2 <- graph_from_data_frame(visits_all_ww2$History10000[[1]], vertices = vertices10000)
net_nw <- graph_from_data_frame(visits_all_nw$History10000[[1]], vertices = vertices10000)

# Check if weighted
is_weighted(net_ww1)
is_weighted(net_ww2)
is_weighted(net_nw)

# Calculate metrics
res_ww1 <- netMet(network = net_ww1, raster = r10000)
res_ww2 <- netMet(network = net_ww2, raster = r10000)
res_nw <- netMet(network = net_nw, raster = r10000)

# Calculate centralization betweenness
centralization.betweenness(net_ww1)$centralization
centralization.betweenness(net_ww2)$centralization
centralization.betweenness(net_nw)$centralization

# Compare results
getwd()
par(mfrow = c(2, 2))
plot(sqrt(res_ww1[["betweenness"]]), col = viridis(20), main = "With Weights fun 1")
plot(sqrt(res_ww2[["betweenness"]]), col = viridis(20), main = "With Weights fun 2")
plot(sqrt(res_nw[["betweenness"]]), col = viridis(20), main = "Without Weights")

par(mfrow = c(2, 2))
plot(sqrt(res_ww1[["degree"]]), col = viridis(20), main = "With Weights fun 1")
plot(sqrt(res_ww2[["degree"]]), col = viridis(20), main = "With Weights fun 2")
plot(sqrt(res_nw[["degree"]]), col = viridis(20), main = "Without Weights")

################################################################################
#### Approach I: Calculate Network Metrics Over All Trajectories
################################################################################
# Bind the visitation histories of all trajectories together and calculate edge
# weights. In igraph, edge weights are comparable to costs. That is, igraph uses
# the weight as a proxy for distance. For our purposes we therefore need to
# convert the "TotalConnections" metric into a "Distance" metric. One way to
# achieve this is to invert values using the formula f(x) = max(x) - x.
weights <- c(1, 2, 3, 4)
max(weights) - weights
1/weights
mean(weights)/weights

visits_all <- tibble(
    History10000  = list(
      do.call(rbind, visits$History10000) %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections)) %>%
      ungroup() %>%
      mutate(weight = max(TotalConnections) - TotalConnections) + 1
  )
  , History5000   = list(
      do.call(rbind, visits$History5000) %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections)) %>%
      ungroup() %>%
      mutate(weight = max(TotalConnections) - TotalConnections) + 1
  )
  , History2500   = list(
      do.call(rbind, visits$History2500) %>%
      group_by(from, to) %>%
      summarize(TotalConnections = sum(TotalConnections)) %>%
      ungroup() %>%
      mutate(weight = max(TotalConnections) - TotalConnections) + 1
  )
)

# We now use the visitation histories as edgelists to create networks / graphs
visits_all$Graph10000 <- list(
  graph_from_data_frame(visits_all$History10000[[1]], vertices = vertices10000)
)
visits_all$Graph5000 <- list(
  graph_from_data_frame(visits_all$History5000[[1]], vertices = vertices5000)
)
visits_all$Graph2500 <- list(
  graph_from_data_frame(visits_all$History2500[[1]], vertices = vertices2500)
)

# Make sure the networks are weighted
is_weighted(visits_all$Graph10000[[1]])
is_weighted(visits_all$Graph5000[[1]])
is_weighted(visits_all$Graph2500[[1]])

# Calculate network metrics
visits_all$Metrics10000 <- list(netMet(visits_all$Graph10000[[1]], r10000))
visits_all$Metrics5000 <- list(netMet(visits_all$Graph5000[[1]], r5000))
visits_all$Metrics2500 <- list(netMet(visits_all$Graph2500[[1]], r2500))

# Visualize them
plot(
    sqrt(visits_all$Metrics10000[[1]][["betweenness"]])
  , col   = viridis(20)
  , main  = "Betweenness"
)
plot(
    visits_all$Metrics10000[[1]][["closeness"]]
  , col   = viridis(20)
  , main  = "Closeness"
)
plot(
    visits_all$Metrics10000[[1]][["degree"]]
  , col   = viridis(20)
  , main  = "Degree"
)

###############################################################################
#### Approach II: Apply to each Trajectory Individually, Then Merge Metrics
################################################################################
# In contrast to be approach above, we know create a network for each dispersal
# trajectory individually.
# backup <- visits
# visits <- visits[sample(1:nrow(visits), 20), ]
visits$Graph10000 <- pbmclapply(1:nrow(visits)
  , mc.cores            = 1
  , ignore.interactive  = T
  , function(x){
    graph_from_data_frame(visits$History10000[[x]], vertices = vertices10000)
  })
visits$Graph5000 <- pbmclapply(1:nrow(visits)
  , mc.cores            = 1
  , ignore.interactive  = T
  , function(x){
    graph_from_data_frame(visits$History5000[[x]], vertices = vertices5000)
  })
visits$Graph2500 <- pbmclapply(1:nrow(visits)
  , mc.cores            = 1
  , ignore.interactive  = T
  , function(x){
    graph_from_data_frame(visits$History2500[[x]], vertices = vertices2500)
  })

# Calculate network metrics for each network
visits$Metrics10000 <- pbmclapply(1:nrow(visits)
  , mc.cores            = 1
  , ignore.interactive  = T
  , function(x){
      netMet(
          network   = visits$Graph10000[[x]]
        , raster    = r10000
        , metrics   = c("betweenness", "degree")
        , tempfile  = T
      )
  }
)
visits$Metrics5000 <- pbmclapply(1:nrow(visits)
  , mc.cores            = 1
  , ignore.interactive  = T
  , function(x){
      netMet(
          network   = visits$Graph5000[[x]]
        , raster    = r5000
        , metrics   = c("betweenness", "degree")
        , tempfile  = T
      )
  }
)
visits$Metrics2500 <- pbmclapply(1:nrow(visits)
  , mc.cores            = 1
  , ignore.interactive  = T
  , function(x){
      netMet(
          network   = visits$Graph2500[[x]]
        , raster    = r2500
        , metrics   = c("betweenness", "degree")
        , tempfile  = T
      )
  }
)

# Put metrics them into single maps
betweenness <- stackMet(
    metrics = visits$Metrics10000
  , metric  = "betweenness"
  , fun     = function(x){mean(x)}
)
degree <- stackMet(
    metrics = visits$Metrics10000
  , metric  = "degree"
  , fun     = function(x){mean(x)}
)

# Visualize
plot(betweenness, col = viridis(20))
plot(degree, col = viridis(20))
