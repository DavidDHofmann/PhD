################################################################################
#### Number of Connections Between Protected Areas
################################################################################
# Description: In this script, we analyze the simulated dispersal trajectories
# and identify the number of connections between protected areas. We do so for
# different points in time as well as for different point sampling regimes.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(pbmcapply)      # To run on multiple cores with progress bar
library(rgeos)          # For spatial ananylsis
library(igraph)         # For network analysis

# Load custom functions
source("Functions.r")

################################################################################
#### Identify Reached Areas
################################################################################
################################################################################
#### NOTE: THIS CAN BE DONE USING IGRAPH RATHER THAN THE STUPID RGEOS
################################################################################

# Load source areas and points
areas <- shapefile("03_Data/03_Results/99_SourceAreas2.shp")
points <- shapefile("03_Data/03_Results/99_SourcePoints2.shp")

# Load the simulated dispersal trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")

# Function to identify all areas reached by a trajectory
getReached <- function(trajectories = NULL, areas = NULL){

    # Identify reached areas
    reached <- trajectories %>%
      gIntersects(areas, ., byid = T) %>%
      as.data.frame() %>%
      setNames(areas$ID) %>%
      cbind(trajectories@data, .) %>%
      dplyr::select(ID, StartPoint, everything()) %>%
      gather(key = EndPoint, value = Reached, 3:ncol(.)) %>%
      subset(Reached) %>%
      dplyr::select(-Reached) %>%
      mutate(EndPoint = as.numeric(EndPoint)) %>%
      arrange(ID, StartPoint, EndPoint)

  # Return reached areas
  return(reached)
}

# Prepare a tibble indicating the points in time and point sampling regimes for
# which we want to identify the areas reached.
reached <- as_tibble(
  expand.grid(
      steps     = c(68, 125, 250, 500, 1000, 2000)
    , sampling  = c("Static", "Random")
    , bootstrap = 1:2
  )
)

# Identify reached areas for all of the above study designs
reached$AreasReached <- pbmclapply(1:nrow(reached)
  , mc.cores            = 1
  , ignore.interactive  = T
  , function(i){

  # Sample 50 trajectories per start point
  sims_sub <- sims %>%
    subset(.
      , StepNumber <= reached$steps[i]
      & PointSampling == reached$sampling[i]
    ) %>%
    group_by(StartPoint, ID) %>%
    nest() %>%
    group_by(StartPoint) %>%
    sample_n(50) %>%
    unnest(data)

  # Prepare trajectories
  trajs <- sims2tracks(
      simulations = sims_sub
    , steps       = reached$steps[i]
    , sampling    = reached$sampling[i]
    , mc.cores    = 2
  )

  # Identify areas reached
  result <- getReached(trajectories = trajs, areas = areas)

  # Collect garbage
  gc()

  # Return it
  return(tibble(result))
})

# Identify the number of times an area was reached from each source point
reached$FreqAreasReached <- lapply(1:nrow(reached), function(x){
  reached$AreasReached[[x]] %>%
    group_by(StartPoint, EndPoint) %>%
    summarize(FreqAreasReached = n())
})

# Identify diversity of areas reached
reached$UniqueAreasReached <- lapply(1:nrow(reached), function(i){
  reached$FreqAreasReached[[i]] %>%
    group_by(StartPoint) %>%
    summarize(UniqueAreasReached = n())
})

# Let's check the results
print(reached)
reached$AreasReached[[1]]
reached$FreqAreasReached[[1]]
reached$UniqueAreasReached[[1]]

# Store the results to file
write_rds(reached, "03_Data/03_Results/99_AreasReached.rds")

################################################################################
#### Create Network
################################################################################
# Reload results
reached <- read_rds("03_Data/03_Results/99_AreasReached.rds")

# Create networks for all study facettes
reached$Networks <- lapply(1:nrow(reached), function(x){

  # Coerce dataframe to graph
  net <- graph_from_data_frame(reached$FreqAreasReached[[x]], vertices = points$ID)

  # Assign edge weights
  E(net)$weight <- reached$FreqAreasReached[[x]]$FreqAreasReached

  # Remove self-loops
  net <- simplify(net)

  # Return result
  return(net)
})

# Function to plot a network
plotNet <- function(net){

  # Identify isolate nodes
  isolates <- degree(net) == 0

  # Prepare color palette
  pal1 <- colorRampPalette(c("orange", "brown"))

  # We can actually plot the network graph on top of other spatial plots
  plot(areas, col = "gray20", bg = "black")
  plot(areas[areas$ID %in% names(isolates[isolates]), ], col = "gray60", add = T)
  plot(net
    , add                 = T
    , layout              = coordinates(points)
    , rescale             = F
    , vertex.size         = 8
    , vertex.label        = NA
    , vertex.color        = pal1(10)[1]
    , vertex.frame.color  = pal1(10)[1]
    , edge.color          = colTrans(pal1(10)[5], percent = 50)
    , edge.arrow.size     = 0
    , edge.curved         = 0.3
    , edge.width          = sqrt(E(net)$weight) / 3
  )
}

# Plot some of the networks
plotNet(reached$Networks[[1]])
plotNet(reached$Networks[[2]])
plotNet(reached$Networks[[3]])
plotNet(reached$Networks[[4]])
plotNet(reached$Networks[[5]])
plotNet(reached$Networks[[6]])

################################################################################
#### Calculate Network Metrics
################################################################################
# Calculate all centralization metrics
metrics <- lapply(1:nrow(reached), function(x){
  degree        <- centralization.degree(reached$Networks[[x]])$centralization
  betweenness   <- centralization.betweenness(reached$Networks[[x]])$centralization
  density       <- edge_density(reached$Networks[[x]])
  diameter      <- diameter(reached$Networks[[x]])
  transitivity  <- transitivity(reached$Networks[[x]])
  return(data.frame(cbind(degree, betweenness, density, diameter, transitivity)))
}) %>% do.call(rbind, .)
reached$Degree        <- metrics$degree
reached$Betweenness   <- metrics$betweenness
reached$Density       <- metrics$density
reached$Diameter      <- metrics$diameter
reached$Transitivity  <- metrics$transitivity

# Subset to a single network
net <- reached$Networks[[6]]

# Plot some node-level metrics
plot(areas, col = "gray20", bg = "black")
plot(net
  , add                 = T
  , layout              = coordinates(points)
  , rescale             = F
  , vertex.size         = degree(net)
  , vertex.label        = NA
  , vertex.color        = pal1(10)[1]
  , vertex.frame.color  = pal1(10)[1]
  , edge.color          = colTrans(pal1(10)[5], percent = 50)
  , edge.arrow.size     = 0
  , edge.curved         = 0.3
  , edge.width          = sqrt(E(net)$weight) / 3
)

# Plot some node-level metrics
plot(areas, col = "gray20", bg = "black")
plot(net
  , add                 = T
  , layout              = coordinates(points)
  , rescale             = F
  , vertex.size         = betweenness(net) / 10
  , vertex.label        = NA
  , vertex.color        = pal1(10)[1]
  , vertex.frame.color  = pal1(10)[1]
  , edge.color          = colTrans(pal1(10)[5], percent = 50)
  , edge.arrow.size     = 0
  , edge.curved         = 0.3
  , edge.width          = sqrt(E(net)$weight) / 3
)
