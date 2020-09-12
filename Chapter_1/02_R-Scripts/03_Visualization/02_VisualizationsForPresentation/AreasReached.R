################################################################################
#### Plot Reached Areas
################################################################################
# Description: Plot the areas reached from a source point

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/media/david/My Passport/Backups/WildDogs/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(pbmcapply)      # To run on multiple cores with progress bar
library(rgeos)          # For spatial ananylsis
library(igraph)         # For network analysis
library(davidoff)       # Custom functions

################################################################################
#### Prepare Data
################################################################################
# Load source areas and points
areas <- shapefile("03_Data/03_Results/99_SourceAreas2.shp")
points <- shapefile("03_Data/03_Results/99_SourcePoints2.shp")

# Visualize
plot(areas)
text(x = points, labels = points$ID, cex = 1)

# Load the simulated dispersal trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Subset to a single start point
sims <- subset(sims, StartPoint == 38)

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
  )
)

# Identify reached areas for all of the above study designs
reached$AreasReached <- pbmclapply(1:nrow(reached)
  , mc.cores            = 1
  , ignore.interactive  = T
  , function(i){

  # Sample 50 trajectories per start point
  sims_sub <- subset(sims
      , StepNumber <= reached$steps[i]
      & PointSampling == reached$sampling[i]
  )

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
}) %>% .[["value"]]

# Identify the number of times an area was reached from each source point
reached$FreqAreasReached <- lapply(1:nrow(reached), function(x){
  reached$AreasReached[[x]] %>%
    group_by(StartPoint, EndPoint) %>%
    summarize(FreqAreasReached = n(), .groups = "drop")
})

# Identify diversity of areas reached
reached$UniqueAreasReached <- lapply(1:nrow(reached), function(i){
  reached$FreqAreasReached[[i]] %>%
    group_by(StartPoint) %>%
    summarize(UniqueAreasReached = n(), .groups = "drop")
})

# Let's check the results
print(reached)
reached$AreasReached[[1]]
reached$FreqAreasReached[[1]]
reached$UniqueAreasReached[[1]]

################################################################################
#### Create Network
################################################################################
# Create networks for all study facettes
reached$Networks <- lapply(1:nrow(reached), function(x){

  # Coerce dataframe to graph
  net <- graph_from_data_frame(reached$AreasReached[[x]][, c(2, 3)], vertices = points$ID)

  # Assign edge weights
  E(net)$weight <- reached$FreqAreasReached[[x]]$FreqAreasReached

  # Remove self-loops
  net <- simplify(net)

  # Return result
  return(net)
})

# Function to plot a network
plotNet <- function(net, main = NULL){

  # Identify isolate nodes
  isolates <- degree(net) == 0

  # Prepare color palette
  pal1 <- colorRampPalette(c("orange", "brown"))

  # We can actually plot the network graph on top of other spatial plots
  plot(areas, col = "gray20", bg = "black", main = main, lwd = 4)
  plot(areas[areas$ID %in% names(isolates[isolates]), ], col = "gray60", add = T, lwd = 4)
  plot(net
    , add                 = T
    , layout              = coordinates(points)
    , rescale             = F
    , vertex.size         = 14
    , vertex.label        = NA
    , vertex.color        = pal1(10)[1]
    , vertex.frame.color  = pal1(10)[1]
    , edge.color          = colTrans(pal1(10)[5], percent = 30)
    , edge.arrow.size     = 0
    , edge.curved         = 0.3
    , edge.width          = sqrt(E(net)$weight)
  )
}

# Plot some of the networks
for (i in 1:nrow(reached)){
  name <- paste0("test", i, ".png")
  png(name, width = 1980, height = 1980, pointsize = 50)
  plotNet(
      reached$Networks[[i]]
    , main = paste0(
       "Steps = "
      , reached$steps[[i]]
      , " , Sampling = "
      , reached$sampling[[i]]
    )
  )
  dev.off()
}
