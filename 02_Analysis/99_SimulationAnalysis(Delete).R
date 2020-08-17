################################################################################
#### Analysis of Simulated Dispersal Events
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(dismo)          # To calculate voronoi polygons
library(rgeos)          # For spatial analysis
library(viridis)        # For nice colors
library(parallel)       # To run on multiple cores
library(pbmcapply)      # To run on multiple cores with progress bar
library(tmap)           # For nice spatial plots
library(igraph)         # For network analysis
library(arrangements)   # To get permutations
library(tictoc)         # To keep track of time

# Load custom functions
source("Functions.r")

# Set a seed
set.seed(12345)

################################################################################
#### Prepare Source Areas and Source Points
################################################################################
# Load source areas and points
areas <- shapefile("03_Data/03_Results/99_SourceAreas2.shp")
points <- shapefile("03_Data/03_Results/99_SourcePoints2.shp")

# Load country borders
africa <- shapefile("03_Data/02_CleanData/00_General_Africa.shp")

# Rename column name (was abbreviated when stored)
names(points)[2] <- "ProtectedArea"

# Make sure that IDs of polygons and points match and visualize them
plot(areas, border = areas$ID)
plot(points, col = points$ID, add = T)
text(x = points, labels = points$ID, cex = 0.8)

# Rasterize source areas onto our reference raster.
r <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Rasterize source areas and make resulting raster factorial. Rasterized values
# will correspond to the area's ID.
areas_r <- rasterize(areas, r, field = "ID", background = NA)
areas_r <- ratify(areas_r)

# Visualize them
plot(areas_r, col = rainbow(68)[sample(68)])
plot(areas, border = "black", add = T)

################################################################################
#### Load and Clean Simulated Trajectories
################################################################################
steps <- c(68, 125, 250, 500, 1000, 2000)
conns <- list()

# Loop through the step numbers
for (i in 1:length(steps)){

  # Load the dispersal trajectories
  sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

  # Subset to data of interest
  sims <- subset(sims
    , StepNumber    <= steps[i]
    & PointSampling == "Static"
  )

  ################################################################################
  #### DELETE
  ################################################################################
  # Subset data for now (5 trajectories per start point)
  set.seed(123)
  sims <- sims %>%
    group_by(StartPoint, ID) %>%
    nest() %>%
    group_by(StartPoint) %>%
    sample_n(5) %>%
    unnest()
  ################################################################################
  #### DELETE
  ################################################################################
  # Nest the data so that each trajectory (unique ID) receives its own row
  sims <- sims %>% group_by(ID) %>% nest()

  # Move from point to step representation (create lines from points)
  sims_traj <- mclapply(1:nrow(sims), mc.cores = detectCores() - 1, function(x){
    p <- SpatialPoints(sims$data[[x]][, c("x", "y")])
    l <- spLines(p)
    l$StartPoint <- sims$data[[x]]$StartPoint[1]
    return(l)
  }) %>% do.call(rbind, .)

  # Assign proper crs
  crs(sims_traj) <- CRS("+init=epsg:4326")

  # Assign original IDs back to lines
  sims_traj$ID <- sims$ID

  ################################################################################
  #### Connectivity Analysis: Approach I, Direct Connections
  ################################################################################
  # Identify StartPoint of each trajectory
  StartPoints <- suppressMessages(
    mclapply(sims$data, mc.cores = detectCores() - 1, function(x){
      x$StartPoint[1]
    }) %>% do.call(rbind, .)
  )

  # Identify EndPoint of each trajectory
  EndPoints <- suppressMessages(
    mclapply(sims$data, mc.cores = detectCores() - 1, function(x){
      last <- tail(x, 1)
      coordinates(last) <- c("x", "y")
      crs(last) <- CRS("+init=epsg:4326")
      endpoint <- raster::extract(areas_r, last)
      return(endpoint)
    }) %>% do.call(rbind, .)
  )

  # Put Start and EndPoints together and assign corresponding simulation ID
  direct_connections <- data.frame(
      TrackID     = sims$ID
    , StartPoint  = StartPoints
    , EndPoint    = EndPoints
  )
  conns[[i]] <- direct_connections
}
names(conns) <- steps
test <- bind_rows(conns, .id = "Steps")
test <- test %>%
  subset(StartPoint != EndPoint) %>%
  group_by(Steps, StartPoint, EndPoint) %>%
  na.omit() %>%
  summarize(TotalConns = n())
test %>% group_by(Steps) %>%
  summarize(MeanConns = mean(TotalConns)) %>%
  ggplot(aes(x = as.numeric(steps), y = MeanConns)) + geom_line() + geom_point()

################################################################################
#### Connectivity Analysis: Approach II, Indirect Connections
################################################################################
# Let's identify through which areas each trajectory moves.
visits <- suppressMessages(
  pbmclapply(
      X                   = 1:nrow(sims_traj)
    , mc.cores            = detectCores() - 1
    , ignore.interactive  = T
    , FUN                 = function(x){

      # Get visitation history of respective trajectory
      visits <- suppressWarnings(
        visitHist(line = sims_traj[x, ], raster = areas_r)
      )

      # In case there is a visitation history, return it
      if (!is.null(visits)){
        visits$ID <- sims_traj$ID[x]
        visits$StartPoint <- sims_traj$StartPoint[x]
        return(visits)

      # In case there is no visitation history, return NULL
      } else {
        return(NULL)
      }

  # Finally, we need to collapse the list
  }) %>% do.call(rbind, .)
)

# Save visitation history to file
# write_rds(visits, "03_Data/03_Results/00_Temporary_VisitationHistory.rds")
visits <- read_rds("03_Data/03_Results/00_Temporary_VisitationHistory.rds")

# We only want to keep unique entries and no trajectories that go to the source
# they came from
indirect_connections <- visits %>%
  distinct() %>%
  subset(StartPoint != to) %>%
  dplyr::select(-StartPoint) %>%
  rename(StartPoint = from, EndPoint = to) %>%
  rename(TrackID = ID) %>%
  dplyr::select(c(TrackID, everything())) %>%
  arrange(TrackID)

################################################################################
#### Combining Data from Approach I and Approach II
################################################################################
# Let's compare the two dataframes
direct_connections
indirect_connections

# Summarize data
connectivity_direct <- direct_connections %>%
  group_by(StartPoint, EndPoint) %>%
  summarize(TotalConnections = n())
connectivity_indirect <- indirect_connections %>%
  group_by(StartPoint, EndPoint) %>%
  summarize(TotalConnections = n())

# Create a dataframe for all potential connections (i.e. for all permutations of
# IDs). We will iteratively add information to this dataframe
connectivity <- points$ID %>%
  permutations(k = 2) %>%
  as.data.frame() %>%
  setNames(c("from", "to"))

# Combine the data from the direct and indirect connectivity approach
connectivity_direct <- left_join(
    connectivity
  , connectivity_direct
  , by = c("from" = "StartPoint", "to" = "EndPoint")
)
connectivity_indirect <- left_join(
    connectivity
  , connectivity_indirect
  , by = c("from" = "StartPoint", "to" = "EndPoint")
)
connectivity <- data.frame(
    from  = connectivity_direct$from
  , to    = connectivity_direct$to
  , DirectConnections = connectivity_direct$TotalConnections
  , IndirectConnections = connectivity_indirect$TotalConnections
)

# Replace NAs with 0s
connectivity <- connectivity %>%
  mutate(DirectConnections = replace_na(DirectConnections, 0)) %>%
  mutate(IndirectConnections = replace_na(IndirectConnections, 0))

################################################################################
#### Identify and Visualize Reached Areas
################################################################################
# Identify all areas that are reached from a specific source point
reached <- sims_traj %>%
  gIntersects(areas, ., byid = T) %>%
  as.data.frame() %>%
  setNames(areas$ID) %>%
  cbind(sims_traj@data, .) %>%
  dplyr::select(ID, everything()) %>%
  gather(key = EndPoint, value = Reached, 3:ncol(.)) %>%
  subset(Reached) %>%
  dplyr::select(-Reached) %>%
  arrange(StartPoint) %>%
  mutate(EndPoint = as.numeric(EndPoint))

# We may want to know how many dispersers actually stayed in each of these areas
stayed <- direct_connections %>%
  group_by(StartPoint, EndPoint) %>%
  summarize(TotalStayed = n()) %>%
  na.omit()

# Put the two dataframes together
reached <- stayed %>%
  left_join(reached, ., by = c("StartPoint", "EndPoint")) %>%
  mutate(TotalStayed = replace_na(TotalStayed, 0))

# Function to visualize reached areas
visCon <- function(n = 1, intensity = F, alpha = 0.5, tracks = F, exclude_source = F){

  # Get the source area
  source <- subset(areas, ID == n)

  # Identify all trajectories that are successful from this source area
  trajs <- subset(sims_traj, ID %in% reached$ID[
    reached$StartPoint == n & reached$StartPoint != reached$EndPoint
    ])

  # Get all areas that are reached from this source area
  reach <- subset(areas, ID %in% reached$EndPoint[reached$StartPoint == n])

  # Maybe the source area should be excluded?
  if (exclude_source){
    reach <- subset(reach, ID != n)
  }

  # Identify how often they were reached
  reach@data <- reached %>%
    subset(StartPoint == n) %>%
    dplyr::select(c(EndPoint, TotalStayed)) %>%
    distinct() %>%
    left_join(reach@data, ., by = c("ID" = "EndPoint"))

  # Visualize everything
  p <- tm_shape(areas) +
      tm_polygons(col = "gray90") +
      tm_text("ID") +
    tm_shape(reach) +
      tm_polygons(col = "TotalStayed", palette = "viridis", style = "cont") +
      tm_text("ID") +
    tm_shape(source) +
      tm_polygons(alpha = 0, border.col = "red", lwd = 2) +
      tm_text("ID")

  # In case tracks should be plotted, add them
  if (tracks){
    p <- p + tm_shape(trajs) + tm_lines(col = "red", alpha = alpha)
  }

  # plot
  p
}

# Try it out
visCon(n = 40, alpha = 0.1, tracks = F, exclude_source = F)
visCon(n = 40, alpha = 0.1, tracks = F, exclude_source = T)
visCon(n = 38, alpha = 0.02, tracks = T, exclude_source = F)

################################################################################
#### Graph Representation: Direct Connections
################################################################################
# Visualize the indirect network. Remove non-connections.
connectivity1 <- subset(connectivity, DirectConnections > 0)

# Identify areas with 0 connections (for plotting only)
islands1 <- subset(connectivity, DirectConnections > 0)
islands1 <- subset(areas, !(ID %in% islands1$from) & !(ID %in% islands1$to))

# Create a graph
net1 <- graph.data.frame(d = connectivity1, vertices = points$ID, directed = T)

# We can use attributes to specify node and edge sizes
E(net1)$width = log(E(net1)$DirectConnections) * 0.5

# Visualize. Note that we use the point-coordinates to get the correct layout
plot(net1
  , layout          = coordinates(points)
  , vertex.size     = 8
  , edge.arrow.size = 0.5
)

# We can actually plot the network graph on top of other spatial plots
plot(areas, col = viridis(20)[5], bg = "gray80")
plot(islands1, col = viridis(20)[15], add = T)
plot(net1
  , add             = T
  , layout          = coordinates(points)
  , rescale         = F
  , vertex.size     = 8
  , vertex.label    = NA
  , vertex.color    = "white"
  , edge.color      = "black"
  , edge.arrow.size = 0.3
  , edge.curved     = 0.3
)

# We could even identify modules (maybe we should use the directed graph?)
netundir1 <- as.undirected(net1, mode = "each")
clusters <- cluster_edge_betweenness(netundir1)
plot(areas, col = viridis(20)[5], bg = "gray80")
plot(islands1, col = viridis(20)[15], add = T)
plot(clusters, net1
  , add             = T
  , layout          = coordinates(points)
  , rescale         = F
  , vertex.size     = 8
  , vertex.label    = NA
  , edge.color      = "black"
  , edge.arrow.size = 0.3
  , edge.curved     = 0.3
)

# Let's check the adjacency matrix
as_adjacency_matrix(net1)

# Calculate the degree (number of connections)
plot(areas, col = "gray40", bg = "black")
plot(africa, add = T, border = "blue")
plot(net1
  , add             = T
  , layout          = coordinates(points)
  , rescale         = F
  , vertex.size     = 50
  , vertex.label    = degree(net1)
  , vertex.color    = c("orange", "red")[ifelse(degree(net1) == 0, 2, 1)]
  , edge.color      = "gray90"
  , edge.arrow.size = 0.1
  , edge.curved     = 0.3
)

# Calculate the closeness
plot(areas, col = "gray40", bg = "black")
plot(net1
  , add               = T
  , layout            = coordinates(points)
  , rescale           = F
  , vertex.size       = 50
  , vertex.label      = round(closeness(net1) / max(closeness(net1)), 2)
  , vertex.label.cex  = 0.7
  , vertex.color      = "orange"
  , edge.color        = "gray90"
  , edge.arrow.size   = 0.1
  , edge.curved       = 0.3
)

# Calculate the betweenness
plot(areas, col = "gray40", bg = "black")
plot(net1
  , add               = T
  , layout            = coordinates(points)
  , rescale           = F
  , vertex.size       = 50
  , vertex.label      = round(betweenness(net1) / max(betweenness(net1)), 2)
  , vertex.label.cex  = 0.7
  , vertex.color      = "orange"
  , edge.color        = "gray90"
  , edge.arrow.size   = 0.1
  , edge.curved       = 0.3
)

# Calculate the hub-score
plot(areas, col = "gray40", bg = "black")
plot(net1
  , add               = T
  , layout            = coordinates(points)
  , rescale           = F
  , vertex.size       = 50
  , vertex.label      = round(hub_score(net1)$vector, 2)
  , vertex.label.cex  = 0.7
  , vertex.color      = "orange"
  , edge.color        = "gray90"
  , edge.arrow.size   = 0.1
  , edge.curved       = 0.3
)

# Calculate the authority-score
plot(areas, col = "gray40", bg = "black")
plot(net1
  , add               = T
  , layout            = coordinates(points)
  , rescale           = F
  , vertex.size       = 50
  , vertex.label      = round(authority_score(net1)$vector, 2)
  , vertex.label.cex  = 0.7
  , vertex.color      = "orange"
  , edge.color        = "gray90"
  , edge.arrow.size   = 0.1
  , edge.curved       = 0.3
)

# Some graph level metrics
centralization.degree(net1)$centralization
centralization.closeness(net1)$centralization
centralization.betweenness(net1)$centralization
edge_density(net1)
reciprocity(net1)
transitivity(net1)
diameter(net1)
modularity(clusters)

################################################################################
#### Graph Representation: Indirect Connections
################################################################################
# Visualize the indirect network. Remove non-connections.
connectivity2 <- subset(connectivity, IndirectConnections > 0)

# Identify areas with 0 connections (for plotting only)
islands2 <- subset(connectivity2, IndirectConnections > 0)
islands2 <- subset(areas, !(ID %in% islands2$from) & !(ID %in% islands2$to))

# Create a graph
net2 <- graph.data.frame(d = connectivity2, vertices = points$ID, directed = T)

# We can use attributes to specify node and edge sizes
E(net2)$width = log(E(net2)$IndirectConnections) * 0.5

# Visualize. Note that we use the point-coordinates to get the correct layout
plot(net2
  , layout          = coordinates(points)
  , vertex.size     = hub_score(net2)$vector
  , edge.arrow.size = 0.5
)

# We can actually plot the network graph on top of other spatial plots
plot(areas, col = viridis(20)[5], bg = "gray80")
plot(islands2, col = viridis(20)[15], add = T)
plot(net2
  , add             = T
  , layout          = coordinates(points)
  , rescale         = F
  , vertex.size     = 8
  , vertex.label    = NA
  , vertex.color    = "white"
  , edge.color      = "black"
  , edge.arrow.size = 0.3
  , edge.curved     = 0.3
)

# We could even identify modules (maybe we should use the directed graph?)
netundir2 <- as.undirected(net2, mode = "each")
clusters <- cluster_edge_betweenness(netundir2)
clusters <- cluster_edge_betweenness(netundir2)
plot(areas, col = viridis(20)[5], bg = "gray80")
plot(islands2, col = viridis(20)[15], add = T)
plot(clusters, net2
  , add             = T
  , layout          = coordinates(points)
  , rescale         = F
  , vertex.size     = 8
  , vertex.label    = NA
  , edge.color      = "black"
  , edge.arrow.size = 0.3
  , edge.curved     = 0.3
)

################################################################################
#### Metrics
################################################################################
# We can derive some interesting metrics for each source area/source point
point_info <- data.frame(
  ID = points$ID
)

# For instance, we can identify how many dispersers leave a specific area
point_info <- direct_connections %>%
  subset(StartPoint != EndPoint) %>%
  group_by(StartPoint) %>%
  summarize(TotalLeft = n()) %>%
  left_join(point_info, ., by = c("ID" = "StartPoint")) %>%
  mutate(TotalLeft = replace_na(TotalLeft, 0))

# Conversely, we can check how many trajectories enter a specific area
point_info <- direct_connections %>%
  subset(StartPoint != EndPoint) %>%
  group_by(EndPoint) %>%
  summarize(TotalJoined = n()) %>%
  left_join(point_info, ., by = c("ID" = "EndPoint")) %>%
  mutate(TotalJoined = replace_na(TotalJoined, 0))

# We can also calculate the difference between the two
point_info$NetImmigration <- point_info$TotalJoined - point_info$TotalLeft

# Let's identify the how many other areas are reached from a specific source
# point or area (excluding the point where the trajectory came from)
point_info <- reached %>%
  dplyr::select(StartPoint, EndPoint) %>%
  group_by(StartPoint) %>%
  subset(StartPoint != EndPoint) %>%
  distinct() %>%
  summarize(NoAreasReached = n()) %>%
  left_join(point_info, ., by = c("ID" = "StartPoint")) %>%
  mutate(NoAreasReached = replace_na(NoAreasReached, 0))

# Let's calculate the total cumulative distance traveled along each trajectory
sims$CumulativeDistance <- suppressMessages(
  mclapply(1:nrow(sims), mc.cores = detectCores() - 1, function(x){
    x <- sims$data[[x]]
    distance <- sum(x$sl_, na.rm = T)
    return(distance)
  }) %>% do.call(rbind, .)
)

# Let's calculate the euclidean distance traveled along each trajectory
sims$EuclideanDistance <- suppressMessages(
  mclapply(1:nrow(sims), mc.cores = detectCores() - 1, function(x){
    x <- sims$data[[x]]
    x <- rbind(head(x, 1), tail(x, 1))
    coordinates(x) <- c("x", "y")
    crs(x) <- CRS("+init=epsg:4326")
    x <- spTransform(x, CRS("+init=epsg:32734"))
    x <- as.data.frame(x, xy = T)
    distance <- sqrt((x$x[2] - x$x[1]) ** 2 + (x$y[2] - x$y[1]) ** 2)
    return(distance)
  }) %>% do.call(rbind, .)
)

# Let's calculate the maximal distance traveled from each source point
sims$MaxDistance <- suppressMessages(
  mclapply(1:nrow(sims), mc.cores = detectCores() - 1, function(x){
    x <- sims$data[[x]]
    coordinates(x) <- c("x", "y")
    crs(x) <- CRS("+init=epsg:4326")
    x <- spTransform(x, CRS("+init=epsg:32734"))
    x <- as.data.frame(x, xy = T)
    distance <- sqrt((x$x - x$x[1]) ** 2 + (x$y - x$y[1]) ** 2)
    distance <- max(distance)
    return(distance)
  }) %>% do.call(rbind, .)
)

# Calculate average cumulative distance traveled by source point/source area
point_info <- data.frame(
      StartPoint          = sapply(sims$data, function(x){x$StartPoint[[1]]})
    , CumulativeDistance  = sims$CumulativeDistance
  ) %>%
  group_by(StartPoint) %>%
  summarize(MeanCumulativeDistance = round(mean(CumulativeDistance))) %>%
  left_join(point_info, ., by = c("ID" = "StartPoint"))

# Calculate average Euclidean distance to source point
point_info <- data.frame(
      StartPoint  = sapply(sims$data, function(x){x$StartPoint[[1]]})
    , EuclideanDistance = sims$EuclideanDistance
  ) %>%
  group_by(StartPoint) %>%
  summarize(MeanEuclideanDistance = round(mean(EuclideanDistance))) %>%
  left_join(point_info, ., by = c("ID" = "StartPoint"))

# Calculate average maximal distance to source point
point_info <- data.frame(
      StartPoint  = sapply(sims$data, function(x){x$StartPoint[[1]]})
    , MaxDistance = sims$MaxDistance
  ) %>%
  group_by(StartPoint) %>%
  summarize(MeanMaxDistance = round(mean(MaxDistance))) %>%
  left_join(point_info, ., by = c("ID" = "StartPoint"))

# Join all metrics to the areas and points
areas@data <- left_join(areas@data, point_info, by = "ID")
points@data <- left_join(points@data, point_info, by = "ID")

# Write stuff to file
write_rds(connectivity, "03_Data/03_Results/00_Temporary_Connectivity.rds")
write_rds(areas, "03_Data/03_Results/00_Temporary_Areas.rds")
write_rds(points, "03_Data/03_Results/00_Temporary_Points.rds")
connectivity <- read_rds("03_Data/03_Results/00_Temporary_Connectivity.rds")
areas <- read_rds("03_Data/03_Results/00_Temporary_Areas.rds")
points <- read_rds("03_Data/03_Results/00_Temporary_Points.rds")

############################################################
#### Visualize Metrics
############################################################
# Plot the number of areas reached
p1 <- tm_shape(areas) +
  tm_polygons("NoAreasReached", palette = "viridis", border.col = "black", style = "cont") +
  tm_text("NoAreasReached", size = 0.7)

# Plot the number of dispersers that left
p2 <- tm_shape(areas) +
  tm_polygons("TotalLeft", palette = "viridis", border.col = "black", style = "cont") +
  tm_text("TotalLeft", size = 0.7)

# Plot the number of dispersers that joined
p3 <- tm_shape(areas) +
  tm_polygons("TotalJoined", palette = "viridis", border.col = "black", style = "cont") +
  tm_text("TotalJoined", size = 0.7)

# Plot the number of dispersers that joined
p4 <- tm_shape(areas) +
  tm_polygons("NetImmigration", palette = "viridis", border.col = "black", style = "cont") +
  tm_text("NetImmigration", size = 0.7)

# Plot the mean cumulative distance
p5 <- tm_shape(areas) +
  tm_polygons("MeanCumulativeDistance", palette = "viridis", border.col = "black", style = "cont") +
  tm_text("MeanCumulativeDistance", size = 0.7)

# Plot the mean mean euclidean distance
p6 <- tm_shape(areas) +
  tm_polygons("MeanEuclideanDistance", palette = "viridis", border.col = "black", style = "cont") +
  tm_text("MeanEuclideanDistance", size = 0.7)

# Plot the mean mean euclidean distance
p7 <- tm_shape(areas) +
  tm_polygons("MeanMaxDistance", palette = "viridis", border.col = "black", style = "cont") +
  tm_text("MeanMaxDistance", size = 0.7)

# Put all plots together
tmap_arrange(p1, p2, p3, p4, p5, p6, p7)

############################################################
#### Connectivity Analysis: Successfull Events
############################################################
# Cut successful connections at the closest point to their destination
dat$LinesCut <- lapply(1:length(dat$Lines), function(x){
  mclapply(1:length(dat$Lines[[x]]), mc.cores = detectCores() - 1, function(y){
    cutLineClose(
        line    = dat$Lines[[x]][y, ]
      , point   = dat$Destin_point[[x]]
      , polygon = dat$Destin_area[[x]]
    )
  }) %>% do.call(rbind, .)
})

# Visualize some of the connections
n <- 150
plot(areas)
plot(dat$LinesCut[[n]], add = T, col = "red")
plot(dat$Origin_point[[n]], add = T, col = "blue")
plot(dat$Destin_point[[n]], add = T, col = "blue")

# Check the number of connections and which points they link
connections_0 <- connectivity %>%
  na.omit() %>%
  group_by(StartPoint, EndPoint) %>%
  summarize(TotalConnections = sum(Reached))
connections_0

# Show only those points that actually have a connection
connections_0 %>%
  subset(TotalConnections > 0) %>%
  arrange(-TotalConnections)

# Check number of connections that each start point has
connections_1 <- connections_0 %>%
  group_by(StartPoint) %>%
  summarize(TotalConnections = sum(TotalConnections))
connections_1 %>% arrange(-TotalConnections)

# Check number of connections that each end point has
connections_2 <- connections_0 %>%
  group_by(EndPoint) %>%
  summarize(TotalConnections = sum(TotalConnections))
connections_2 %>% arrange(-TotalConnections)

# Identify all areas that neither produce nor receive connections
loosers_1 <- subset(connections_1, TotalConnections == 0)
loosers_2 <- subset(connections_2, TotalConnections == 0)
isolated <- subset(areas
  , ID %in% loosers_1$StartPoint
  & ID %in% loosers_2$EndPoint
)

# Visualize them
plot(areas)
plot(isolated, add = T, col = "red")

# Could also identify sources and sinks
sources <- areas
sources@data <- connections_1 %>%
  rename(ID = StartPoint) %>%
  left_join(areas@data, ., by = "ID")
sinks <- areas
sinks@data <- connections_2 %>%
  rename(ID = EndPoint) %>%
  left_join(areas@data, ., by = "ID")

# Visualize
spplot(sources, "TotalConnections")
spplot(sinks, "TotalConnections")
