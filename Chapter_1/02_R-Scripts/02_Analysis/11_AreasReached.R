################################################################################
#### Number of Connections Between Protected Areas
################################################################################
# Description: In this script, we analyze the simulated dispersal trajectories
# and identify the number of connections between protected areas

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(pbmcapply)      # To run on multiple cores with progress bar
library(rgeos)          # For spatial ananylsis
library(igraph)         # For network analysis
library(davidoff)       # Custom functions
library(rgdal)          # To load spatial data

################################################################################
#### Rasterize National Parks
################################################################################
# Load protected areas
prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Load reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Keep only national parks
prot <- subset(prot, Desig == "National Park")

# Assign a unique ID to each of the areas
prot$ID <- 1:nrow(prot)

# Rasterize IDs
prot_r <- raster(terra::rasterize(x = vect(prot), y = rast(r), field = "ID"))

# Visualize them
plot(prot_r, main = "National Parks")
text(prot, "Name", cex = 0.5)

################################################################################
#### Identify Connections
################################################################################
# Load simulations
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")

# Create SpatialPoints for first location of each trajectory
first <- sims %>%
  dplyr::select("x", "y", "TrackID") %>%
  group_by(TrackID) %>%
  slice(1) %>%
  SpatialPointsDataFrame(
      coords      = cbind(.[["x"]], .[["y"]])
    , proj4string = CRS("+init=epsg:4326")
  )

# Assess from which protected area each trajectory leaves
first$From <- as.numeric(over(first, prot)$ID)
first <- first@data[, c("TrackID", "From")]

# Join information to simulated tracks
sims <- left_join(sims, first, by = "TrackID")

# We only care about simulations leaving from national parks
sims <- subset(sims, !is.na(From))

# Identify number of trajectories leaving from each national park
nsims <- sims %>%
  group_by(TrackID, From) %>%
  nest() %>%
  ungroup() %>%
  count(From) %>%
  arrange(From) %>%
  setNames(c("From", "Simulations"))

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# Identify through which national parks the dispersers moved
visits <- data.frame(
    TrackID    = sims$TrackID
  , StepNumber = sims$StepNumber
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , Prot       = raster::extract(prot_r, sims)
)

# Add this information to the simulations
sims$To <- visits$Prot

# Remove spatial data
sims <- as.data.frame(sims, xy = T)
sims$xy <- NULL

# Identify how long it takes on average to reach different areas
visits <- sims %>%
  group_by(TrackID, From, To) %>%
  summarize(StepNumber = min(StepNumber), .groups = "drop") %>%
  subset(!is.na(From) & !is.nan(To)) %>%
  arrange(TrackID, StepNumber) %>%
  group_by(From, To) %>%
  summarize(
      MeanStepNumber = mean(StepNumber)
    , SDStepNumber   = sd(StepNumber)
    , Frequency      = n()
    , .groups        = "drop"
  ) %>%
  subset(From != To)

# Join the number of dispersers
visits <- left_join(visits, nsims, by = "From")

# Calculate relative frequency
visits$RelFrequency <- visits$Frequency / visits$Simulations

# Store areas reached to file
write_rds(visits, "03_Data/03_Results/99_AreasReached.rds")

################################################################################
#### Network View
################################################################################
# Coerce the visitations to a graph
net <- graph_from_data_frame(
    d        = visits
  , vertices = unique(prot$ID)
  , directed = T
)
vertex_attr(net)
edge_attr(net)

# Add area as vertex information
V(net)$Area <- prot$Area

# Prepare layout
lay <- coordinates(gCentroid(prot, byid = T))

# Prepare plot for ggplotting
net_p <- ggnetwork(net, layout = lay, arrow.gap = 0)

# Plot
ggplot(net_p, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(col = MeanStepNumber, size = RelFrequency), curvature = 0.2) +
  geom_nodes(color = "orange") +
  scale_color_continuous() +
  theme_blank()
