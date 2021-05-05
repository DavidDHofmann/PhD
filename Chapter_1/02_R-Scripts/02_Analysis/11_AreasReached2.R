################################################################################
#### Cluster Analysis of Simulated Connections
################################################################################
# Description: In this script, we overlay the study area with a regular raster
# and generate a network which allows us to do a simple cluster analysis

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
library(dismo)          # To calculate voronoi polygons

################################################################################
#### Create Regular Grid
################################################################################
# Load areas from which we released dispersers
areas <- readOGR("03_Data/03_Results/99_SourceAreas.shp")

# Load reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Create points that are spaced regularly at 50 km. We use a raster with desired
# grid-size for this purpose
points <- r %>%
  aggregate(., fact = (50000) / 250, fun = mean) %>%
  rasterToPoints(., spatial = TRUE)

# Keep only those points that lie within protected areas
points <- points[areas, ]

# Find all source areas that do not contain any points yet, we'll need to assign
# seperate points to them
missing <- areas %>%
  gIntersects(., points, byid = TRUE) %>%
  colSums(.) %>%
  as.vector(.)

# Identify those polygons with zero points
missing <- missing == 0

# Create centroids for these polygons, but make sure they are WITHIN the
# polygon. We can use our custom function for this.
centroids <- gCentroidWithin(areas[missing, ])
points <- rbind(SpatialPoints(points), SpatialPoints(centroids))
crs(points) <- crs(areas)

# Plot them
plot(areas)
plot(points, add = T, pch = 20, cex = 0.4)
text(points, "ID", cex = 0.7, halo = T, pos = 3)

# Also indicate to which area they belong to
points$ID <- over(points, areas)$ID

# Plot to verify that IDs match protected area
plot(areas, border = areas$ID)
plot(points, col = points$ID, add = T, pch = 20, cex = 0.4)

# We now want to cut the areas into catchment areas such that each source point
# is allocated to one polygon. We can create voronoi polygons for this purpose.
areas <- lapply(1:length(unique(areas$ID)), function(x){

  # Subset to respective areas and points
  areas_sub <- subset(areas, ID == x)
  points_sub <- subset(points, ID == x)

  # Can only tesselate the source area if it contains more than 1 point
  if (nrow(points_sub) > 1){

    # Create voronoi polygons based on the selected points. We make the extent
    # larger to make sure the entire source_area is covered by the tesselated
    # polygon.
    voris <- voronoi(points_sub, ext = extent(areas_sub) + c(-1, 1, -1, 1))

    # For each voronoi polygon we now create a clipped area
    areas_sub <- gIntersection(voris, areas_sub, byid = T)
  }

  # Return the clipped area
  return(as(areas_sub, "SpatialPolygons"))

}) %>% do.call(rbind, .)

# Visualize the cut areas
plot(areas)
plot(points, add = T, pch = 20, cex = 0.4)

# Assign ID to each area, corresponding to the ID of the point that it contains.
areas$ID <- areas %>%
  over(points) %>%
  dplyr::select(ID) %>%
  as.matrix() %>%
  as.vector()

# Rasterize the area IDs
areas_r <- raster(terra::rasterize(x = vect(areas), y = rast(r), field = "ID"))

# Visualize the results and visually verify that IDs were assigned correctly
plot(areas, border = areas$ID, col = adjustcolor(areas$ID, alpha.f = 0.3))
plot(points, col = points$ID, add = T, pch = 20, cex = 0.4)

################################################################################
#### Identify Connections
################################################################################
# Load simulations
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")

# Create a spatial point for the first location of each trajectory
first <- sims %>%
  dplyr::select("x", "y", "TrackID") %>%
  group_by(TrackID) %>%
  slice(1) %>%
  SpatialPointsDataFrame(
      coords      = cbind(.[["x"]], .[["y"]])
    , proj4string = CRS("+init=epsg:4326")
  )

# Assess the ID of the area from which each trajectory leaves
first$From <- as.numeric(over(first, areas)$ID)
first <- first@data[, c("TrackID", "From")]

# Check it
head(first, 10)

# Join this information to the simulated tracks
sims <- left_join(sims, first, by = "TrackID")

# Check it
head(sims, 10)

# We only care about simulations leaving from one of these areas (i.e. come from
# the main area)
sims <- subset(sims, Area == "Main")
sims <- subset(sims, !is.na(From))

# Identify number of trajectories leaving from each area
nsims <- sims %>%
  group_by(TrackID, From) %>%
  nest() %>%
  ungroup() %>%
  count(From) %>%
  arrange(From) %>%
  setNames(c("From", "Simulations"))

# We can now verify that we have all simulations
sum(nsims$Simulations)

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# Identify through which areas each disperser moved
visits <- data.frame(
    TrackID    = sims$TrackID
  , StepNumber = sims$StepNumber
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , Area       = raster::extract(areas_r, sims)
)

# Add this information to the simulations
sims$To <- visits$Area

# Remove spatial data
sims <- as.data.frame(sims, xy = T)
sims$xy <- NULL

# Identify how long it takes on average to reach the different areas
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

# Add the information about how many dispersers we simulated from each area
visits <- left_join(visits, nsims, by = "From")

# Calculate relative frequency
visits$RelFrequency <- visits$Frequency / visits$Simulations

# Store the visits to file
write_rds(visits, "03_Data/03_Results/99_AreasReached2.rds")
writeOGR(areas
  , dsn       = "03_Data/03_Results"
  , layer     = "99_NetworkAreas"
  , driver    = "ESRI Shapefile"
  , overwrite = T
)
writeOGR(points
  , dsn       = "03_Data/03_Results"
  , layer     = "99_NetworkPoints"
  , driver    = "ESRI Shapefile"
  , overwrite = T
)

################################################################################
#### Network View
################################################################################
# Coerce the visitations to a graph
net <- graph_from_data_frame(
    d        = visits
  , vertices = unique(areas$ID)
  , directed = T
)

# Check attributes
vertex_attr(net)
edge_attr(net)

# Prepare layout
lay <- coordinates(gCentroid(areas, byid = T))

# Plot
plot(net
  , layout             = lay
  , vertex.label       = NA
  , vertex.size        = 1
  , edge.arrow.size    = 0
  , edge.curved        = 0.1
  , edge.width         = E(net)$RelFrequency
  , edge.color         = "orange"
  , vertex.frame.color = NA
  , vertex.color       = "black"
)

# Main Plot
library(ggnetwork)
library(viridis)
r <- as.data.frame(r, xy = T)
pal <- colorRampPalette(plasma(100, begin = 0.8, end = 0))
net_p <- ggnetwork(net, layout = lay, arrow.gap = 0.1, scale = F)
ggplot() +
  geom_edges(
      data      = net_p
    , mapping   = aes(x = x, y = y, xend = xend, yend = yend
      , size = RelFrequency, col = MeanStepNumber)
    , curvature = 0.2
    , arrow     = arrow(length = unit(0, "pt"), type = "closed", angle = 10)
  ) +
  scale_size_area(
      name     = "Visitation Frequency"
    , max_size = 0.15
  ) +
  scale_color_gradientn(
      colors  = pal(100)
    , guide   = guide_colorbar(
        title          = "Number of Steps"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(3.0, "cm")
    )
  ) +
  scale_fill_manual(
    values = c("#70ab70", "#d9f0d3")
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    , title    = "Areas Reached and Visitation Frequency"
    , subtitle = "In Relation to Number of Steps"
  ) +
  guides(
    size = guide_legend(title.position = "top")
  ) +
  theme(
      legend.position      = "bottom"
    , legend.box           = "horizontal"
    , legend.title.align   = 0.5
    , panel.background     = element_blank()
    , panel.border         = element_rect(colour = "black", fill = NA, size = 1)
    , legend.title         = element_text(size = 10),
    , legend.text          = element_text(size = 8)
    , legend.margin        = margin(c(0, 0, 0, 0))
  )

# Try to cluster the nodes
cluster_betwee <- cluster_edge_betweenness(
    graph   = as.undirected(net)
  , weights = E(net)$RelFrequency
)

# Alternatively, we can identify communities based on greedy optimization
cluster_greedy <- cluster_fast_greedy(
    graph   = as.undirected(net)
  , weights = E(net)$RelFrequency
)

# Visualize the community
plot(areas, main = "Cluster Betweenness")
  plot(cluster_betwee
    , net
    , layout             = lay
    , vertex.label       = NA
    , vertex.size        = 1
    , edge.arrow.size    = 0
    , edge.curved        = 0.1
    , edge.color         = NA
    , vertex.frame.color = NA
    , add                = T
    , rescale            = F
  )

plot(areas, main = "Cluster Greedy")
  plot(cluster_greedy
    , net
    , layout             = lay
    , vertex.label       = NA
    , vertex.size        = 1
    , edge.arrow.size    = 0
    , edge.curved        = 0.1
    , edge.color         = NA
    , vertex.frame.color = NA
    , add                = T
    , rescale            = F
  )







# Prepare plot for ggplotting
net_p <- ggnetwork(net, layout = lay, arrow.gap = 0)

# Plot
ggplot(net_p, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(col = MeanStepNumber, size = RelFrequency), curvature = 0.2) +
  geom_nodes(color = "orange") +
  scale_color_continuous() +
  theme_blank()
