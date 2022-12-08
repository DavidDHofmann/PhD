############################################################
#### Calculating Least Cost Paths and Corridors
############################################################
# Description: Based on the permeability map derived in the previous script we
# can now identify least cost paths and corridors between sites of interst. We
# will follow two different approaches two select these source points. (1) We
# will define source points within protected areas large enough to sustain
# viable wild dog populations. (2) We will generate source points distributed
# along the borders of the study extent.

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Set a seed for reproducability
set.seed(12345)

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(raster)
library(data.table)
library(rgeos)
library(tidyverse)
library(viridis)
library(rasterVis)
library(gdistance)
library(rgdal)
library(arrangements)
library(tictoc)
library(measurements)
library(tictoc)
library(parallel)
library(davidoff)
library(igraph)

################################################################################
#### Define Start Points: Approach I
################################################################################
# Reload the permeability map that we derived in the last script
permeability <- raster("03_Data/03_Results/99_PermeabilityMap.tif")

# Load kaza shapefile for plotting
kaza <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")

# Load protected areas
prot <- shapefile(
  "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(1Class).shp"
)

# Dissolve the borders
prot <- prot %>%

  # Flatten all polygons into a single object
  aggregate(., dissolve = TRUE) %>%

  # The geometry is now invalid. We apply a tiny buffer (1m) to make it valid
  gBuffer(., width = 1 / 110 * 0.001) %>%

  # Disaggregate polygons again, so that each object gets its own slot
  disaggregate(.) %>%

  # Convert to SpatialPolygonsDataFrame
  as(., "SpatialPolygonsDataFrame") %>%

  # Crop to original extent
  crop(., permeability)

# Calculate area of each polygon
prot$Area <- prot %>%

  # Need to transform to utm to be able to calculate area in km
  spTransform(., CRS("+init=epsg:32734")) %>%

  # Get the area
  gArea(., byid = TRUE) %>%

  # Convert to km**2
  conv_unit(., from = "m2", to = "km2")

# Retain only polygons with an area of over 700km2
plot(subset(prot, Area >= 700))
prot <- subset(prot, Area >= 700)

# Assign unique ID to each polygon
prot$ID <- 1:nrow(prot)
prot$ID <- as.numeric(prot$ID)

# Rearrange columns
prot@data <- dplyr::select(prot@data, c(ID, Area))

# Create points spaced regularly at 100 km. We use a raster with desired
# grid-size for this purpose
points1 <- permeability %>%
  aggregate(., fact = (100 * 1000) / 250, fun = mean) %>%
  rasterToPoints(., spatial = TRUE)

# Keep only those points that lie within protected areas
index <- gContains(prot, points1, byid = TRUE) %>% rowSums(.)
points1 <- points1[index > 0, ]

# Find all polygons in which there is no point yet
missing <- prot %>%

  # Check intersecting points and polygons
  gIntersects(., points1, byid = TRUE) %>%

  # Identify number of points in each polygon
  colSums(.) %>%

  # Convert result to vector
  as.vector(.)

# Identify those polygons with zero points
missing <- missing == 0

# Create centroids for these polygons, but make sure they are WITHIN the
# polygon. We can use our custom function for this.
centroids <- gCentroidWithin(prot[missing, ])

# Put the centroids and the other points together
points1 <- rbind(SpatialPoints(points1), SpatialPoints(centroids))
crs(points1) <- crs(prot)

# Assign unique IDs to each point
points1$ID <- 1:length(points1)
points1$ID <- as.numeric(points1$ID)

# Also indicate to which protected area they belong to
together <- gContains(prot, points1, byid = T)
points1$ProtectedArea <- apply(together, 1, which)
points1$ProtectedArea <- as.numeric(points1$ProtectedArea)

# Plot to verify that IDs match protected area
plot(prot, border = rainbow(length(unique(prot$ID)))[prot$ID])
plot(points1, col = rainbow(length(unique(prot$ID)))[points1$ProtectedArea], add = T)

# Let's check the number of combinations we need to run
nrow(combinations(nrow(points1), 2))

# Store the source points and areas
writeOGR(prot
  , dsn       = "03_Data/03_Results"
  , layer     = "99_SourceAreas"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
writeOGR(points1
  , dsn       = "03_Data/03_Results"
  , layer     = "99_SourcePoints1"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

################################################################################
#### Define Start Points: Approach II
################################################################################
# We want to distribute source points as close to the borders as we can, yet we
# need to make sure they lie WITHIN the study extent. To achieve this, we will
# distribute them along the extent of a raster that was cropped minimally.
r <- permeability
r[1, ] <- NA
r[nrow(r), ] <- NA
r[, 1] <- NA
r[, ncol(r)] <- NA

# Trim the raster and retrieve the extent
r <- trim(r)
ext <- as(extent(r), "SpatialLines")
crs(ext) <- crs(r)

# Distribute points along the extent
points2 <- spsample(ext, type = "regular", n = 68)
plot(points2)

# Add point ids
points2$ID <- 1:length(points2)

# Visualize them
plot(prot)
plot(points2, add = T)

# Store them
writeOGR(points2
  , dsn       = "03_Data/03_Results"
  , layer     = "99_SourcePoints2"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

################################################################################
#### Reload Points
################################################################################
points1 <- readOGR("03_Data/03_Results/99_SourcePoints1.shp")
points2 <- readOGR("03_Data/03_Results/99_SourcePoints2.shp")
plot(prot)
plot(points1, add = T, col = "blue")
plot(points2, add = T, col = "purple")

################################################################################
#### Visualize Connections
################################################################################
# Identify all possible connections for "points1"
combis1 <- nrow(points1) %>% combinations(., 2) %>%
  as.data.frame(.) %>%
  rename(., Origin = V1, Destin = V2) %>%
  mutate(Connectiosn = 1)

# Identify all possible connections for "points2"
combis2 <- nrow(points2) %>% combinations(., 2) %>%
  as.data.frame(.) %>%
  rename(., Origin = V1, Destin = V2) %>%
  mutate(Connectiosn = 1)

# Let's use a network graph to show the connections
graph1 <- graph_from_data_frame(combis1, vertices = points1$ID)
plot(graph1
  , vertex.size     = 4
  , vertex.label    = NA
  , edge.width      = 0.1
  , edge.arrow.size = 0.1
  , layout          = coordinates(points1)
)
graph2 <- graph_from_data_frame(combis2, vertices = points2$ID)
plot(graph2
  , vertex.size     = 4
  , vertex.label    = NA
  , edge.width      = 0.1
  , edge.arrow.size = 0.1
  , layout          = coordinates(points2)
)

############################################################
#### Prepare Transition Layer
############################################################
# For least cost path/corridor modelling we now need to convert the permeability
# map to a "transition" layer that allows us to apply graph theory. To calculate
# the edge weights between adjacent cells we use their mean, i.e. f(x, y) =
# mean(cost(x), cost(y))
trans <- transition(permeability, directions = 8, transitionFunction = mean)

# We also need to apply a geometric correction that accounts for the distances
# from the cell center to the cell corners. For details and explanations you can
# check in the book "Spatial Ecology and Conservation Modeling" chapter 9.3
trans <- geoCorrection(trans, type = "c", multpl = FALSE)

# For some reason this operation leaves a ton of junk in the ram. Let's clean it
# by calling gc() (garbage-collector)
gc()

############################################################
#### Least Cost Paths
############################################################
# Now find the shortest paths for all possible connections. I follow a slightly
# weird procedure since I try to avoid looping over all origins separately and
# rather cluster them (that is, I calculate from point 1 to all other points in
# one single step). This strictly increases performance. I tried to run it in
# parallel, but it keeps causing a memory overflow.
shortest1 <- suppressMessages(
  lapply(1:length(unique(combis1$Origin)), function(x){
    shortestPath(trans
      , origin  = points1[x, ]
      , goal    = points1[combis1$Destin[combis1$Origin == x], ]
      , output  = "SpatialLines"
    )
  })
)
shortest1 <- do.call(rbind, shortest1)
shortest1$ID <- 1:length(shortest1)

# Repeat for the second set of source points
shortest2 <- suppressMessages(
  lapply(1:length(unique(combis2$Origin)), function(x){
    shortestPath(trans
      , origin  = points2[x, ]
      , goal    = points2[combis2$Destin[combis2$Origin == x], ]
      , output  = "SpatialLines"
    )
  })
)
shortest2 <- do.call(rbind, shortest2)
shortest2$ID <- 1:length(shortest2)

# Plot the least cost paths on top of the permeability map
plot(permeability, col = viridis(50))
lines(shortest1, col = "white")
points(points1, col = "white")

plot(permeability, col = viridis(50))
lines(shortest2, col = "white")
points(points2, col = "white")

# Write the least cost paths to file
writeOGR(shortest1
  , "03_Data/03_Results"
  , "99_LeastCostPaths1"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

# Write the least cost paths to file
writeOGR(shortest2
  , "03_Data/03_Results"
  , "99_LeastCostPaths2"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

# Read the paths back in
shortest1 <- readOGR("03_Data/03_Results/99_LeastCostPaths1.shp")
shortest2 <- readOGR("03_Data/03_Results/99_LeastCostPaths2.shp")

################################################################################
#### Rasterize LCPs
################################################################################
# Buffer paths
shortest_buff1 <- gBuffer(shortest1
  , byid  = TRUE
  , width = metersToDegrees(2300)
)

# Buffer paths
shortest_buff2 <- gBuffer(shortest2
  , byid  = TRUE
  , width = metersToDegrees(2300)
)

# Rasterize paths
shortest_rast1 <- rasterizeTerra(
    vect(shortest_buff1)
  , rast(permeability)
)

# Rasterize paths
shortest_rast2 <- rasterizeTerra(
    vect(shortest_buff2)
  , rast(permeability)
)

# Plot
plot(shortest_rast1, col = viridis(50))
plot(shortest_rast2, col = viridis(50))
plot(points2, add = T, col = "white", pch = 20)

# Store the rasterized LCPs
writeRaster(
    shortest_rast1
  , "03_Data/03_Results/99_LeastCostPaths1.tif"
  , overwrite = TRUE
)

# Store the rasterized LCPs
writeRaster(
    shortest_rast2
  , "03_Data/03_Results/99_LeastCostPaths2.tif"
  , overwrite = TRUE
)

############################################################
#### Least Cost Corridors: Points I
############################################################
# Create a list of points
points_l <- splitShape(points1, n = length(points1))

# Compute the cumulative cost map for each point
cost <- lapply(points_l, function(x){
  costi <- accCost(trans, x)
  costi <- writeRaster(costi, tempfile())
  return(costi)
  gc()
})

# Load all costmaps into a single stack
cost <- stack(cost)
cost <- stackSave(cost, tempfile())

# Prepare an empty raster onto which we will add all least cost corridors
r0 <- raster(cost[[1]])
values(r0) <- 0

# Loop through all combinations and calculate a least cost corridor
for (i in 1:nrow(combis1)){ <- 0

  # Sum cost maps of the two poitns
  cost1 <- cost[[combis1[i, "Origin"]]]
  cost2 <- cost[[combis1[i, "Destin"]]]
  corr <- calc(stack(cost1, cost2), sum)

  # Calculate threshold below which the costs are not higher than + x% of the
  # minimal costs and reclassify raster
  threshold <- minValue(corr) * 1.05
  corr <- reclassify(corr, c(threshold, Inf, NA))

  # Scale the costs to 0 and 1
  min <- minValue(corr)
  max <- maxValue(corr)
  corr <- calc(corr, function(x){(x - min) / (max - min)})

  # Assign all NAs a value of 1
  corr <- reclassify(corr, c(NA, NA, 1), right = FALSE)

  # Add the corridor to the r0 raster
  r0 <- calc(stack(r0, corr), sum)

  # Collect garbage
  gc()
}

# Rescale values to between 0 and 1 and invert them
min <- minValue(r0)
max <- maxValue(r0)
r0 <- calc(r0, function(x){(x - min) / (max - min)})
r0 <- r0 * (-1) + 1

# Plot the result
plot(r0, col = viridis(50))
plot(kaza, add = T, border = "white")

# Write result to file
writeRaster(r0
  , "03_Data/03_Results/99_LeastCostCorridors1.tif"
  , overwrite = TRUE
)

############################################################
#### Least Cost Corridors: Points II
############################################################
# Create a list of points
points_l <- splitShape(points2, n = length(points2))

# Compute the cumulative cost map for each point
cost <- lapply(points_l, function(x){
  costi <- accCost(trans, x)
  costi <- writeRaster(costi, tempfile())
  return(costi)
  gc()
})

# Load all costmaps into a single stack
cost <- stack(cost)
cost <- stackSave(cost, tempfile())

# Prepare an empty raster onto which we will add all least cost corridors
r0 <- raster(cost[[1]])
values(r0) <- 0

# Loop through all combinations and calculate a least cost corridor
for (i in 1:nrow(combis2)){

  # Sum cost maps of the two poitns
  cost1 <- cost[[combis2[i, "Origin"]]]
  cost2 <- cost[[combis2[i, "Destin"]]]
  corr <- calc(stack(cost1, cost2), sum)

  # Calculate threshold below which the costs are not higher than + x% of the
  # minimal costs and reclassify raster
  threshold <- minValue(corr) * 1.05
  corr <- reclassify(corr, c(threshold, Inf, NA))

  # Scale the costs to 0 and 1
  min <- minValue(corr)
  max <- maxValue(corr)
  corr <- calc(corr, function(x){(x - min) / (max - min)})

  # Assign all NAs a value of 1
  corr <- reclassify(corr, c(NA, NA, 1), right = FALSE)

  # Add the corridor to the r0 raster
  r0 <- calc(stack(r0, corr), sum)

  # Collect garbage
  gc()
}

# Rescale values to between 0 and 1 and invert them
min <- minValue(r0)
max <- maxValue(r0)
r0 <- calc(r0, function(x){(x - min) / (max - min)})
r0 <- r0 * (-1) + 1

# Plot the result
plot(r0, col = viridis(50))
plot(kaza, add = T, border = "white")

# Write result to file
writeRaster(r0
  , "03_Data/03_Results/99_LeastCostCorridors2.tif"
  , overwrite = TRUE
)

################################################################################
#### Pearson Correlation
################################################################################
# Reload LCC maps
lcc1 <- raster("03_Data/03_Results/99_LeastCostCorridors1.tif")
lcc2 <- raster("03_Data/03_Results/99_LeastCostCorridors2.tif")

# Plot both
par(mfrow = c(2, 1))
plot(lcc1, col = viridis(50))
plot(kaza, add = T, border = "white")
plot(lcc2, col = viridis(50))
plot(kaza, add = T, border = "white")

# Check correlation
cor(values(lcc1), values(lcc2))
