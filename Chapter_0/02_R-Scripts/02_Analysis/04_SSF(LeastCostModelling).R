############################################################
#### Calculating Least Cost Paths and Corridors
############################################################
# Description: Based on the permeability map derived in the previous script we
# can now identify least cost paths or corridors between sites of interst. In
# our case, we are interested in assessing the connectivity between national
# parks

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

# Reload the permeability map that we derived in the last script
permeability <- "03_Data/03_Results/99_PermeabilityMap.tif" %>% raster()

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
#### Define Start and Endpoints: Method I
############################################################
# We want to identify dispersal corridors between national parks. Let's
# therefore load the shapefiles containing all national parks in the KAZA
prot <- "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(3Classes).shp" %>%
  shapefile()

# Subset to national parks only
prot <- subset(prot, Desig == "National Park")

# Remove some of the national parks that we don't want to consider
prot <- prot[-which(prot$Name == "Lusaka"), ]
prot <- prot[-which(prot$Name == "Kafue")[1], ]
prot <- prot[-which(prot$Name == "Northern Tuli"), ]
prot <- prot[-which(prot$Name == "Ngonye Falls"), ]
prot <- prot[-which(prot$Name == "CT/12 (Soad Ash Botswana & Nata Sanctuary)"), ]
prot <- prot[-which(prot$Name == "CT/22"), ]

# Let's plot the remaining parks
plot(prot)

# Specify points that we will connect and calculate least cost paths/corridors
points <- gCentroid(prot, byid = TRUE)

# Note that we could also sample multiple points within each protected area
# However, this would severely increase the computational requirements. For
# instance, if we were to sample 5 points in each national park, this would give
# us 31 * 5 = 155 points an therefore 11935 possible (distinct) connections
length(prot) * 5
nrow(combinations(155, 2))

# Plot everything
plot(prot)
points(points)

# For each point we want to know to which national park it belongs. To achieve
# this we convert the SpatialPoints object to a SpatialPointsDataFrame and
# assigne the data from the protected areas layer as new data
points <- SpatialPointsDataFrame(points, data = over(points, prot))

############################################################
#### Define Start and Endpoints: Method II
############################################################
# Load protected areas
prot <- "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(1Class).shp" %>%
  shapefile()

# Create points spaced out by a desired distance. We use a raster with desired
# grid-size for this purpose
points <- permeability %>%
  aggregate(., fact = (100 * 1000)/250, fun = mean) %>%
  rasterToPoints(., spatial = TRUE)

# Keep only those points that lie within protected areas
index <- gContains(prot, points, byid = TRUE) %>% rowSums(.)
points <- points[index > 0, ]

# Plot the points
plot(prot)
plot(points, add = TRUE, pch = 1, col = "red")

# Let's check the number of combinations we need to run
nrow(combinations(nrow(points), 2))

############################################################
#### Define Start and Endpoints: Method III
############################################################
# Load protected areas
prot <- "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(1Class).shp" %>%
  shapefile()

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

# Calculate area of each of the polygons
prot$Area <- prot %>%

  # Need to transform to utm
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

# Create points spaced out by a desired distance. We use a raster with desired
# grid-size for this purpose
points <- permeability %>%
  aggregate(., fact = (100 * 1000) / 250, fun = mean) %>%
  rasterToPoints(., spatial = TRUE)

# Keep only those points that lie within protected areas
index <- gContains(prot, points, byid = TRUE) %>% rowSums(.)
points <- points[index > 0, ]

# Find all polygons in which there is no point yet
missing <- prot %>%

  # Check intersecting points and polygons
  gIntersects(., points, byid = TRUE) %>%

  # Identify number of points in each polygon
  colSums(.) %>%

  # Convert result to vector
  as.vector(.)

# Identify those polygons with zero points
missing <- missing == 0

# Create centroids for these polygons, but make sure they are WITHIN the polygon
centroids <- gCentroidWithin(prot[missing, ])
# centroids <- gCentroid(prot[missing, ], byid = T)

# Put the centroids and the other points together
points <- rbind(SpatialPoints(points), SpatialPoints(centroids))
crs(points) <- crs(prot)

# Assign unique IDs to each point
points$ID <- 1:length(points)
points$ID <- as.numeric(points$ID)

# Also indicate to which protected area they belong to
together <- gContains(prot, points, byid = T)
points$ProtectedArea <- apply(together, 1, which)
points$ProtectedArea <- as.numeric(points$ProtectedArea)

# Plot to verify that IDs match protected area
plot(prot, border = rainbow(length(unique(prot$ID)))[prot$ID])
plot(points, col = rainbow(length(unique(prot$ID)))[points$ProtectedArea], add = T)

# Let's check the number of combinations we need to run
nrow(combinations(nrow(points), 2))

# Store the source points and areas
writeOGR(prot
  , dsn       = "03_Data/03_Results"
  , layer     = "99_SourceAreas"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
writeOGR(points
  , dsn       = "03_Data/03_Results"
  , layer     = "99_SourcePoints"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

############################################################
#### Least Cost Paths
############################################################
# We want to calculate least cost paths between all points. Unfortunately, the
# package gdistance can't handle to calculate all possible connections in one
# step. Let's find all possible pairwise connections ourselves. Note that we
# only look at distinct connections (i.e. only at 1->2 and not 2->1)
combis <- nrow(points) %>% combinations(., 2) %>%

  # Convert the result to a dataframe
  as.data.frame(.) %>%

  # Rename the columns more nicely
  rename(., Origin = V1, Destin = V2)

# Now find the shortest paths for all possible connections. I follow a slightly
# weird procedure since I try to avoid looping over all origins separately and
# rather cluster them (that is, I calculate from point 1 to all other points in
# one single step). This strictly increases performance
shortest <- suppressMessages(
  lapply(1:length(unique(combis$Origin)), function(x){
    shortestPath(trans
      , origin  = points[x, ]
      , goal    = points[combis$Destin[combis$Origin == x], ]
      , output  = "SpatialLines"
    )
  })
)

# Bind all resulting paths together
shortest <- do.call(rbind, shortest)

# Plot the least cost paths on top of the permeability map
plot(permeability, col = viridis(50))
lines(shortest, col = "white")
points(points, col = "white")

# Make the object a SpatialLinesDataFrame and add some artificial data
shortest <- as(shortest, "SpatialLinesDataFrame")
shortest@data <- data.frame(id = 1:length(shortest))

# Write the least cost paths to file
writeOGR(shortest
  , "03_Data/03_Results"
  , "99_LeastCostPaths"
  , driver = "ESRI Shapefile"
  , overwrite = TRUE
)

# Read the paths back in
shortest <- readOGR("03_Data/03_Results/99_LeastCostPaths.shp")

# Let's bufferthe LCPs
shortest_buff <- gBuffer(shortest
  , byid  = TRUE
  , width = metersToDegrees(2650)
)

# Rasterize them
tic()
shortest_rast <- rasterizeTerra(
    vect(shortest_buff)
  , rast(permeability)
)
toc()

# Plot
plot(shortest_rast, col = viridis(50))

# Store the rasterized LCPs
writeRaster(
    shortest_rast
  , "03_Data/03_Results/99_LeastCostPaths.tif"
  , overwrite = TRUE
)

############################################################
#### Least Cost Corridors
############################################################
# Create a list of points
points_l <- list()
for (i in 1:length(points)){
  points_l[[i]] <- points[i, ]
}

# Keep track of the processing time
tic()

# Run a function to calculate the cumulative cost map for each point
cost <- lapply(points_l, function(x){

  # Calculate the costmap for the respective point
  costi <- accCost(trans, x)

  # Store the layer to a temporary file to avoid memory overflow
  costi <- writeRaster(costi, tempfile())

  # Clean garbage (althouth R should do that automatically)
  gc()

  # Return the costmap
  return(costi)
})

# Load all costmaps into a single stack
cost <- stack(cost)

# Write the stack to a temporary file
cost <- stackSave(cost, tempfile())

# Prepare an empty raster onto which we will add all least cost corridors
r0 <- cost[[1]]
values(r0) <- 0

# Loop through all combinations and calculate a least cost corridor
for (i in 1:nrow(combis)){

  # Get the cost maps of the two points
  cost1 <- cost[[combis[i, "Origin"]]]
  cost2 <- cost[[combis[i, "Destin"]]]

  # Sum the layers
  corr <- calc(stack(cost1, cost2), sum)

  # Calculate threshold below which the costs are not higher than + x% of the
  # minimal costs
  # threshold <- minValue(corr) * 1.05
  threshold <- minValue(corr) * 1.05

  # Make values above the threshold NAs
  corr <- reclassify(corr, c(threshold, Inf, NA))

  # Scale the costs to 0 and 1
  min <- minValue(corr)
  max <- maxValue(corr)
  corr <- calc(corr, function(x){(x - min) / (max - min)})

  # assign all NAs a value of 1
  corr <- reclassify(corr, c(NA, NA, 1), right = FALSE)

  # Add the corridor to the r0 raster
  r0 <- calc(stack(r0, corr), sum)

  # Collect garbage
  gc()

}

# Stop the timer
toc()

# Rescale values to between 0 and 1 and invert them
min <- minValue(r0)
max <- maxValue(r0)
r0 <- calc(r0, function(x){(x - min) / (max - min)})
r0 <- r0 * (-1) + 1

# Plot the result
plot(r0, col = viridis(50))

# Write result to file
writeRaster(r0
  , "03_Data/03_Results/99_LeastCostCorridors.tif"
  , overwrite = TRUE
)
