################################################################################
#### Figures for the Graphical Abstract
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse) # To wrangle data
library(raster)    # To handle spatial data
library(terra)     # To handle spatial data
library(rgeos)     # To manipulate spatial data
library(igraph)    # To create networks
library(NLMR)      # To create random landscapes
library(gdistance) # To calculate least-cost paths
library(pbmcapply) # For multicore abilities with progress bar
library(spatstat)  # To quickly rasterize lines
library(maptools)  # To quickly rasterize lines
library(viridis)   # For nice colors

################################################################################
#### Simulate Permeability Surface (could be any set of covariates)
################################################################################
# Create a permeability surface
set.seed(1)
perm <- nlm_gaussianfield(ncol = 100, nrow = 100)
names(perm) <- "permeability"

# Visualize
plot(perm)

# Function to extend a raster artificially
extendRaster <- function(x, y){
  na_mask <- is.na(x)
  vals <- values(x)
  vals <- na.omit(vals)
  r <- extend(x, y)
  na_mask <- extend(na_mask, y, value = 0)
  indices <- which(is.na(values(r)))
  values(r)[indices] <- vals[runif(length(indices), min = 1, max = length(vals))]
  r <- mask(r, na_mask, maskvalue = 1, updatevalue = NA)
  return(r)
}

# Expand permeability surface artificially (so that dispersers can move outside
# too)
perm <- extendRaster(perm, extent(perm) + c(-10, 10, -10, 10))

# Distribute 5 source areas in suitable habitat patches
source_areas <- SpatialPoints(rbind(
    c(40, 95)
  , c(85, 82)
  , c(20, 60)
  , c(58, 62)
  , c(65, 40)
  , c(56, 05)
))
source_areas$ID <- 1:length(source_areas)

# Create buffers
source_areas <- gBuffer(source_areas, byid = T, width = 4)

# Visualize them
plot(perm)
plot(source_areas, add = T)
text(source_areas, "ID")

# Currently, our permeability surface is completely random and does not contain
# natural "corridors". However, to better illustrate the point, we'd like to
# have a permeability surface that "channels" our virtual dispersers, such that
# we can get some nice maps afterwards. For this, I'll first generate a set of
# least cost corridors below which I will increase permeability artificially.
# Hence, let's first define connections between which we want to have such
# corridors.
conns <- rbind(
    c(1, 2)
  , c(1, 3)
  , c(1, 4)
  , c(3, 5)
  , c(3, 6)
  , c(4, 2)
  , c(4, 5)
  , c(5, 6)
)

# Create transition matrix
t <- transition(perm ** 2, directions = 8, transitionFunction = mean)
t <- geoCorrection(t, type = "c")

# Calculate least-cost paths for those connections
paths <- list()
for (i in 1:nrow(conns)){
  paths[[i]] <- shortestPath(t
    , origin  = gCentroid(source_areas[conns[i, 1], ])
    , goal    = gCentroid(source_areas[conns[i, 2], ])
    , output  = "SpatialLines"
  )
}
paths <- do.call(rbind, paths)

# Let's buffer the paths so that we get some nice corridors below which we can
# increase permeability
paths_buff <- gBuffer(paths, width = 3, byid = F)

# Visualize them
plot(perm)
plot(paths_buff, add = T)
plot(source_areas, add = T)
text(source_areas, "ID")

# Make those areas more permeable. We will also make all source areas are more
# permeable
add1 <- mask(perm, paths_buff)
add1[is.na(add1)] <- 0
add2 <- mask(perm, source_areas)
add2[is.na(add2)] <- 0
plot(perm + 2 * add1 + 0.3 * add2)
perm <- perm + 1 * add1 + 0.3 * add2

# Visualize final permeability map
plot(perm)
plot(source_areas, add = T)
plot(paths_buff, add = T)
text(source_areas, "ID")

################################################################################
#### Simulate Movement
################################################################################
# Function to simulate movement
move <- function(xy
    , covars   = NULL    # Covariate layer
    , prefs    = NULL    # Preferences towards covariate layer
    , sl_dist  = NULL    # Step length distribution
    , n_steps  = 10      # Number of simulated steps
    , n_rsteps = 25      # Number of proposed random steps
    , stop     = TRUE    # Should the algorithm stop if a boundary is hit?
  ){

  # Create a new dataframe based on the source point. Note that we draw random
  # turning angles to start off
  track <- data.frame(
      x       = coordinates(xy)[, 1]
    , y       = coordinates(xy)[, 2]
    , absta_  = runif(1, min = 0, max = 2 * pi)
    , ta_     = runif(1, min = -pi, max = pi)
    , sl_     = NA
  )

  # Simulate random steps
  for (i in 1:n_steps){

    # Prepare an empty list in which we can store the random steps
    rand <- list()

    # Draw random turning angles
    ta_new <- runif(n_rsteps
      , min = -pi
      , max = +pi
    )

    # Draw random step lengths
    sl_new <- rgamma(n_rsteps
      , shape = sl_dist$shape
      , scale = sl_dist$scale
    )

    # Make sure that the steps cover at least a minimal distance
    sl_new[sl_new < 0.001] <- 0.001

    # Put the step lengths and turning angles into a new dataframe. These are
    # our proposed random steps.
    rand <- data.frame(
        absta_  = track$absta_[i] + ta_new
      , ta_     = ta_new
      , sl_     = sl_new
    )

    # We need to make sure that the absolute turning angle ranges from 0 to 2 *
    # pi
    rand$absta_[rand$absta_ > 2 * pi] <-
      rand$absta_[rand$absta_ > 2 * pi] - 2 * pi
    rand$absta_[rand$absta_ < 0] <-
      rand$absta_[rand$absta_ < 0] + 2 * pi

    # Calculate new endpoints
    rand$x <- track$x[i] + sin(rand$absta_) * rand$sl_
    rand$y <- track$y[i] + cos(rand$absta_) * rand$sl_

    # Create spatial points from endpoints
    coordinates(rand) <- c("x", "y")

    # Depending on the answer in the beginning, the loop breaks if one of the
    # new coordinates is outside the map boundaries
    if (stop){

      # If one of the random steps leaves the extent, we break the loop
      extent  <- as(extent(covars), "SpatialPolygons")
      inside  <- as.vector(gContainsProperly(extent, rand, byid = TRUE))
      if (sum(!inside) > 0) break

    } else {

      # Remove any point/step that does not fully lie withing the boundaries of
      # our map (otherwise we can't extract the covariates below). One could
      # also resample those steps until they fully lie within the study area.
      extent  <- as(extent(covars), "SpatialPolygons")
      keep    <- as.vector(gContainsProperly(extent, rand, byid = TRUE))
      rand    <- rand[keep, ]

    }

    # Coerce back to regular dataframe
    rand <- as.data.frame(rand, xy = T)

    # Prepare a "line" for each random step. We first need the coordinates of
    # the steps for this
    begincoords <- track[i, c("x", "y")] # Start point
    endcoords   <- rand[, c("x", "y")]  # Possible endpoints

    # Create spatial lines
    steps <- list()
    for (j in 1:nrow(endcoords)){
      steps[[j]] <- spLines(SpatialPoints(rbind(begincoords, endcoords[j, ])))
    }
    steps <- do.call(rbind, steps)

    # Extract covariates along steps
    extracted <- terra::extract(rast(covars), vect(steps), fun = mean)[, 2]

    # Bind with other data
    rand <- cbind(rand, extracted)
    names(rand)[ncol(rand)] <- "permeability"

    # Calculate cos_ta and log_sl
    rand$cos_ta <- cos(rand$ta_)
    rand$log_sl <- log(rand$sl_)

    # Prepare model matrix
    mat <- model.matrix(~ permeability + sl_ + log_sl + cos_ta, rand)
    mat <- mat[ , 2:ncol(mat)]

    # Calculate selection scores
    score <- exp(mat %*% prefs)

    # Convert scores to probabilities
    probs <- score / sum(score)

    # Keep only the step with the highest score
    rand <- rand[sample(nrow(rand), 1, prob = probs), ]

    # Add the step to our track
    track <- rbind(
        track[, c("x", "y", "absta_", "ta_", "sl_")]
      , rand[, c("x", "y", "absta_", "ta_", "sl_")]
    )
  }
  return(track)
}

# Specify movement parameters
sl_dist <- list(
    scale = 3
  , shape = 0.5
)
preferences <- c(
    Permeability = 5
  , sl_          = 0
  , log_sl       = 0
  , cos_ta       = 5
)
n_disp <- 1000
n_steps <- 150
n_rsteps <- 25

# Initiate 100 dispersers in each source area
source_points <- lapply(1:length(source_areas), function(x){
  source_point <- spsample(source_areas[x, ], n = n_disp, type = "random")
  source_point$SourceArea <- rep(source_areas$ID[x], n_disp)
  return(source_point)
}) %>% do.call(rbind, .)

# Visualize the source points
plot(perm)
plot(paths_buff, add = T)
plot(source_points, col = source_points$SourceArea, add = T, pch = 16, cex = 0.3)
plot(source_areas, add = T, border = "black", lwd = 2)
text(source_areas, "ID", col = "white")

# Initiate a disperser at each source point
sims <- pbmclapply(
    X                  = 1:length(source_points)
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(x){
    sim <- move(
        xy       = source_points[x, ]
      , covars   = perm
      , prefs    = preferences
      , n_steps  = n_steps
      , n_rsteps = n_rsteps
      , stop     = F
      , sl_dist  = sl_dist
    )
    rownames(sim) <- NULL
    sim <- sim[, c("x", "y")]
    sim$ID <- x
    sim$SourceArea <- source_points$SourceArea[x]
    sim$StepNumber <- 1:nrow(sim)
    return(sim)
})

# Store the simulations to file
write_rds(sims, "test_sims.rds")

# Coerce simulated tracks into spatial lines
tracks <- lapply(1:length(sims), function(x){
  track <- sims[[x]]
  coordinates(track) <- c("x", "y")
  track <- spLines(track)
  track$ID <- sims[[x]]$ID[1]
  track$SourceArea <- sims[[x]]$SourceArea[1]
  return(track)
}) %>% do.call(rbind, .)

# Visualize them
plot(perm)
plot(tracks, add = T)

################################################################################
#### Heatmap
################################################################################
# Function to rasterize spatial lines (l) quickly onto a raster (r)
rasterizeSpatstat <- function(l, r){

  # Create im layer
  values(r) <- 0
  im <- as.im.RasterLayer(r)
  summed <- im

  # Prepare progress bar
  pb <- txtProgressBar(
      min     = 0
    , max     = length(l)
    , initial = 0
    , style   = 3
    , width   = 55
  )

  # Go through all lines and rasterize them
  for (y in 1:length(l)){
    line    <- as.psp(l[y, ], window = im)
    line    <- as.mask.psp(line)
    line_r  <- as.im.owin(line, na.replace = 0)
    summed  <- Reduce("+", list(summed, line_r))
    setTxtProgressBar(pb, y)
  }

  # Return heatmap as a raster
  return(raster(summed))

}

# Rasterize simulated trajctories
heatmap <- rasterizeSpatstat(
    l = as(tracks, "SpatialLines")
  , r = disaggregate(perm, fact = 2)
)

# Visualize them
plot(heatmap, col = viridis(50))

################################################################################
#### Betweenness
################################################################################
# Overlay the study area with a regular raster grid
regular <- raster(extent(perm), ncol = 50, nrow = 50, vals = 1:50 ** 2)

# Bind all simulations into a single dataframe
sims <- do.call(rbind, sims)

# We want to determine all connections from one raster cell on the regular grid
# to another cell as per our simulations. For this, we'll identify the
# "visitation history" of each trajectory, i.e. the sequence of raster cells
# across which it moves.
visits <- data.frame(
    ID         = sims$ID
  , StepNumber = sims$StepNumber
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , cell       = raster::extract(regular, sims[, c("x", "y")])
)

# Function to retrieve the visitation history from a sequence of values
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

# Try it
visitHist(c(1, 2, 2, 2, 3, 3, 1))
visitHist(c(1, 2, 2, 2, 3, 3, 1), singlecount = T)

# Apply it to each simulated trajectory separately
visits <- nest(visits, data = -ID)
history <- pbmclapply(
    X                  = visits$data
  , ignore.interactive = T
  , mc.cores           = 1
  , FUN                = function(y){
      visitHist(y$cell, singlecount = T)
  }) %>%
  do.call(rbind, .) %>%
  group_by(from, to) %>%
  summarize(TotalConnections = sum(TotalConnections), .groups = "drop") %>%
  ungroup() %>%
  mutate(weight = mean(TotalConnections) / TotalConnections)

# Function to calculate network metrics
netMet <- function(
    network   = NULL    # The network based on which the metrics are calculated
  , raster    = NULL    # The raster onto which the metrics are calculated
  , tempfile  = F       # Should the resulting raster go to a temporary file?
  , metrics = c("betweenness", "closeness", "degree")
  ){

    # Calculate the desired network metrics and put them as raster values
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

# Coerce visits into a network (each cell in the regular raster serves as
# potential node/vertex)
vertices <- 1:ncell(regular)
net <- graph_from_data_frame(history, vertices = vertices)

# Calculate Betweenness
betweenness <- netMet(
    network  = net
  , raster   = regular
  , metrics  = "betweenness"
  , tempfile = T
)

# Visualize it
plot(betweenness, col = magma(50))

################################################################################
#### Interpatch Connectivity
################################################################################
# Determine for each trajectory into which areas it moves
ints <- gIntersects(source_areas, tracks, byid = T)
visits <- cbind(tracks$SourceArea, ints) %>%
  as.data.frame() %>%
  gather(key = To, value = Reached, 2:ncol(.)) %>%
  setNames(c("From", "To", "Reached")) %>%
  subset(From != To & Reached == 1) %>%
  group_by(From, To) %>%
  count() %>%
  mutate(n = n / n_disp)

# Create network
net <- graph_from_data_frame(visits, vertices = source_areas$ID)
lay <- coordinates(gCentroid(source_areas, byid = T))

################################################################################
#### Visualize All Maps
################################################################################
# Visualize
par(mfrow = c(2, 2))
plot(perm)
plot(paths_buff, add = T)
plot(source_areas, add = T)
text(source_areas, "ID")
plot(heatmap, col = viridis(50))
plot(betweenness, col = magma(50))
plot(perm, col = "white")
plot(net
  , layout      = lay
  , rescale     = F
  , add         = T
  , vertex.size = 500
  , edge.label  = round(E(net)$n, 2)
  , edge.width  = E(net)$n * 20
  , edge.curved = 0.2
  , edge.arrow.size = 0.5
)
