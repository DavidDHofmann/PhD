############################################################
#### Comparing Least-Cost to Simulated Paths
############################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load packages
library(tidyverse)    # For data wrangling
library(raster)       # For raster manipulation
library(terra)        # For raster manipulation
library(viridis)      # For nicer colors
library(rgeos)        # To manipulate spatial data
library(gdistance)    # To calculate least-cost paths
library(tmap)         # To prepare nice maps
library(glmmTMB)      # Modelling
library(parallel)     # To use multiple cores

# Load custom functions
source("Functions.r")

############################################################
#### Prepare Dispersal Data
############################################################
# Load dispersal trajectories of individuals that we did not use to train our
# iSSF models and coerce their relocations to spatial lines
disp <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv() %>%
  subset(
    State == "Disperser" &
    DogName %in% c("Pula", "Rattan", "Appalachia", "Karisimbi")
  ) %>%
  group_by(DogName) %>%
  nest() %>%
  mutate(Tracks = map(data, function(x){
    coordinates(x) <- c("x", "y")
    x <- createSegments(x)
    crs(x) <- CRS("+init=epsg:4326")
    return(x)
  }))

# Identify origin, destination, and furthest relocation of each track. We also
# want to know to know the number of realized steps and the total time (in days)
# travelled
disp <- disp %>%

  # Identification of origin
  mutate(Origin = map(data, function(x){
    reloc <- head(x, 1)
    coordinates(reloc) <- c("x", "y")
    crs(reloc) <- CRS("+init=epsg:4326")
    return(reloc)
  })) %>%

  # Identification of destination
  mutate(Destin = map(data, function(x){
    reloc <- tail(x, 1)
    coordinates(reloc) <- c("x", "y")
    crs(reloc) <- CRS("+init=epsg:4326")
    return(reloc)
  })) %>%

  # Calculate buffer around destination
  mutate(DestinBuffer = map(Destin, function(x){
    gBuffer(x, width = metersToDegrees(2650))
  })) %>%

  # Identification of furthest point
  mutate(Furthest = map(data, function(x){
    origin <- head(x, 1)
    distances <- sqrt((x$x - origin$x) ** 2 + (x$y - origin$y) ** 2)
    index <- which(distances == max(distances))[1]
    furthest <- x[index, ]
    coordinates(furthest) <- c("x", "y")
    crs(furthest) <- CRS("+init=epsg:4326")
    return(furthest)
  })) %>%

  # Calculate buffer around furthest point
  mutate(FurthestBuffer = map(Furthest, function(x){
    gBuffer(x, width = metersToDegrees(2650))
  })) %>%

  # Calculate number of steps
  mutate(NrSteps = map(Tracks, function(x){
    length(x)
  }) %>% do.call(rbind, .)) %>%

  # Calculate total dispersal time
  mutate(DaysTravelled = map(data, function(x){
    max(x$Timestamp) - min(x$Timestamp)
  }) %>% do.call(rbind, .))

# Show tibble
print(disp)

# Visualize all
par(mfrow = c(2, 2))
  for (n in 1:4){
    plot(disp$Tracks[[n]], main = disp$DogName[[n]])
    plot(disp$Origin[[n]], add = T, col = "red", pch = 20)
    plot(disp$Destin[[n]], add = T, col = "green", pch = 20)
    plot(disp$DestinBuffer[[n]], add = T, border = "green", pch = 20)
    plot(disp$Furthest[[n]], add = T, col = "purple", pch = 20)
    plot(disp$FurthestBuffer[[n]], add = T, border = "purple", pch = 20)
  }

############################################################
#### Prepare Permeability Surfaces
############################################################
# Read in the model results from the analysis of the 4-hourly fixes
models <- readRDS("03_Data/03_Results/99_ModelSelection.rds")

# Subset to the best model and identify the model estimates
coeffs <- models$Model[[1]] %>% getCoeffs()

# We want to use the derived model to predict a resistance surface.
# Let's create a simplified dataframe of the model coefficients
pred <- as.data.frame(dplyr::select(coeffs, Coefficient))

# For easier indexing we assign the covariates as row-names
rownames(pred) <- coeffs$Covariate

# Load required covariate layers
water <- stack("03_Data/02_CleanData/01_LandCover_Water(Merged).tif")
trees <- stack("03_Data/02_CleanData/01_LandCover_TreeCover_MODIS.tif")
shrubs <- stack("03_Data/02_CleanData/01_LandCover_NonTreeVegetation_MODIS.tif")
humans <- raster("03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Buffer5000).tif")

# Identify the mid-date of each dispersal event. We will use it to figure out
# which dynamic maps we should use to prepare the permeability surface
dates_disp <- lapply(disp$data, function(x){
  as.Date(mean(range(x$Timestamp)))
})

# Extract dates of all watermaps (dates are the same for water, trees and shrubs)
dates <- names(water) %>%
  substr(start = 2, stop = nchar(.)) %>%
  as.Date(format = "%Y.%m.%d")

# For each dispersal event we now identify the closest watermap
water <- lapply(dates_disp, function(x){
  water[[which.min(abs(x - dates))]]
}) %>% stack()

# For each dispersal event we now identify the closest treemap
trees <- lapply(dates_disp, function(x){
  trees[[which.min(abs(x - dates))]]
}) %>% stack()

# For each dispersal event we now identify the closest shrubmap
shrubs <- lapply(dates_disp, function(x){
  shrubs[[which.min(abs(x - dates))]]
}) %>% stack()

# Visualize the watermaps and corresponding trajectories
par(mfrow = c(2, 2))
  for (n in 1:4){
    plot(crop(water[[n]], disp$Tracks[[n]])
      , col = c("white", "blue")
      , main = disp$DogName[[n]])
    plot(disp$Tracks[[n]], add = T)
  }

# Visualize the treemaps and corresponding trajectories
par(mfrow = c(2, 2))
  for (n in 1:4){
    plot(crop(trees[[n]], disp$Tracks[[n]])
      , main = disp$DogName[[n]])
    plot(disp$Tracks[[n]], add = T)
  }

# Visualize the shrubmaps and corresponding trajectories
par(mfrow = c(2, 2))
  for (n in 1:4){
    plot(crop(shrubs[[n]], disp$Tracks[[n]])
      , main = disp$DogName[[n]])
    plot(disp$Tracks[[n]], add = T)
  }

# Define the extents of the analysis. We will first use a larger extent to crop
# the water layer as this will allow us to calculate "DistanceToWater" more
# reliably
extent_1 <- extent(do.call(rbind, disp$Tracks)) + c(-1, 1, -1, 1)
extent_2 <- extent(do.call(rbind, disp$Tracks)) + c(-0.5, 0.5, -0.5, 0.5)

# Crop water layers
water <- crop(water, extent_1)

# Calculate distace to water
distancetowater <- suppressMessages(
  mclapply(1:nlayers(water), mc.cores = detectCores() - 1, function(x){
    distanceTo(water[[x]], value = 1)
  }) %>% stack()
)

# Crop all layers to the small extent
water <- crop(water, extent_2)
distancetowater <- crop(distancetowater, extent_2)
trees <- crop(trees, extent_2)
shrubs <- crop(shrubs, extent_2)
humans <- crop(humans, extent_2)

# Specify layernames
names <- c("Water", "DistanceToWater", "Trees", "Shrubs", "HumansBuff5000")

# We also need to scale the covariates using the same scaling parameters that we
# used to scale the covariates during the modelling part
scaling <- read_rds("03_Data/03_Results/99_ScalingSSF.rds")

# Make some nice rownames
rownames(scaling) <- scaling$ColumnName

# Write a function to scale a layer using the table above
scaleLayer <- function(layer, table){
  scale(layer,
      center  = table[names(layer), ]$Center
    , scale   = table[names(layer), ]$Scale
  )
}

# We will also need layers for the turning angles and step lengths. Actually,
# they simply inflate the probabilities. For the step length layer we can use
# the average step length (around 2650)
`log(sl_)` <- water[[1]]
values(`log(sl_)`) <- log(2650)
names(`log(sl_)`) <- "log(sl_)"

# For the turning angles we will use zero
`cos(ta_)` <- water[[1]]
values(`cos(ta_)`) <- cos(0)
names(`cos(ta_)`) <- "cos(ta_)"

# Put covariates into the tibble and scale them
disp$Covars <- lapply(1:nrow(disp), function(x){
  covars <- list(water[[x]], sqrt(distancetowater[[x]]), trees[[x]], shrubs[[x]], humans)
  covars <- stack(covars)
  names(covars) <- names
  for (i in 1:nlayers(covars)){
    covars[[i]] <- scaleLayer(covars[[i]], scaling)
  }
  covars <- stack(`cos(ta_)`, `log(sl_)`, covars)
  return(covars)
})

# Visualize
plot(disp$Covars[[1]])
plot(disp$Covars[[2]])
plot(disp$Covars[[3]])
plot(disp$Covars[[4]])

# Prepare permeability suraface for each individual
disp <- disp %>%
  mutate(Permeability = map(Covars, function(x){
    permeability <- exp(sum(
        pred["cos(ta_)", ]        * x[["cos.ta_."]]
      , pred["log(sl_)", ]        * x[["log.sl_."]]
      , pred["Water", ]           * x[["Water"]]
      , pred["DistanceToWater", ] * x[["DistanceToWater"]]
      , pred["Shrubs", ]          * x[["Shrubs"]]
      , pred["HumansBuff5000", ]  * x[["HumansBuff5000"]]
      , pred["Trees", ]           * x[["Trees"]]
    ))
    upper <- quantile(values(permeability), 0.99, na.rm = TRUE)
    lower <- quantile(values(permeability), 0.01, na.rm = TRUE)
    values(permeability)[values(permeability) > upper] <- upper
    values(permeability)[values(permeability) < lower] <- lower
    return(permeability)
  }))

# Visualize all permeability surfaces
par(mfrow = c(2, 2))
  for (n in 1:4){
    plot(disp$Permeability[[n]], col = viridis(50), main = disp$DogName[[n]])
    plot(disp$Tracks[[n]], add = T, main = disp$DogName[[n]], col = "white")
    plot(disp$Origin[[n]], add = T, col = "red")
    plot(disp$Destin[[n]], add = T, col = "green")
    plot(disp$Furthest[[n]], add = T, col = "purple")
  }

# ############################################################
# #### Least-Cost Paths
# ############################################################
# # Calculate transition layer
# disp$Transition <- suppressMessages(
#   mclapply(1:nrow(disp), mc.cores = detectCores() - 1, function(x){
#     trans <- transition(
#         disp$Permeability[[x]]
#       , directions = 8
#       , transitionFunction = mean
#     )
#     trans <- geoCorrection(trans, type = "c", multpl = FALSE)
#     gc()
#     return(trans)
#   })
# )
#
# # Calculate Least-Cost Paths to destination of dispersal
# disp$LeastCostPathsLast <- lapply(1:nrow(disp), function(x){
#   path <- shortestPath(disp$Transition[[x]]
#     , origin  = disp$Origin[[x]]
#     , goal    = disp$Destin[[x]]
#     , output  = "SpatialLines"
#   )
#   return(path)
# })
#
# # Calculate Least-Cost Paths to furthest point
# disp$LeastCostPathsFurthest <- lapply(1:nrow(disp), function(x){
#   path <- shortestPath(disp$Transition[[x]]
#     , origin  = disp$Origin[[x]]
#     , goal    = disp$Furthest[[x]]
#     , output  = "SpatialLines"
#   )
#   return(path)
# })
#
# # Visualize LCPs to destination
# par(mfrow = c(2, 2))
# for (n in 1:4){
#   plot(disp$Permeability[[n]], col = viridis(50), main = disp$DogName[[n]])
#   plot(disp$Tracks[[n]], add = T, main = disp$DogName[[n]], col = "white")
#   plot(disp$Origin[[n]], add = T, col = "red")
#   plot(disp$Destin[[n]], add = T, col = "green")
#   plot(disp$Furthest[[n]], add = T, col = "purple")
#   plot(disp$LeastCostPathsLast[[n]], add = T, col = "red")
# }
#
# # Visualize LCPs to furthest point
# par(mfrow = c(2, 2))
# for (n in 1:4){
#   plot(disp$Permeability[[n]], col = viridis(50), main = disp$DogName[[n]])
#   plot(disp$Tracks[[n]], add = T, main = disp$DogName[[n]], col = "white")
#   plot(disp$Origin[[n]], add = T, col = "red")
#   plot(disp$Destin[[n]], add = T, col = "green")
#   plot(disp$Furthest[[n]], add = T, col = "purple")
#   plot(disp$LeastCostPathsFurthest[[n]], add = T, col = "red")
# }
#
# ############################################################
# #### Least-Cost Corridors
# ############################################################
# # Calculate Least-Cost Corridor to destination
# disp$LeastCostCorridorsLast <- suppressMessages(
#   mclapply(1:nrow(disp), mc.cores = detectCores() - 1, function(x){
#     cost1 <- accCost(disp$Transition[[x]], disp$Origin[[x]])
#     cost2 <- accCost(disp$Transition[[x]], disp$Destin[[x]])
#     corr <- calc(stack(cost1, cost2), sum)
#     threshold <- minValue(corr) * 1.50
#     corr <- reclassify(corr, c(threshold, Inf, NA))
#     corr <- (corr - minValue(corr)) / (maxValue(corr) - minValue(corr))
#     return(corr)
#   })
# )
#
# # Calculate Least-Cost Corridor to furthest point
# disp$LeastCostCorridorsFurthest <- suppressMessages(
#   mclapply(1:nrow(disp), mc.cores = detectCores() - 1, function(x){
#     cost1 <- accCost(disp$Transition[[x]], disp$Origin[[x]])
#     cost2 <- accCost(disp$Transition[[x]], disp$Furthest[[x]])
#     corr <- calc(stack(cost1, cost2), sum)
#     threshold <- minValue(corr) * 1.50
#     corr <- reclassify(corr, c(threshold, Inf, NA))
#     corr <- (corr - minValue(corr)) / (maxValue(corr) - minValue(corr))
#     return(corr)
#   })
# )
#
# # Visualize LCCs to destination
# par(mfrow = c(2, 2))
# for (n in 1:4){
#   plot(disp$LeastCostCorridorsLast[[n]], col = rev(magma(50)))
#   plot(disp$Tracks[[n]], add = T, main = disp$DogName[[n]], col = "black")
#   plot(disp$Origin[[n]], add = T, col = "red")
#   plot(disp$Destin[[n]], add = T, col = "green")
#   plot(disp$Furthest[[n]], add = T, col = "purple")
# }
#
# # Visualize LCCs to furthest point
# par(mfrow = c(2, 2))
# for (n in 1:4){
#   plot(disp$LeastCostCorridorsFurthest[[n]], col = rev(magma(50)))
#   plot(disp$Tracks[[n]], add = T, main = disp$DogName[[n]], col = "black")
#   plot(disp$Origin[[n]], add = T, col = "red")
#   plot(disp$Destin[[n]], add = T, col = "green")
#   plot(disp$Furthest[[n]], add = T, col = "purple")
# }

############################################################
#### Simulations
############################################################
# Load the movement model
move.mod <- "03_Data/03_Results/99_MovementModel.rds" %>% readRDS()

# Select the most parsimonious model
move.mod <- move.mod$Model[[1]]

# Get the coefficients
coeffs <- getCoeffs(move.mod)

# We want to use the derived model to predict a selection score. Let's convert
# the coefficients into a simpler dataframe
pred <- dplyr::select(coeffs, Coefficient) %>% as.matrix() %>% as.vector()

# For easier indexing we assign the covariates as row-names
names(pred) <- coeffs$Covariate

# We will need to simulate random step lengths. For this we can use the
# gamma distribution that we fitted to create random steps earlier.
gamma <- "03_Data/03_Results/99_GammaDistribution.rds" %>% readRDS()


# Run the simulation multiple times at each start location
disp$Simulations <- lapply(1:nrow(disp), function(x, times = 250){

  # Duplicate the start point multiple times
  point <- disp$Origin[[x]]
  point <- point[rep(1, times), ]

  # Simulate dispersal for each point
  sims <- suppressMessages(
    mclapply(1:length(point), mc.cores = detectCores() - 1, function(y){

      # Run simulation
      sim <- disperse(
          x                   = point[y, ]
        , y                   = disp$Covars[[x]]
        , z                   = pred
        , w                   = gamma
        , n_rsteps            = 25
        , n_steps             = 200
        , stop                = T
      )

      # Assign ID that corresponds to the iteration number
      sim$ID <- y

      # Return the simulation
      return(sim)
    })
  )

  # Bind all simulations
  sims <- do.call(rbind, sims)

  # Return the dataframe
  return(sims)
})

# Backup the data we created so far
write_rds(disp, "03_Data/03_Results/99_LeastCostvsSimulations.rds")

############################################################
#### Simulate with Attraction
############################################################
# Run simulation
n <- 2
sim <- disperse(
    x                   = disp$Origin[[n]]
  , y                   = disp$Covars[[n]]
  , z                   = pred
  , w                   = gamma
  , AttractionPoint     = disp$Destin[[n]]
  , AttractionStrength  = -0.09
  , n_rsteps            = 25
  , n_steps             = 50
  , stop                = T
)
coordinates(sim) <- c("x", "y")
sim <- spLines(sim)

# Visualize
plot(disp$Permeability[[n]], col = viridis(50), main = disp$DogName[[n]])
plot(sim, add = T, col = "orange")
plot(disp$Tracks[[n]], add = T, main = disp$DogName[[n]], col = "white")
plot(disp$Origin[[n]], add = T, col = "red")
plot(disp$Destin[[n]], add = T, col = "green")
plot(disp$Furthest[[n]], add = T, col = "purple")

############################################################
#### Simulations Rasterization: All Trajectories
############################################################
# Coerce simulated trajectories to proper lines
disp <- disp %>% mutate(SimulationsLines = map(Simulations, function(x){

  # Nest the data so that each trajectory receives its own row
  sims <- x %>% group_by(ID) %>% nest()

  # Move from point to step represenation (create lines from points)
  sims_traj <- lapply(1:nrow(sims), function(x){
    p <- SpatialPoints(sims$data[[x]][, c("x", "y")])
    l <- spLines(p)
    return(l)
  }) %>% do.call(rbind, .)

  # Create unique ID
  sims_traj$ID <- 1:length(sims_traj)

  # Return the lines
  return(sims_traj)
}))

# Create raster to rasterize simulated tracks onto
r <- rast(raster(disp$Covars[[1]]))
values(r) <- 0

# Rasterize all simulated trajectories
disp$RasterizedSimulations <- lapply(1:nrow(disp), function(x){
  rasterized <- suppressMessages(
    rasterizeTerra(vect(disp$SimulationsLines[[x]]), r)
  )
  return(rasterized)
})

# Visualize them
par(mfrow = c(2, 2))
for (n in 1:4){
  plot(sqrt(disp$RasterizedSimulations[[n]]), col = viridis(50))
  plot(disp$Tracks[[n]], add = T, col = "white")
}

############################################################
#### Simulations Rasterization: Successfull Trajectories
############################################################
# Identify all simulated trajectories that reach the destination and crop the
# lines after the intersection
disp$SimulationsToDestin <- lapply(1:nrow(disp), function(x){
  inter <- gIntersects(disp$DestinBuffer[[x]], disp$SimulationsLines[[x]], byid = T)
  inter <- which(inter)
  sub <- disp$SimulationsLines[[x]][inter, ]
  if (nrow(sub) > 0){
    cut <- list()
    for (i in 1:nrow(sub)){
      cut[[i]] <- cutLineClose(
          line      = sub[i, ]
        , point     = disp$Destin[[x]]
        , polygon   = disp$DestinBuffer[[x]]
        , precision = 250
      )
    }
    sub <- do.call(rbind, cut)
  }
  return(sub)
})

# Identify all simulated trajectories that reach the furthest point and crop the
# lines after the intersection
disp$SimulationsToFurthest <- lapply(1:nrow(disp), function(x){
  inter <- gIntersects(disp$FurthestBuffer[[x]], disp$SimulationsLines[[x]], byid = T)
  inter <- which(inter)
  sub <- disp$SimulationsLines[[x]][inter, ]
  if (nrow(sub) > 0){
    cut <- list()
    for (i in 1:nrow(sub)){
      cut[[i]] <- cutLineClose(
          line      = sub[i, ]
        , point     = disp$Furthest[[x]]
        , polygon   = disp$FurthestBuffer[[x]]
        , precision = 250
      )
    }
    sub <- do.call(rbind, cut)
  }
  return(sub)
})

# Rasterize simulations that reached the destination
disp$RasterizedSimulationsToDestin <- lapply(1:nrow(disp), function(x){
  if (length(disp$SimulationsToDestin[[x]]) > 0){
    rasterized <- suppressMessages(
      rasterizeTerra(vect(disp$SimulationsToDestin[[x]]), r)
    )
  } else {
    rasterized <- raster(r)
    values(rasterized) <- 0
  }
  return(rasterized)
})

# Rasterize simulations that reached the furthest point
disp$RasterizedSimulationsToFurthest <- lapply(1:nrow(disp), function(x){
  if (length(disp$SimulationsToFurthest[[x]]) > 0){
    rasterized <- suppressMessages(
      rasterizeTerra(vect(disp$SimulationsToFurthest[[x]]), r)
    )
  } else {
    rasterized <- raster(r)
    values(rasterized) <- 0
  }
  return(rasterized)
})

# Visualize simulations to destination
par(mfrow = c(2, 2))
for (n in 1:4){
  plot(sqrt(disp$RasterizedSimulationsToDestin[[n]]), col = viridis(50))
  plot(disp$Tracks[[n]], add = T, col = "white")
}

# Visualize simulations to furthest
par(mfrow = c(2, 2))
for (n in 1:4){
  plot(sqrt(disp$RasterizedSimulationsToFurthest[[n]]), col = viridis(50))
  plot(disp$Tracks[[n]], add = T, col = "white")
}
