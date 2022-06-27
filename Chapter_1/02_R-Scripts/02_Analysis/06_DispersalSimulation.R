################################################################################
#### Simulation of Dispersers Based on the Movement Model
################################################################################
# Description: Here, we will use the movement model that we derived in the
# previous script and simulate dispersers across our study extent

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(davidoff)     # Custom Functions
library(raster)       # For raster manipulation
library(tidyverse)    # For data wrangling
library(glmmTMB)      # For handling logistic model
library(rgdal)        # To handle spatial data
library(velox)        # For quick extraction
library(lubridate)    # To handle timestamps
library(tictoc)       # For simple benchmarking
library(pbmcapply)    # For multicore lapply with progress bar
library(rgeos)        # For manipulation of spatial data

# How many steps per disperser do you want to simulate?
n_steps <- 2000

# How many dispersers per iteration do you want to simulate?
n_points <- 1000

# How many random steps do you want to simulate per realized step?
n_rsteps <- 25

# # How many dispersers do you want to simulate per source point?
# n_disp <- 100

# Do you want to break the simulation of a track if it hits a boundary?
stop <- FALSE

# What is the largest step possible in 4 hours (in meters)?
sl_max <- 35000

################################################################################
#### Load and Prepare the Parametrized Movement Model
################################################################################
# Load the movement model
move.mod <- read_rds("03_Data/03_Results/99_MovementModel.rds")

# Take a look at the models
print(move.mod)

# Select the most parsimonious model
move.mod <- move.mod$Model[[1]]

# Summary
summary(move.mod)

# Visualize it
showCoeffs(getCoeffs(move.mod)[-1, ], xlim = c(-2, 2))

# We will need to simulate random step lengths. For this we will use the gamma
# distribution that we fitted to create random steps in the iSSF model.
sl_dist <- read_rds("03_Data/03_Results/99_GammaDistribution.rds")

################################################################################
#### Load Covariate Layers
################################################################################
# Specify the filenames of the layers we need for prediction. Note that we will
# use averaged layers whenever there is a dynamic counterpart.
layers <- c(
      "01_LandCover_WaterCoverAveraged_MERGED.tif"
    , "01_LandCover_DistanceToWaterAveraged_MERGED.tif"
    , "01_LandCover_TreeCoverAveraged_MODIS.tif"
    , "01_LandCover_NonTreeVegetationAveraged_MODIS.tif"
    , "04_AnthropogenicFeatures_HumanInfluenceBuff_FACEBOOK.grd"
  ) %>%

  # Combine with souce directory
  paste0("03_Data/02_CleanData/", .) %>%

  # Load each file as rasterlayer into a list
  lapply(., stack) %>%

  # Set correct layernames
  set_names(
    c(
        "Water"
      , "SqrtDistanceToWater"
      , "Trees"
      , "Shrubs"
      , "HumansBuff5000"
    )
  )

# Note that we only need to keep one of the human influence layers
layers$HumansBuff5000 <- layers$HumansBuff5000[["Buffer_5000"]]

# We also need to take the sqrt for the "DistanceTo"-layers since we fitted the
# model using the sqrt of these distances
layers[["SqrtDistanceToWater"]] <- sqrt(layers[["SqrtDistanceToWater"]])

# Assign correct layernames
names <- names(layers)
for (i in 1:length(layers)){
  names(layers[[i]]) <- names[i]
}

# Stack the layers
layers <- stack(layers)

# Extent to which we want to expand covariates artifically
ext1 <- extent(layers)
ext2 <- extent(layers) + c(-1, 1, -1, 1) * metersToDegrees(100000)
ext1 <- as(ext1, "SpatialPolygons")
ext2 <- as(ext2, "SpatialPolygons")

# Loop through the layers and extend them
layers <- lapply(1:nlayers(layers), function(x){
  extendRaster(layers[[x]], ext2)
}) %>% stack()

# Visualize the covariates
plot(layers)

# Visualize the original extent
plot(layers[[2]])
plot(ext1, add = T)
plot(ext2, add = T)

################################################################################
#### Define Source Areas
################################################################################
# Load protected areas
# prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
#
# # Visualize them
# u <- unique(prot$Desig)
# m <- match(prot$Desig, u)
# n <- length(unique(prot$Desig))
# pal <- c("#E5F5E0", "#A1D99B", "#31A354")
# plot(ext2)
# plot(ext1, add = T)
# plot(prot, col = pal[m], add = T)
# legend("topleft", legend = u, col = pal, pch = 19)
#
# # Dissolve the borders
# source_areas <- prot %>%
#   aggregate(., dissolve = TRUE) %>%
#   gBuffer(., width = 1 / 110 * 0.001) %>%
#   disaggregate(.) %>%
#   as(., "SpatialPolygonsDataFrame") %>%
#   crop(., ext1)
#
# # Calculate area of each polygon (in km2)
# source_areas$Area <- source_areas %>%
#   spTransform(., CRS("+init=epsg:32734")) %>%
#   gArea(., byid = TRUE) %>%
#   "*" (1e-6)
#
# # Retain only polygons with an area of over 700km2
# plot(subset(source_areas, Area >= 700))
# source_areas <- subset(source_areas, Area >= 700)
#
# # Assign unique ID to each polygon
# source_areas$ID <- 1:nrow(source_areas)
# source_areas$ID <- as.numeric(source_areas$ID)
#
# # Visualize them
# plot(source_areas)
# text(source_areas, source_areas$ID)
#
# # Store them
# writeOGR(source_areas
#   , dsn       = "03_Data/03_Results"
#   , layer     = "99_SourceAreas"
#   , driver    = "ESRI Shapefile"
#   , overwrite = TRUE
# )
#
# # Add a source area around the map border
# buffer <- ext2 - ext1
# crs(buffer) <- crs("+init=epsg:4326")
#
# # Add some info
# buffer$ID <- 1
#
# # Visualize it
# plot(buffer)
# plot(source_areas, add = T)
#
# # Store it
# writeOGR(buffer
#   , dsn       = "03_Data/03_Results"
#   , layer     = "99_BufferArea"
#   , driver    = "ESRI Shapefile"
#   , overwrite = TRUE
# )

# Reload them
buffer <- readOGR("03_Data/03_Results/99_BufferArea.shp")
source_areas <- readOGR("03_Data/03_Results/99_SourceAreas.shp")

# ################################################################################
# #### Sampling Source Points
# ################################################################################
# # We want to sample source points following two different approaches. First, we
# # want to create source points at the same 68 locations that we used to
# # calculate least-cost paths. Second, we want to randomize source points within
# # the catchment area of each of these start points. Let's load the source points
# # and source areas for this first.
# source_areas <- shapefile("/home/david/ownCloud/University/15. PhD/Chapter_0/03_Data/03_Results/99_SourceAreas.shp")
# source_points <- shapefile("/home/david/ownCloud/University/15. PhD/Chapter_0/03_Data/03_Results/99_SourcePoints.shp")
#
# # Rename column names that were abbreviated when stored
# names(source_points)[2] <- "ProtectedArea"
#
# # We now want to cut protected areas into separated polygons such that each
# # source point is allocated to one polygon. We will call this area catchement
# # area. We can create voronoi polygons for this purpose.
# source_areas <- lapply(1:length(unique(source_areas$ID)), function(x){
#
#   # Subset to respective areas and points
#   areas_sub <- subset(source_areas, ID == x)
#   points_sub <- subset(source_points, ProtectedArea == x)
#
#   # Can only tesselate the source area if it contains more than 1 point
#   if (nrow(points_sub) > 1){
#
#     # Create voronoi polygons based on the selected points. We make the extent
#     # larger to make sure the entire source_area is covered by the tesselated
#     # polygon.
#     voris <- voronoi(points_sub, ext = extent(areas_sub) + c(-1, 1, -1, 1))
#
#     # For each voronoi polygon we now create a clipped area
#     areas_sub <- gIntersection(voris, areas_sub, byid = T)
#   }
#
#   # Return the clipped area
#   return(as(areas_sub, "SpatialPolygons"))
#
# }) %>% do.call(rbind, .)
#
# # Assign ID to each area, corresponding to the ID of the point that it contains.
# source_areas$ID <- source_areas %>%
#   over(source_points) %>%
#   dplyr::select(ID) %>%
#   as.matrix() %>%
#   as.vector()
#
# # Add a source area around the map border
# buffer <- ext2 - ext1
# crs(buffer) <- crs("+init=epsg:4326")
#
# # Add relevant information
# buffer$ID <- max(source_areas$ID) + 1
# buffer_point <- spsample(buffer, n = 1, type = "random")
# buffer_point$ID <- max(source_points$ID) + 1
# buffer_point$ProtectedArea <- buffer$ID
# source_areas <- rbind(source_areas, buffer)
# source_points <- rbind(source_points, buffer_point)
#
# # Visualize the results and visually verify that IDs were assigned correctly
# plot(source_areas, border = source_areas$ID)
# plot(source_points, col = source_points$ID, add = T)
# text(x = source_points, labels = source_points$ID, cex = 1)
#
# # Overlay the areas with all the protected areas
# prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
# u <- unique(prot$Desig)
# m <- match(prot$Desig, u)
# n <- length(unique(prot$Desig))
# pal <- c("#E5F5E0", "#A1D99B", "#31A354")
# plot(ext2)
# plot(ext1, add = T)
# plot(prot, col = pal[m], add = T)
# plot(source_points, add = T, col = "red", pch = 20)
# plot(source_areas, add = T, border = "black", lwd = 0.1, lty = 2)
# legend("topleft", legend = u, col = pal, pch = 19)
#
# # Store the source points and source areas
# writeOGR(source_areas
#   , dsn       = "03_Data/03_Results"
#   , layer     = "99_SourceAreas"
#   , driver    = "ESRI Shapefile"
#   , overwrite = TRUE
# )
# writeOGR(source_points
#   , dsn       = "03_Data/03_Results"
#   , layer     = "99_SourcePoints"
#   , driver    = "ESRI Shapefile"
#   , overwrite = TRUE
# )
#
# # Now we can use our custom function to sample random points at those locations
# for (x in c(T, F)){
#   points <- createPoints(
#       points    = source_points
#     , areas     = source_areas
#     , randomize = x
#     , n         = 10
#   )
#   plot(ext2, main = paste0("Randomized = ", x))
#   plot(ext1, add = T)
#   plot(source_areas, add = T, border = "gray")
#   plot(source_points, add = T, col = "red", pch = 20)
#   plot(points, add = T, col = "blue", pch = 20, cex = 0.2)
# }
#
# points <- createPoints(
#     points    = source_points
#   , areas     = source_areas
#   , randomize = T
#   , n         = 100
# )

################################################################################
#### Source Areas
################################################################################
# # Alternative way to sample source areas
# prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
#
# # Subset to national parks only, we'll use those as source areas
# source_areas <- subset(prot, Desig == "National Park")
# source_areas$AreaID <- 1:nrow(source_areas)
#
# # Let's add a source area around the map border
# buffer <- ext2 - ext1
# crs(buffer) <- crs("+init=epsg:4326")
# buffer$Name <- "Buffer"
# buffer$IUCN <- NA
# buffer$Country <- NA
# buffer$Desig <- NA
# buffer$Values <- NA
# buffer$AreaID <- max(source_areas$AreaID + 1)
# source_areas <- rbind(source_areas, buffer)
#
# # Exemplify how we're going to sample source points
# source_points <- lapply(1:nrow(source_areas), function(x){
#   spsample(source_areas[x, ], n = n_disp, type = "random")
# }) %>% do.call(rbind, .)
# source_points$AreaID <- over(source_points, source_areas)$AreaID
# source_points$PointID <- 1:nrow(source_points)
#
# # Visualize the results and visually verify that IDs were assigned correctly
# plot(source_areas, border = source_areas$AreaID)
# plot(source_points, col = source_points$AreaID, add = T)
# text(source_areas, labels = source_areas$AreaID, cex = 0.7)
#
# # Store the source points and areas
# writeOGR(source_areas
#   , dsn       = "03_Data/03_Results"
#   , layer     = "99_SourceAreas"
#   , driver    = "ESRI Shapefile"
#   , overwrite = TRUE
# )

################################################################################
#### Preparation for Simulation
################################################################################
# Prepare raster layers for simulation. This step is for performance reasons
# only. Basically, we calculate some stuff that would otherwise need to be
# calculated in each iteration of the simulation loop. By doing this outside the
# loop we can severely speed up the simulation.
covars <- prepareCovars(layers)

# Prepare movement model for simulation (same reasoning as above)
model <- prepareModel(move.mod)

# Prepare scaling parameters. We will need them to scale covariates extracted
# along steps.
scaling <- read_rds("03_Data/03_Results/99_Scaling.rds")

# Function to simulate dispersal based on a step selection model that was fitted
# in the glmmTMB framework
disperse <- function(
    source              = NULL    # Start Coordinates
  , covars              = NULL    # Spatial Covariates, prepared with our funct.
  , model               = NULL    # iSSF Model, prepared with our funct.
  , sl_dist             = NULL    # Step Length Distribution
  , sl_max              = Inf     # What is the largest possible step?
  , date                = as.POSIXct("2015-06-15 07:00:00", tz = "UTC")
  , n_steps             = 10      # Number of steps simulated
  , n_rsteps            = 25      # Number of random steps proposed
  , scaling             = NULL    # Dataframe to scale extracted covariates
  , stop                = F){     # Should the simulation stop at a boundary?

  # Create a new dataframe indicating the first location. Note that we draw
  # random turning angles to start off
  track <- data.frame(
      x           = coordinates(source)[, 1]
    , y           = coordinates(source)[, 2]
    , absta_      = runif(1, min = 0, max = 2 * pi) # tentative
    , ta_         = NA
    , sl_         = NA
    , Timestamp   = date
    , BoundaryHit = FALSE
    , inactive    = NA
  )

  # Simulate random steps
  for (i in 1:n_steps){

    # Check if the timestamp corresponds to low or high activity
    inactive <- strftime(date, tz = "UTC", format = "%H:%M:%S")
    inactive <- ifelse(inactive %in% c("07:00:00"), 1, 0)

    # Draw random turning angles
    ta_new <- runif(n_rsteps
      , min = - pi
      , max = + pi
    )

    # Draw random step lengths
    sl_new <- rgamma(n_rsteps
      , shape = sl_dist$params$shape
      , scale = sl_dist$params$scale
    )

    # In case the sl_ should be capped, do so
    if (sl_max != Inf){
      sl_new <- pmin(sl_new, sl_max)
    }

    # Identify origin of track
    begincoords <- track[i, c("x", "y")]

    # Calculate new absolute turning angles
    absta_new <- getAbsNewC(
        absta = track$absta_[i]
      , ta    = ta_new
    )

    # Calculate new endpoints
    endpoints_new <- calcEndpointsC(
        xy    = as.matrix(track[i, c("x", "y")])
      , absta = absta_new
      , sl    = sl_new
    )

    # Check which endpoints leave the study extent
    inside <- pointsInside(
        xy     = endpoints_new
      , extent = covars$extent
    )

    # In case some steps are not inside the study area and we want the loop to
    # break
    if (sum(!inside) > 0 & stop){

      # Break the loop
      break

    # In case some steps are not inside the study area and we DONT want the loop
    # to break
    } else if (sum(!inside) > 0 & !stop){

      # Keep only steps inside the study area
      endpoints_new <- endpoints_new[inside, ]
      absta_new     <- absta_new[inside]
      ta_new        <- ta_new[inside]
      sl_new        <- sl_new[inside]

    }

    # Create spatial lines from origin to new coordinates
    l <- vector("list", nrow(endpoints_new))
    for (j in seq_along(l)){
        l[[j]] <- Lines(
          list(
            Line(
              rbind(
                  begincoords[1, ]
                , endpoints_new[j,]
              )
            )
          ), as.character(j)
        )
    }

    # Coerce to spatial lines
    steps <- SpatialLines(l)

    # Extract covariates along each step
    extracted <- extrCov2(covars$covars, steps)

    # Put some nice column names
    names(extracted) <- covars$covar_names

    # Put everything into a dataframe
    rand <- data.frame(
        x           = endpoints_new[, 1]
      , y           = endpoints_new[, 2]
      , absta_      = absta_new
      , ta_         = ta_new
      , sl_         = sl_new
      , BoundaryHit = sum(!inside) > 0
    )

    # Put all covariates into a dataframe. We will use this to calculate
    # selection scores
    covariates <- data.frame(
        extracted
      , cos_ta_  = cos(ta_new)
      , log_sl_  = log(sl_new)
      , sl_      = sl_new
    )

    # Scale covariates
    covariates <- scaleCovars(covariates, scaling)

    # Put the activity phase into the covariate table as well
    covariates$inactive <- inactive

    # Calculate selection scores
    SelectionScore <- as.numeric(predictScore(
        coefficients  = model$coefficients
      , formula       = model$formula
      , data          = covariates
    ))

    # Update date
    date <- date + hours(4)

    # Note that we assume that no fix exists at 11:00. In this case we add
    # another 4 hours
    if(strftime(date, tz = "UTC", format = "%H:%M:%S") == "11:00:00"){
      date <- date + hours(4)
    }

    # Coerce selection scores to probabilities
    Probs <- SelectionScore / sum(SelectionScore)

    # Sample a step according to the above predicted probabilities
    rand <- rand[sample(1:nrow(rand), size = 1, prob = Probs), ]

    # Add updated values to current step
    track$absta_[i]       <- rand$absta_[1]
    track$ta_[i]          <- rand$ta_[1]
    track$sl_[i]          <- rand$sl_[1]
    track$BoundaryHit[i]  <- rand$BoundaryHit[1]
    track$inactive[i]     <- inactive

    # Add new endpoints and new (tentative) absolute turning angle to the next
    # step
    track[i + 1, "x"]         <- rand$x[1]
    track[i + 1, "y"]         <- rand$y[1]
    track[i + 1, "absta_"]    <- rand$absta_[1]
    track[i + 1, "Timestamp"] <- date
  }
  return(track)
}

###############################################################################
### Test of Function
###############################################################################
# # Sample source points
# source_points <- spsample(source_areas, n = 1000, type = "random")
#
# # Visualize
# plot(source_areas)
# plot(source_points, add = T, col = "red", cex = 0.1)
#
# # Try out the function on the first 50 source points
# tic()
# sim <- pbmclapply(1:100
#   , mc.cores           = detectCores() - 1
#   , ignore.interactive = T
#   , FUN                = function(x){
#     si <- disperse(
#         source    = source_points[x, ]
#       , covars    = covars
#       , model     = model
#       , sl_dist   = sl_dist
#       , sl_max    = sl_max
#       , n_rsteps  = n_rsteps
#       , n_steps   = 200
#       , stop      = F
#       , scaling   = scaling
#       , date      = as.POSIXct("2015-06-15 03:00:00", tz = "UTC")
#     )
#     si$TrackID <- x
#     return(si)
# })
# toc()
#
# # Coerce the simulations to tracks
# sims <- lapply(sim, function(x){
#   coordinates(x) <- c("x", "y")
#   l <- spLines(x)
#   return(l)
# }) %>% do.call(rbind, .)
#
# # Bind simulations
# sim <- do.call(rbind, sim)
#
# # Visualize simulations
# plot(crop(layers[[1]], sims))
# plot(sims, add = T)
#
# # Load observed dispersers
# obs <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL(SSF_Extracted).csv")
# obs <- subset(obs, case_)
#
# # Check the distribution of step lengths and turning angles and compare them to
# # the observed data
# cbind(
#     Simulation  = summary(sim$sl_)
#   , Observerd   = summary(obs$sl_)
#   , Tentative   = summary(rgamma(10000
#       , shape = sl_dist$params$shape
#       , scale = sl_dist$params$scale
#     ))
# )
#
# # Visualize
# turns <- data.frame(
#     Group        = c(rep("Simulation", nrow(sim)), rep("Observed", nrow(obs)))
#   , TurningAngle = c(sim$ta_, obs$ta_)
# )
# lengths <- data.frame(
#     Group        = c(rep("Simulation", nrow(sim)), rep("Observed", nrow(obs)))
#   , StepLength   = c(sim$sl_, obs$sl_)
# )
# ggplot(turns, aes(x = cos(TurningAngle), fill = Group)) +
#   geom_density(alpha = 0.2)
# ggplot(lengths, aes(x = StepLength, fill = Group)) +
#   geom_density(alpha = 0.2)
#
# # Compare steps lengths during activity and inactivity
# ggplot(sim, aes(x = factor(inactive), y = sl_)) + geom_boxplot()
#
################################################################################
#### Setting up the Simulation
################################################################################
# Function to easily run and store the simulation
simulateDispersal <- function(filename = NULL){

  # Sample fresh source points
  source_points <- spsample(source_areas, n = n_points, type = "random")

  # Run the simulation for each source point
  tracks <- pbmclapply(
      X                   = 1:length(source_points)
    , mc.cores            = detectCores() - 1
    , ignore.interactive  = T
    , FUN                 = function(x){
      sim <- suppressWarnings(
        disperse(
            source    = source_points[x, ]
          , covars    = covars
          , model     = model
          , sl_dist   = sl_dist
          , sl_max    = sl_max
          , date      = as.POSIXct("2021-01-01 07:00:00", tz = "UTC")
          , n_rsteps  = n_rsteps
          , n_steps   = n_steps
          , scaling   = scaling
          , stop      = stop
        )
      )

      # Assign ID that corresponds to the iteration number
      sim$TrackID <- x

      # Return the simulation
      return(sim)
  })

  # Put tracks together
  tracks <- do.call(rbind, tracks)

  # Convert to tibble (this saves an amazing amount of space)
  tracks <- as_tibble(tracks)

  # Store the simulations
  write_rds(tracks, filename)

}

################################################################################
#### Run Simulation
################################################################################
# Prepare directories into which we can store the results
dir.create("03_Data/03_Results/99_Simulations", showWarnings = F)

# Prepare subfolders
dir.create("03_Data/03_Results/99_Simulations/Main", showWarnings = F)
dir.create("03_Data/03_Results/99_Simulations/Buffer", showWarnings = F)

# Initiate a file to keep track of simulation progress
report <- data.frame(
    iteration    = NA
  , source_areas = NA
  , n_steps      = NA
  , n_points     = NA
  , sl_max       = NA
  , stop         = NA
  , duration     = NA
  , filename     = NA
)
write_csv(report, "03_Data/03_Results/99_Simulations/Report.csv")

Run simulation 10 times
for (i in 1:50){

  # Keep track of duration
  start <- Sys.time()

  # Run simulation
  simulateDispersal(filename  = paste0(
        "03_Data/03_Results/99_Simulations/Main/Iteration_"
      , i
      , ".rds"
    )
  )

  # Update report
  report <- read_csv("03_Data/03_Results/99_Simulations/Report.csv")
  report <- drop_na(rbind(report, data.frame(
      iteration    = i
    , source_areas = "Main"
    , n_steps      = n_steps
    , n_points     = n_points
    , sl_max       = sl_max
    , stop         = stop
    , duration     = difftime(Sys.time(), start, units = "mins")
    , filename = paste0("03_Data/03_Results/99_Simulations/Main/Iteration_", i, ".rds")
  )))
  write_csv(report, "03_Data/03_Results/99_Simulations/Report.csv")

  # Print the duration
  print(Sys.time() - start)

}

# Load Buffer
source_areas <- readOGR("03_Data/03_Results/99_BufferArea.shp")

# Run simulation again
for (i in 1:30){

  # Keep track of duration
  start <- Sys.time()

  # Run simulation
  simulateDispersal(filename  = paste0(
        "03_Data/03_Results/99_Simulations/Buffer/Iteration_"
      , i
      , ".rds"
    )
  )

  # Update report
  report <- read_csv("03_Data/03_Results/99_Simulations/Report.csv")
  report <- drop_na(rbind(report, data.frame(
      iteration    = i
    , source_areas = "Buffer"
    , n_steps      = n_steps
    , n_points     = n_points
    , sl_max       = sl_max
    , stop         = stop
    , duration     = difftime(Sys.time(), start, units = "mins")
    , filename = paste0("03_Data/03_Results/99_Simulations/Main/Iteration_", i, ".rds")
  )))
  write_csv(report, "03_Data/03_Results/99_Simulations/Report.csv")

  # Print the duration
  print(Sys.time() - start)
}
