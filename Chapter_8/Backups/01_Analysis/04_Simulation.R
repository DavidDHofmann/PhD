################################################################################
#### Simulation
################################################################################
# Description: Simulation of dispersal through seasonal landscapes

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(velox)          # For quick extraction
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(glmmTMB)        # To handle the movement model
library(Rcpp)           # To load C++ functions
library(pbmcapply)      # For multicore use with progress bar

# Load custom functions
source("02_R-Scripts/00_Functions.R")
sourceCpp("02_R-Scripts/00_Functions.cpp")

################################################################################
#### Load Required Data
################################################################################
# Load the calibrated movement model
mod <- read_rds("03_Data/02_CleanData/MovementModel.rds")

# Load the scaling parameters
scaling <- read_rds("03_Data/02_CleanData/Scaling.rds")

# Load the fitted gamma distribution (for sampling step lengths)
sl_dist <- read_rds("03_Data/02_CleanData/GammaDistribution.rds")

# Load the habitat layers
cat("Preparing covariate layers...\n")
water <- rast("03_Data/02_CleanData/WaterCover.tif")
dista <- rast("03_Data/02_CleanData/DistanceToWater.tif")
trees <- rast("03_Data/02_CleanData/TreeCover.tif")
shrub <- rast("03_Data/02_CleanData/ShrubCover.tif")
human <- rast("03_Data/02_CleanData/HumanInfluence.tif")
prote <- vect("03_Data/02_CleanData/Protected.shp")
areas <- vect("03_Data/02_CleanData/SourceAreas.shp")

# We are interested in the square rooted distance to water
dista_sqrt <- sqrt(dista)

# Extend the covariates artificially
ext <- ext(water) * 1.2
water <- extendRaster(water, ext)
trees <- extendRaster(trees, ext)
shrub <- extendRaster(shrub, ext)
human <- extendRaster(human, ext)
dista_sqrt <- extendRaster(dista_sqrt, ext)

# Combine the rasters of the same flood-level together
cov <- list(
    Min  = c(water[["min"]], dista_sqrt[["min"]], trees[["min"]], shrub[["min"]], human)
  , Mean = c(water[["mean"]], dista_sqrt[["mean"]], trees[["mean"]], shrub[["mean"]], human)
  , Max  = c(water[["max"]], dista_sqrt[["max"]], trees[["max"]], shrub[["max"]], human)
)

# Make nicer layernames
names(cov[[1]]) <- names(cov[[2]]) <- names(cov[[3]]) <-
  c("Water", "SqrtDistanceToWater", "Trees", "Shrubs", "HumansBuff5000")

# Apply our custom functions to make everything ready for the simulation
mod <- prepareModel(mod)
cov[["Min"]]  <- prepareCovars(cov[["Min"]])
cov[["Mean"]] <- prepareCovars(cov[["Mean"]])
cov[["Max"]]  <- prepareCovars(cov[["Max"]])

################################################################################
#### Simulation Function
################################################################################
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
  , stop                = F) {    # Should the simulation stop at a boundary?

  # # For testing
  # source <- cbind(x = c(23), y = c(-19))
  # covars <- cov$Min
  # model <- mod
  # date <- as.POSIXct("2015-06-15 07:00:00", tz = "UTC")
  # sl_max <- 35000
  # n_rsteps <- 25
  # stop <- F

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
  for (i in 1:n_steps) {

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
    if (sl_max != Inf) {
      sl_new <- pmin(sl_new, sl_max)
    }

    # Identify origin of track
    begincoords <- track[i, c("x", "y")]

    # Calculate new absolute turning angles
    absta_new <- getAbsNew(
        absta = track$absta_[i]
      , ta    = ta_new
    )

    # Calculate new endpoints
    endpoints_new <- calcEndpoints(
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
    if (sum(!inside) > 0 & stop) {

        # Break the loop
        break

      # In case some steps are not inside the study area and we DONT want the
      # loop to break
      } else if (sum(!inside) > 0 & !stop) {

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
    extracted <- extrCov(covars$covars, steps)

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

    # Prepare model matrix (and remove intercept)
    mat <- model.matrix(model$formula, covariates)
    mat <- mat[ , -1]

    # Calculate selection scores
    score <- exp(mat %*% model$coefficients)

    # Update date
    date <- date + hours(4)

    # Note that we assume that no fix exists at 11:00. In this case we add
    # another 4 hours
    if(strftime(date, tz = "UTC", format = "%H:%M:%S") == "11:00:00") {
      date <- date + hours(4)
    }

    # Coerce selection scores to probabilities
    Probs <- score / sum(score)

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

################################################################################
#### Simulation Setup
################################################################################
# Prepare directories into which we can store the results
dir.create("03_Data/03_Results/99_Simulations", showWarnings = F)

# Number of simulated steps
n_steps <- 2000

# How many dispersers do you want to simulate per iteration?
n_points <- 100

# How many iterations (per flood level and source area) do you want to run?
iterations <- 10

# How many random steps do you want to simulate per realized step?
n_rsteps <- 25

# Do you want to break the simulation of a track if it hits a boundary?
stop <- F

# What is the largest step possible in 4 hours (in meters)?
sl_max <- 35000

# Prepare design through which we can loop
design <- expand_grid(
    FloodLevel = c("Min", "Mean", "Max")
  , Iteration  = 1:iterations
  , SourceArea = unique(areas$ID)
)

# Also prepare filenames for the stored simulations
design$Filename <- with(design, paste0(
    "03_Data/03_Results/99_Simulations/Flood_"
  , FloodLevel
  , "_Source_"
  , sprintf("%02d", SourceArea)
  , "_Rep_"
  , sprintf("%02d", Iteration)
  , ".rds"
))
design$Done <- file.exists(design$Filename)

# Initiate a file to keep track of simulation progress
if (!file.exists("03_Data/03_Results/99_Simulations/Report.csv")) {
  report <- data.frame(
      FloodLevel   = NA
    , SourceArea   = NA
    , Iteration    = NA
    , n_steps      = NA
    , n_points     = NA
    , sl_max       = NA
    , stop         = NA
    , duration_min = NA
    , filename     = NA
  )
  write_csv(report, "03_Data/03_Results/99_Simulations/Report.csv")
}

################################################################################
#### Run Simulation
################################################################################
# Go through the design and run the simulation
cat("Simulating dispersal...\n")
for (i in 1:nrow(design)) {

  # Run through the design
  if (!design$Done[i]) {

    # Keep track of duration
    start <- Sys.time()

    # Generate source points
    source_points <- spatSample(areas[areas$ID == design$SourceArea[i]]
      , size   = n_points
      , method = "random"
    )
    source_points <- as(source_points, "Spatial")

    # Run the simulation for each source point
    tracks <- pbmclapply(
        X                   = 1:length(source_points)
      , mc.cores            = detectCores() - 1
      , ignore.interactive  = T
      , FUN                 = function(x) {
        sim <- suppressWarnings(
          disperse(
              source    = source_points[x, ]
            , covars    = cov[[design$FloodLevel[i]]]
            , model     = mod
            , sl_dist   = sl_dist
            , sl_max    = sl_max
            , date      = as.POSIXct("2021-01-01 07:00:00", tz = "UTC")
            , n_rsteps  = n_rsteps
            , n_steps   = n_steps
            , scaling   = scaling
            , stop      = stop
          )
        )

        # Assign some more information
        sim$TrackID <- x

        # Return the simulation
        return(sim)
    })

    # Put tracks together
    tracks <- do.call(rbind, tracks)

    # Convert to tibble (this saves an amazing amount of space)
    tracks <- as_tibble(tracks)

    # Store the simulations
    write_rds(tracks, design$Filename[i])

    # Update the report
    report <- suppressMessages(
      read_csv("03_Data/03_Results/99_Simulations/Report.csv")
    )
    report <- drop_na(rbind(report, data.frame(
        FloodLevel   = design$FloodLevel[i]
      , SourceArea   = design$SourceArea[i]
      , Iteration    = design$Iteration[i]
      , n_steps      = n_steps
      , n_points     = n_points
      , sl_max       = sl_max
      , stop         = stop
      , duration_min = difftime(Sys.time(), start, units = "mins")
      , filename     = design$Filename[i]
    )))
    write_csv(report, "03_Data/03_Results/99_Simulations/Report.csv")
  }

  # Print progress
  cat("Run", i, "out of", nrow(design), "done...\n")

}

################################################################################
#### Combine Simulations
################################################################################
# Once the simulations are done, let's combine them into a single big dataframe
sims <- lapply(1:nrow(design), function(x) {
  dat <- read_rds(design$Filename[x])
  dat$SimID      <- x
  dat$FloodLevel <- design$FloodLevel[x]
  dat$SourceArea <- design$SourceArea[x]
  return(dat)
}) %>% do.call(rbind, .)

# Because in each simulation we start off with new IDs for trajectories, they
# are not unique across simulations. We thus combine ID and SimID to create an
# ID that is unqiue to each simulated path, across all simulations
sims <- sims %>%
  group_by(SimID, TrackID) %>%
  mutate(TrackID = cur_group_id())

# Make sure it worked
table(table(sims$TrackID))

# Collect garbage
gc()

# Let's also create a step counter, indicating the number of the step in its
# respective trajectory
sims <- sims %>%
  group_by(TrackID) %>%
  mutate(StepNumber = row_number())

# Check object size
format(object.size(sims), units = "Gb")

# Ungroup
sims <- ungroup(sims)

# Write to an rds
write_rds(sims, "03_Data/03_Results/DispersalSimulation.rds")

# Remove separate files to free some space
# file.remove(design$Filename)
