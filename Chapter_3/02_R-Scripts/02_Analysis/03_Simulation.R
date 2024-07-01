################################################################################
#### Simulation
################################################################################
# Description: Dispersal simulation across static and dynamic landscapes

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_3")

# Load required packages
# library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(Rcpp)           # To load C++ functions
library(pbmcapply)      # For multicore use with progress bar
library(beepr)          # For beep sounds

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Define simulation parameters
n_steps  <- 2000
n_rsteps <- 25
n_indiv  <- 500
sl_max   <- 42000
ext      <- ext(c(21, 27, -21, -17))

################################################################################
#### Preparation
################################################################################
# Load all elements needed for the simulation
gamma  <- read_rds("03_Data/03_Results/StepLengthDistribution.rds")
forms  <- read_rds("03_Data/03_Results/Formula.rds")
scal   <- read_rds("03_Data/03_Results/Scaling.rds")
model  <- read_rds("03_Data/03_Results/MovementModels.rds") %>% subset(NumberRandomSteps == 100)
light  <- read_rds("03_Data/02_CleanData/Moonlight.rds") %>% dplyr::select(c(Timestamp, LightType)) %>% mutate(LightType = as.factor(LightType))
source <- vect("03_Data/02_CleanData/Sources.gpkg")

# Load covariates and associated lookup table
covs <- read_rds("03_Data/02_CleanData/Covariates.rds") %>% dplyr::select(-Dates)
look <- read_rds("03_Data/02_CleanData/LookupTable.rds")

# We'll use the dynamic_aggregated data this time, not the dynamic one
covs <- subset(covs, Type != "Dynamic")
look <- subset(look, Type != "Dynamic")

# Dates for which we'll simulate diseprsal
simdates <- look %>%
  subset(year(Timestamp) == 2000) %>%
  pull(Timestamp)

# Cretae a covariates object from the covariates and lookup table
covs <- covariates(covs, look)
show(covs)
rm(look)

# Drop covariates that we don't need
covs <- subset(covs, !(Covariate %in% c("NDVI", "DistanceToPans", "Forest")))

# There are a couple of useful functions that we can apply to this object. Let's
# try them.
test <- subset(covs, Type == "Static" & Covariate %in% c("Humans", "Trees"))
test <- loadCovariates(test, readall = F)
extractCovariates(test, t(c(24, -19)))
getCovariates(covs, type = "DynamicAggregated", covariate = "Water", timestamps = simdates[1])

# Sample points within each source area. Start dates will be distributed equally
# across the year 2000
set.seed(12345)
start <- source %>% spatSample(size = n_indiv, strata = "Name")
start <- tibble(Source = start$Name, as.data.frame(crds(start))) %>%
  group_by(Source) %>%
  mutate(Replicate = 1:n()) %>%
  mutate(Timestamp = simdates[seq(1, length(simdates), length.out = length(Source))]) %>%
  ungroup()

# Design through which we'll loop
design <- tibble(
      FittingCovariates    = c("Static", "Dynamic")
    , ModelSeasons         = c("Single", "Multi")
    , PredictionCovariates = c("Static", "DynamicAggregated")
    , ModelCode            = c("SSS", "DMD")
  ) %>% expand_grid(
      Source    = unique(source$Name)
    , Formula   = c("Simple", "Full")
    , Replicate = 1:n_indiv
    , .
  ) %>%
  mutate(ModelCode = paste0(ModelCode, "_", substr(Formula, start = 1, stop = 1))) %>%
  left_join(., start, by = c("Source", "Replicate")) %>%
  mutate(Filename = paste0(
      "03_Data/03_Results/Simulations/"
    , ModelCode, "_"
    , Source, "_"
    , sprintf("%04d", Replicate)
    , ".rds")
  ) %>%
  mutate(Done = file.exists(Filename))

# Write the design to file
write_rds(design, "03_Data/03_Results/Simulations.rds")

# Create folder
if (!dir.exists("03_Data/03_Results/Simulations/")) {
  dir.create("03_Data/03_Results/Simulations/")
}

################################################################################
#### Simulation Function
################################################################################
# Function to simulate dispersal from a fitted iSSF model
simulateDispersal <- function(
      x                   # Start longitude
    , y                   # Start latitude
    , timestamps          # Timestamps for which to simulate
    , covars              # Spatial covariates (a covariates object)
    , light               # Light covariates
    , formula             # Prediction formula
    , betas               # Selection coefficients (for both seasons)
    , shape               # Shape of step length distribution
    , scale               # Scale of step length distribution
    , scaling             # Dataframe to scale extracted covariates
    , n_rsteps = 25       # Number of random steps proposed
    , n_steps  = 100      # Number of steps simulated
    , sl_max   = Inf      # What is the largest possible step?
    , stop     = F        # Should the simulation stop when a boundary is reached?
    , type     = "Static" # If static or dynamic covariates should be used
    , seasons  = "Single" # If single-season or multi-season models should be used
    , panex    = F        # If extent should be taken from the pan maps
    , progress = F        # If progress bar should be printed
    , debug    = F        # Debug mode
    , exte     = NULL     # If you want to provide a simulation extent
  ) {

  # For testing
  # k          <- 2
  # x          <- (xmin(ext) + xmax(ext)) / 2
  # y          <- ymin(ext)
  # timestamps <- dispdates
  # covars     <- covs
  # light      <- light
  # formula    <- model_sub$formula
  # betas      <- model_sub$betas
  # shape      <- gamma$shape
  # scale      <- gamma$scale
  # scaling    <- scal_sub
  # n_rsteps   <- n_rsteps
  # n_steps    <- n_steps
  # sl_max     <- sl_max
  # stop       <- F
  # type       <- design$PredictionCovariates[k]
  # seasons    <- design$ModelSeasons[k]
  # panex      <- F
  # progress   <- T
  # debug      <- T
  # exte       <- ext

  # Prepare covariates for prediction
  covars_sub <- prepareCovariates(
      covars     = covars
    , timestamps = timestamps
    , type       = type
    , readall    = F
  )

  # If no extent is provided, obtain it from the covariate layers
  if (is.null(exte)) {
    exte <- covariatesExtent(covars_sub)
  }
  exte <- list(
      xmin = xmin(exte)
    , xmax = xmax(exte)
    , ymin = ymin(exte)
    , ymax = ymax(exte)
  )

  # Initiate a track of the desired length. This will preallocate the necessary
  # memory (at least once we compute step metrics)
  n_steps <- length(timestamps) + 1
  track <- data.frame(
      Timestamp = c(ymd_hms(NA), timestamps)
    , x         = as.numeric(rep(NA, n_steps))
    , y         = as.numeric(rep(NA, n_steps))
  )

  # Prepare a track with a random origin (which is needed so we have a heading)
  # We'll also project the coordinates for now.
  track$x[2] <- x
  track$y[2] <- y
  track$x[1] <- x + runif(1, -1, 1)
  track$y[1] <- y + runif(1, -1, 1)
  track      <- projectCoords(track)

  # Compute step metrics and add them to the steps
  track <- cbind(track, stepMetrics(track))

  # Add light statistics
  track <- left_join(track, light, by = "Timestamp")

  # Determine season
  if (seasons == "Single") {
      track$Season <- "All"
    } else {
      track$Season <- getRainySeason(track$Timestamp)
  }

  # Give simulation summary and print progress bar
  if (progress) {
    cat("\n")
    cat("-------------------------------------------\n")
    cat("Running simulation with following settings:\n")
    cat("-------------------------------------------\n")
    cat("- Prediction covariates   :     ", type, "\n")
    cat("- Model seasons           :     ", seasons, "\n")
    pb <- txtProgressBar(0, n_steps, style = 3)
  }

  # Loop through the desired number of steps and simulate movement
  for (i in 2:n_steps) {

    # Testing
    # i <- 2

    # Draw random turning angles
    ta_new <- runif(n_rsteps
      , min = - pi
      , max = + pi
    )

    # Draw random step lengths
    sl_new <- rgamma(n_rsteps
      , shape = gamma$shape
      , scale = gamma$scale
    )

    # In case the sl_ should be capped, do so
    if (sl_max != Inf) {
      sl_new <- pmin(sl_new, sl_max)
    }

    # Put the step lengths and turning angles into a new dataframe. These are
    # our proposed random steps.
    rand <- data.frame(
        sl     = sl_new
      , absta  = newabsTA(track$absta[i - 1], ta_new)
      , ta     = ta_new
    )

    # Calculate endpoints of each random step
    rand$x_to <- track$x[i] + sin(rand$absta) * rand$sl
    rand$y_to <- track$y[i] + cos(rand$absta) * rand$sl
    # rand      <- unprojectCoords(rand)

    # Subset to the appropriate covariate layers
    covars_sub_sub <- subset(covars_sub, Timestamp == track$Timestamp[i])

    # Check if any endpoints leave the study extent
    if (panex) {

        # Get date of the pan-map used
        pan_date <- covars_sub_sub@lookup %>%
          subset(Covariate == "DistanceToPans") %>%
          pull(Layerdate)

        # Get associated moving window
        wind <- windows %>%
          subset(Year == year(pan_date) & Month == month(pan_date)) %>%
          pull(Window) %>%
          do.call(rbind, .) %>%
          vect()

        # Check which points lie within the moving window
        inside <- rand[, c("x_to", "y_to")] %>%
          unprojectCoords() %>%
          as.matrix() %>%
          vect() %>%
          relate(wind, ., "contains") %>%
          as.vector()

      } else {
        inside <- pointsInside(
            xy     = unprojectCoords(rand[, c("x_to", "y_to")])
          , extent = exte
        )
    }

    # In case some steps are not inside the study area and we want the loop to
    # break
    if (sum(!inside) > 0 & stop) {
        break
      } else if (sum(inside) == 0) {
        if (debug) {
          cat("No random steps inside the extent left. Breaking the simulation.\n")
        }
        break
      } else {
        rand <- rand[inside, ]
    }
    if (debug) {
      cat(nrow(rand), "random steps left\n")
    }

    # Generate interpolated points between start and end of each step
    interpolated <- lapply(seq_len(nrow(rand)), function(j) {
      p <- interpolatePoints(
          x1 = track$x[i]
        , x2 = rand$x_to[j]
        , y1 = track$y[i]
        , y2 = rand$y_to[j]
        , by = 250
      )
      p <- as.data.frame(p)
      p$step <- j
      return(p)
    }) %>% do.call(rbind, .) %>% unprojectCoords()

    # Extract covariates and average along each step
    extr <- extractCovariates(
          covariates = covars_sub_sub
        , xy         = as.matrix(interpolated[, c("x", "y")])
      )
    if (debug & (nrow(extr) == 0 | !is.data.frame(extr))) {
      browser()
    }
    extr <- cbind(step = interpolated$step, extr) %>%
      group_by(step) %>%
      summarise_all(mean, na.rm = T)

    # Compute sqrt distances
    extr$SqrtDistanceToWater <- sqrt(extr$DistanceToWater)
    # extr$SqrtDistanceToPans  <- sqrt(extr$DistanceToPans)

    # Find entries that contain NAs and drop them
    indices <- which(is.na(extr), arr.ind = T)[, 1]
    if (length(indices) > 0) {
      extr <- extr[-indices, ]
      rand <- rand[-indices, ]
    }

    # Put all covariates into a dataframe. We will use this to calculate
    # selection scores
    covariates <- data.frame(
        extr[, -1]
      , cos_ta    = cos(rand$ta)
      , log_sl    = log(rand$sl)
      , sl        = rand$sl
    )

    # Scale covariates
    if (debug & (nrow(covariates) == 0 | !is.data.frame(covariates))) {
      beep()
      Sys.sleep(1)
      beep()
      Sys.sleep(1)
      beep()
      browser()
    }
    covariates <- scaleCovars(covariates, scaling)

    # Add the lighttype
    covariates$LightType <- track$LightType[i]

    # Prepare model matrix (and remove intercept)
    mat <- model.matrix(formula, covariates)
    mat <- mat[ , -1]

    # Calculate selection scores
    score <- exp(mat %*% betas[[track$Season[i]]])

    # Coerce selection scores to probabilities
    probs <- score / sum(score)

    # Sample a step according to the above predicted probabilities
    rand <- rand[sample(seq_len(nrow(rand)), size = 1, prob = probs), ]

    # Add updated values to current step
    track$x_to[i]  <- rand$x_to[1]
    track$y_to[i]  <- rand$y_to[1]
    track$sl[i]    <- rand$sl[1]
    track$absta[i] <- rand$absta[1]
    track$relta[i] <- rand$ta[1]

    # Add new endpoints and new (tentative) absolute turning angle to the next
    # step
    if (i == n_steps) break
    track$x[i + 1] <- track$x_to[i]
    track$y[i + 1] <- track$y_to[i]

    # Update progress
    if (progress) {
      setTxtProgressBar(pb, i)
    }

  }

  # Drop the initial pseudofix and unproject data
  track           <- track[-1, ]
  track           <- unprojectCoords(track)
  track$LightType <- NULL

  # Return the final track
  return(track)
}

################################################################################
#### Simulation
################################################################################
# Subset to combinations that haven't been run yet and randomize the design
design <- subset(design, !Done)
design <- slice_sample(design, n = nrow(design), replace = F)
cat(nrow(design), "Simulations to go through\n")

# Loop through the design and simulate dispersers for the remaining combinations
pbmclapply(
    X                  = seq_len(nrow(design))
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , mc.set.seed        = T
  , FUN                = function(k) {
    # for (k in seq_len(nrow(design))) {

  # TESTING
  # k <- 2

  # Prepare vector of dispersal dates for which we'll simulate movement
  dispdates <- tibble(Timestamp = seq(
        from       = design$Timestamp[k]
      , by         = "4 hours"
      , length.out = n_steps * 1.5
    )) %>%
    subset(hour(Timestamp) %in% c(3, 7, 15, 19, 23)) %>%
    slice(1:n_steps) %>%
    pull(Timestamp) %>%
    update(year = 2000)

  # Select correct formula
  form <- forms[[design$Formula[k]]]

  # Subset to the appropriate model and prepare it for prediction
  model_sub <- subset(model
    , FittingCovariates == design$FittingCovariates[k]
    & ModelSeasons      == design$ModelSeasons[k]
    & Formula           == design$Formula[k]
  ) %>% prepareModels(form, .)

  # Subset the scaling table to the relevant entries
  if (design$PredictionCovariates[k] == "Static") {
      keep <- !grepl(scal$Covariate, pattern = "Dynamic")
      scal_sub <- scal[keep, ]
    } else {
      keep <- !grepl(scal$Covariate, pattern = "Static")
      scal_sub <- scal[keep, ]
  }
  scal_sub$Covariate <- gsub(scal_sub$Covariate, pattern = "Dynamic|Static", replacement = "")

  # Simulate individual
  # i <- 1
  # while (i < 100) {
    # tic <- Sys.time()
    sim <- simulateDispersal(
        x          = design$x[k]
      , y          = design$y[k]
      #   x          = (xmin(ext) + xmax(ext)) / 2
      # , y          = ymin(ext)
      , timestamps = dispdates
      , covars     = covs
      , formula    = model_sub$formula
      , betas      = model_sub$betas
      , shape      = gamma$shape
      , scale      = gamma$scale
      , scaling    = scal_sub
      , n_rsteps   = n_rsteps
      , sl_max     = sl_max
      , stop       = F
      , type       = design$PredictionCovariates[k]
      , light      = light
      , seasons    = design$ModelSeasons[k]
      , progress   = F
      , debug      = F
      , exte       = ext
    )
    # toc <- Sys.time()
  #   cat("\n Processing time:", round(toc - tic, 3), "\n")
  #   i <- i + 1
  # }

  # Store it to file
  write_rds(sim, design$Filename[k])
  return(design$Filename[k])
})

# # Extend the covariates artificially (this will load them into memory)
# ext <- ext(water) * 1.2
# water      <- extendRaster(water, ext)
# trees      <- extendRaster(trees, ext)
# shrub      <- extendRaster(shrub, ext)
# human      <- extendRaster(human, ext)
# dista_sqrt <- extendRaster(dista_sqrt, ext)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/02_Analysis/03_Simulation.rds")
cat("Done :)\n")
