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

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Define simulation parameters
n_steps  <- 10
n_rsteps <- 100
n_indiv  <- 10
sl_max   <- 42000

################################################################################
#### Preparation
################################################################################
# Load all elements needed for the simulation
disps <- read_csv("03_Data/02_CleanData/SSF.csv")
valid <- read_rds("03_Data/02_CleanData/ValidationDispersers.rds")
gamma <- read_rds("03_Data/03_Results/StepLengthDistribution.rds")
form  <- read_rds("03_Data/03_Results/ModelFormula.rds")
scal  <- read_rds("03_Data/03_Results/Scaling.rds")
model <- read_rds("03_Data/03_Results/MovementModels.rds")
light <- read_rds("03_Data/02_CleanData/Moonlight.rds") %>% dplyr::select(c(Timestamp, LightType)) %>% mutate(LightType = as.factor(LightType))

# Use the covariates and lookup table to create a "covariates" object, which is
# a custom class that I created to efficienctly work with such raster data.
covs <- read_rds("03_Data/02_CleanData/Covariates.rds") %>% dplyr::select(-Dates)
look <- read_rds("03_Data/02_CleanData/LookupTable.rds")

# Cretae a covariates object
covs <- covariates(covs, look)
show(covs)
rm(look)

# There are a couple of useful functions that we can apply to this object. Let's
# try them.
test <- subset(covs, Type == "Static" & Covariate %in% c("Humans", "Forest"))
test <- loadCovariates(test)
extractCovariates(test, t(c(24, -19)))

# We'll use only the first location of each validation disperser as source point
# Therefore, let's extract them
sourc <- disps %>%
  subset(ID %in% valid) %>%
  subset(case == 1) %>%
  group_by(ID) %>%
  slice(1) %>%
  dplyr::select(ID, x, y, TimestampRounded)

# Specify the design through which we will loop
if (!dir.exists("03_Data/03_Results/Simulations")) {
  dir.create("03_Data/03_Results/Simulations")
}
outfile <- "03_Data/03_Results/Simulations.rds"
if (!file.exists(outfile)) {
    design <- expand_grid(
        FittingCovariates    = c("Static", "Dynamic")
      , ModelSeasons         = c("Single", "Multi")
      , PredictionCovariates = c("Static", "Dynamic")
      , ID                   = valid
      ) %>%
      subset(!(FittingCovariates == "Static" & PredictionCovariates == "Dynamic")) %>%
      left_join(., sourc, by = "ID") %>%
      mutate(ModelCode = paste0(
          substr(FittingCovariates, 1, 1)
        , substr(ModelSeasons, 1, 1)
        , substr(PredictionCovariates, 1, 1)
        , "_"
        , ID
      )) %>%
      mutate(Filename = paste0("03_Data/03_Results/Simulations/", ModelCode, ".rds")) %>%
      mutate(Done = file.exists(Filename))
    write_rds(design, outfile)
  } else {
    design <- read_rds("03_Data/03_Results/Simulations.rds")
    design$Done <- file.exists(design$Filename)
}

################################################################################
#### Simulation Function
################################################################################
# Function to simulate dispersal from a fitted iSSF model
simulateDispersal <- function(
      x                   # Start longitude
    , y                   # Start latitude
    , timestamps          # Timestamps for which to simulate (gives the number of steps simulated)
    , covars              # Spatial covariates
    , light               # Light covariates
    , formula             # Prediction formula
    , betas               # Selection coefficients (for both seasons)
    , shape               # Shape of step length distribution
    , scale               # Scale of step length distribution
    , scaling             # Dataframe to scale extracted covariates
    , n_rsteps = 25       # Number of random steps proposed
    , sl_max   = Inf      # What is the largest possible step?
    , stop     = F        # Should the simulation stop when a boundary is reached?
    , type     = "Static" # If static or dynamic covariates should be used
    , seasons  = "Single" # If single-season or multi-season models should be used
  ) {

  # # For testing
  # k <- 1
  # x          <- design$x[k]
  # y          <- design$y[k]
  # timestamps <- dispdates
  # covars     <- covs
  # formula    <- model_sub$formula
  # betas      <- model_sub$betas
  # shape      <- gamma$shape
  # scale      <- gamma$scale
  # scaling    <- scal_sub
  # n_rsteps   <- 25
  # sl_max     <- 40000
  # stop       <- F
  # type       <- design$PredictionCovariates[k]
  # light      <- light
  # seasons    <- design$ModelSeasons[k]

  # Prepare covariates for prediction
  covars_sub <- prepareCovariates(
      covars     = covars
    , timestamps = timestamps
    , type       = type
  )

  # Obtain an extent from the covariates
  exte <- ext(covariatesExtent(covars_sub))
  exte <- project(exte, "epsg:4326", "epsg:32734")
  exte <- list(
      xmin = xmin(exte)
    , xmax = xmax(exte)
    , ymin = ymin(exte)
    , ymax = ymax(exte)
  )

  # Initiate a track of the desired length. This will preallocate the necessary
  # memory
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
  track[1:2, c("x", "y")] <- projectCoords(track[1:2, c("x", "y")])

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

  # Loop through the desired number of steps and simulate movement
  for (i in 2:n_steps) {

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

    # Check if any endpoints leave the study extent
    inside <- pointsInside(
        xy     = rand[, c("x_to", "y_to")]
      , extent = exte
    )

    # In case some steps are not inside the study area and we want the loop to
    # break
    if (sum(!inside) > 0 & stop) {
        break
      } else if (sum(!inside) > 0 & !stop) {
        rand <- rand[inside, ]
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
          covariates = subset(covars_sub, Timestamp == track$Timestamp[i])
        , xy         = as.matrix(interpolated[, c("x", "y")])
      ) %>%
      cbind(step = interpolated$step, .) %>%
      group_by(step) %>%
      summarise_all(mean, na.rm = T)
    #
    # # Compute sqrt distances
    # extr$SqrtDistanceToWater <- sqrt(extr$DistanceToWater)
    # extr$SqrtDistanceToPans  <- sqrt(extr$DistanceToPans)
    #
    # # Find entries that contains NAs and drop them
    # indices <- which(is.na(extr), arr.ind = T)[, 1]
    # if (length(indices) > 0) {
    #   extr <- extr[-indices, ]
    #   rand <- rand[-indices, ]
    # }
    #
    # # Put all covariates into a dataframe. We will use this to calculate
    # # selection scores
    # covariates <- data.frame(
    #     extr[, -1]
    #   , cos_ta    = cos(rand$ta)
    #   , log_sl    = log(rand$sl)
    #   , sl        = rand$sl
    # )
    #
    # # Scale covariates
    # covariates <- scaleCovars(covariates, scaling)
    #
    # # Add the lighttype
    # covariates$LightType <- track$LightType[i]
    #
    # # Prepare model matrix (and remove intercept)
    # mat <- model.matrix(formula, covariates)
    # mat <- mat[ , -1]
    #
    # # Calculate selection scores
    # score <- exp(mat %*% betas[[track$Season[i]]])
    #
    # # Coerce selection scores to probabilities
    # probs <- score / sum(score)
    #
    # # Sample a step according to the above predicted probabilities
    # rand <- rand[sample(seq_len(nrow(rand)), size = 1, prob = probs), ]
    #
    # ################################################################################
    # #### REMOVE THIS!!!
    # ################################################################################
    rand <- rand[sample(seq_len(nrow(rand)), size = 1), ]

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
  }

  # Drop the initial pseudofix
  track <- track[-1, ]

  # Return the final track
  return(track)
}

################################################################################
#### Simulation
################################################################################
# Subset to combinations that haven't been run yet
# design <- subset(design, !Done)
design <- subset(design, PredictionCovariates == "Dynamic" & ModelSeasons == "Multi" & ID == "Liuwa")

# Loop through the design and simulate dispersers for the remaining combinations
for (k in seq_len(nrow(design))) {
  k <- 1

  # Print update
  cat("Running simulations for", design$ModelCode[k], "\n")

  # Prepare vector of dispersal dates for which we'll simulate movement
  dispdates <- tibble(Timestamp = seq(
        from       = design$TimestampRounded[k]
      , by         = "4 hours"
      , length.out = n_steps * 1.5
    )) %>%
    subset(hour(Timestamp) %in% c(3, 7, 15, 19, 23)) %>%
    slice(1:n_steps) %>%
    pull(Timestamp)

  # Subset to the appropriate model and prepare it for prediction
  model_sub <- subset(model
    , FittingCovariates == design$FittingCovariates[k]
    & ModelSeasons      == design$ModelSeasons[k]
  ) %>% prepareModel(form, .)

  # Subset the scaling table to the relevant entries
  scal_sub <- lapply(scal, function(y) {
    nams <- names(y)
    if (design$PredictionCovariates[k] == "Static") {
        keep <- !sapply(nams, function(y) {grepl(y, pattern = "Dynamic")})
        subi <- y[keep]
      } else {
        keep <- !sapply(nams, function(y) {grepl(y, pattern = "Static")})
        subi <- y[keep]
    }
    names(subi) <- gsub(nams[keep], pattern = "Dynamic|Static", replacement = "")
    return(subi)
  })

  # Simulate movement on the selected covariates
  # sapply(sims, is.data.frame)
  sims <- tryCatch({
    pbmclapply(
        X                  = seq_len(n_indiv)
      , ignore.interactive = T
      , mc.cores           = detectCores() - 1
      , FUN                = function(l) {
        sim <- simulateDispersal(
            x          = design$x[k]
          , y          = design$y[k]
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
        )
        sim$ID <- l
        return(sim)
    })
  }, error = function(e) {return(e)})

  # # sims <- pbmclapply(
  # sims <- lapply(
  #     X                  = seq_len(n_indiv)
  #   # , ignore.interactive = T
  #   # , mc.cores           = detectCores() - 1
  #   , FUN                = function(l) {
  #     sim <- simulateDispersal(
  #         x          = design$x[k]
  #       , y          = design$y[k]
  #       , timestamps = dispdates
  #       , covars     = covs
  #       , formula    = model_sub$formula
  #       , betas      = model_sub$betas
  #       , shape      = gamma$shape
  #       , scale      = gamma$scale
  #       , scaling    = scal_sub
  #       , n_rsteps   = n_rsteps
  #       , sl_max     = sl_max
  #       , stop       = F
  #       , type       = design$PredictionCovariates[k]
  #       , light      = light
  #       , seasons    = design$ModelSeasons[k]
  #     )
  #     sim$ID <- l
  #     return(sim)
  # }) %>% do.call(rbind, .)
  #
  # # Write the simulations to file
  # write_rds(sims, design$Filename[k])

}

# # Extend the covariates artificially (this will load them into memory)
# ext <- ext(water) * 1.2
# water      <- extendRaster(water, ext)
# trees      <- extendRaster(trees, ext)
# shrub      <- extendRaster(shrub, ext)
# human      <- extendRaster(human, ext)
# dista_sqrt <- extendRaster(dista_sqrt, ext)

################################################################################
#### Combine Simulations
################################################################################
# Once the simulations are done, let's combine them into a single big dataframe
sims <- lapply(1:nrow(design), function(x) {
  dat            <- read_rds(design$Filename[x])
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
