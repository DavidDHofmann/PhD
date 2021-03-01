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
library(davidoff)   # Custom Functions
library(raster)     # For raster manipulation
library(tidyverse)  # For data wrangling
library(glmmTMB)    # For handling logistic model
library(rgdal)      # To handle spatial data
library(velox)      # For quick extraction
library(lubridate)  # To handle timestamps
library(tictoc)     # For simple benchmarking

# How many steps do you want to simulate?
n_steps <- 100

# How many random steps do you want to simulate per realized step?
n_rsteps <- 25

# How many dispersers do you want to simulate per source point?
n_disp <- 1

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

# Visualize all layers
plot(layers)

################################################################################
#### Sampling Source Points
################################################################################
# Load protected areas
prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Plot them nicely
u <- unique(prot$Desig)
m <- match(prot$Desig, u)
n <- length(unique(prot$Desig))
pal <- c("#E5F5E0", "#A1D99B", "#31A354")
plot(prot, col = pal[m])
  legend("topleft", legend = u, col = pal, pch = 19)

# Subset to national parks only, we'll use those as source areas
source_areas <- subset(prot, Desig == "National Park")
source_areas$AreaID <- 1:nrow(source_areas)

# Sample points within each national park
source_points <- lapply(1:nrow(source_areas), function(x){
  spsample(source_areas[x, ], n = 10, type = "random")
}) %>% do.call(rbind, .)
source_points$AreaID <- over(source_points, source_areas)$AreaID
source_points$PointID <- 1:nrow(source_points)

# Visualize the results and visually verify that IDs were assigned correctly
plot(source_areas, border = source_areas$AreaID)
  plot(source_points, col = source_points$AreaID, add = T)
  text(x = source_points, labels = source_points$PointID, cex = 0.6)

# Store the source points and areas
writeOGR(source_areas
  , dsn       = "03_Data/03_Results"
  , layer     = "99_SourceAreas"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
writeOGR(source_points
  , dsn       = "03_Data/03_Results"
  , layer     = "99_SourcePoints"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

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

################################################################################
#### TESTING
################################################################################
# Function to simulate dispersal based on a step selection model that was fitted
# in the glmmTMB framework
disperse <- function(
    source              = NULL    # Start Coordinates
  , covars              = NULL    # Spatial Covariates, prepared with our funct.
  , model               = NULL    # iSSF Model, prepared with our funct.
  , sl_dist             = NULL    # Step Length Distribution
  , sl_max              = Inf     # What is the largest possible step?
  , date                = as.POSIXct("2015-06-15 03:00:00", tz = "UTC")
  , n_steps             = 10      # Number of steps simulated
  , n_rsteps            = 25      # Number of random steps proposed
  , scaling             = NULL    # Dataframe to scale extracted covariates
  , stop                = F){     # Should the simulation stop at a boundary?

  # Create a new dataframe indicating the first location. Note that we draw
  # random turning angles to start off
  track <- data.frame(
      x           = coordinates(source)[, 1]
    , y           = coordinates(source)[, 2]
    , absta_      = runif(1, min = 0, max = 2 * pi)
    , ta_         = runif(1, min = -pi, max = pi)
    , sl_         = NA
    , Timestamp   = date
    , BoundaryHit = FALSE
  )

  # Simulate random steps
  for (i in 1:n_steps){

    # Draw random turning angles
    ta_new <- runif(n_rsteps
      , min = -pi
      , max = +pi
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

    # Check if the timestamp corresponds to low or high activity
    inactive <- strftime(date, tz = "UTC", format = "%H:%M:%S")
    inactive <- ifelse(inactive %in% c("03:00:00", "15:00:00"), T, F)
    inactive <- factor(inactive, levels = c(T, F))

    # Put all covariates into a dataframe. We will use this to calculate
    # selection scores
    covariates <- data.frame(
        extracted
      , cos_ta_   = cos(ta_new)
      , log_sl_   = log(sl_new)
      , sl_       = sl_new
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

    # Add the new date to the dataframe
    rand$Timestamp <- date

    # Coerce selection scores to probabilities
    Probs <- SelectionScore / sum(SelectionScore)

    # Sample a step according to the above predicted probabilities
    rand <- rand[sample(1:nrow(rand), size = 1, prob = Probs), ]

    # Add the step to our track
    track <- rbind(
        track[, c("x", "y", "absta_", "ta_", "sl_", "Timestamp", "BoundaryHit")]
      , rand[, c("x", "y", "absta_", "ta_", "sl_", "Timestamp", "BoundaryHit")]
    )
  }
  return(track)
}

################################################################################
#### TESTING
################################################################################
# Try out the function
# tic()
# test <- mclapply(1:nrow(points), mc.cores = detectCores() - 1, function(x){
#   disperse(
#       source    = points[x, ]
#     , covars    = covars
#     , model     = model
#     , sl_dist   = gamma
#     , sl_max    = sl_max
#     , n_rsteps  = n_rsteps
#     , n_steps   = n_steps
#     , stop      = F
#     , scaling   = scaling
#     , date      = as.POSIXct("2015-06-15 03:00:00", tz = "UTC")
#   )
# })
# toc()

tic()
test <- disperse_old(
    source    = source_points[1, ]
  , covars    = covars
  , model     = model
  , sl_dist   = sl_dist
  , sl_max    = sl_max
  , n_rsteps  = n_rsteps
  , n_steps   = 100
  , stop      = F
  , scaling   = scaling
  , date      = as.POSIXct("2015-06-15 03:00:00", tz = "UTC")
)
toc()

tic()
test <- disperse_new(
    source    = source_points[1, ]
  , covars    = covars
  , model     = model
  , sl_dist   = sl_dist
  , sl_max    = sl_max
  , n_rsteps  = n_rsteps
  , n_steps   = 100
  , stop      = F
  , scaling   = scaling
  , date      = as.POSIXct("2015-06-15 03:00:00", tz = "UTC")
)
toc()

# ################################################################################
# #### TESTING
# ################################################################################
# plot(source_areas, border = source_areas$ID)
# plot(source_points, col = source_points$ID, add = T)
# text(x = source_points, labels = source_points$ID, cex = 1)
#
# # Try out the function
# test <- list()
# for (j in 1:20){
#   i <- sample(1:68, 1)
#   random <- disperse(
#       source    = points[i, ]
#     , covars    = covars
#     , model     = model
#     , sl_dist   = gamma
#     , n_rsteps  = n_rsteps
#     , n_steps   = 100
#     , stop      = F
#     , date      = as.POSIXct("2015-06-15 03:00:00", tz = "UTC")
#   )
#   test[[j]] <- random
# }
# dat <- do.call(rbind, test)
#
# # Load true dispersal data
# observed <- read_csv("03_Data/02_CleanData/00_General_Dispersers_Popecol(SSF4Hours).csv")
# observed <- subset(observed, case_)
#
# # Compare step lengths
# dat <- data.frame(
#     sl = c(random$sl_, observed$sl_)
#   , type = c(rep("random", nrow(random)), rep("observed", nrow(observed)))
# )
#
# # Compare data
# tapply(dat$sl, dat$type, summary)
#
# # Visualize
# ggplot(dat, aes(x = sl, col = type)) + geom_density()
#
# ################################################################################
# #### TESTING
# ################################################################################

# We will have to run the simulations using two different approaches. First, we
# will select static source points. In a second run, source points will be
# randomized. Let's therefore write a wrapper function that allows to achieve
# this.
simulateDispersal <- function(randomize = F, filename = NULL){

  # Sample start points
  points <- createPoints(
      points    = source_points
    , areas     = source_areas
    , n         = n_disp
    , randomize = randomize
  )

  # Run the simulation according to the preliminary inputs for each point
  tracks <- pbmclapply(
      X                   = 1:length(points)
    , mc.cores            = detectCores() - 1
    , ignore.interactive  = T
    , FUN                 = function(x){
      sim <- suppressWarnings(disperse(
          source    = points[x, ]
        , covars    = covars
        , model     = model
        , sl_dist   = gamma
        , n_rsteps  = n_rsteps
        , n_steps   = n_steps
        , stop      = F
      ))

      # Assign ID that corresponds to the iteration number
      sim$ID <- x

      # Assign Point ID from which the disperser originated
      sim$StartPoint <- points$ID[x]

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

# Prepare directories into which we can store the results
dir.create("03_Data/03_Results/99_Simulations"
  , showWarnings = F
)
dir.create("03_Data/03_Results/99_Simulations/BigComputer"
  , showWarnings = F
)
dir.create("03_Data/03_Results/99_Simulations/BigComputer/Static"
  , showWarnings = F
)
dir.create("03_Data/03_Results/99_Simulations/BigComputer/Random"
  , showWarnings = F
)

################################################################################
#### Run Simulation
################################################################################
# Run each simulation 10 times
tic()
for (i in 1:1){
  simulateDispersal(
      randomize = F
    , filename  = paste0(
        "03_Data/03_Results/99_Simulations/Static/Iteration_"
      , i
      , ".rds"
    )
  )
  simulateDispersal(
      randomize = T
    , filename  = paste0(
        "03_Data/03_Results/99_Simulations/Random/Iteration_"
      , i
      , ".rds"
    )
  )
}
toc()
