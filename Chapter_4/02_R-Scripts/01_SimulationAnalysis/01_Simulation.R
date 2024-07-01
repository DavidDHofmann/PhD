################################################################################
#### Simulating Movement Data with Known Preferences
################################################################################
# Description: Simulate movement with known preferences on virtual landscape

# Clear R's brain
rm(list = ls())

# Load required packages
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(raster)        # To handle spatial data
library(sf)            # For plotting spatial features
library(lubridate)     # To handle dates and times

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_4"
setwd(wd)

# Load custom functions
source("02_R-Scripts/Functions.R")

################################################################################
#### Simulation Parameters
################################################################################
n              <- 300                                                          # Resolution of the covariate layers
n_rsteps       <- 1000                                                         # Number of random steps to be generated
n_steps        <- 1000                                                         # Number of consecutive steps to be simulated
formula        <- ~ dist + elev + forest                                       # Formula used to predict step-selection scores
prefs          <- c(-20, 0.5, -1)                                              # Preferences of the simulated individuals
sl_dist        <- list(name = "gamma", params = list(shape = 3, scale = 1))    # Step-length parameters
ta_dist        <- list(name = "vonmises", params = list(kappa = 0.5, mu = 0))  # Turning-angle parameters
stop           <- F                                                            # Should the simulation terminate at boundaries?
autocorr_range <- c(5, 10, 20)                                                 # Autocorrelation ranges to generate landscapes
n_replicates   <- 100                                                          # Number of replications
extract_where  <- "end"                                                        # Whether to extract covariates at the end or along steps

# Let's specify the extent on which animals are allowed to move
ext <- extent(0, n, 0, n)
ext <- as(ext, "SpatialPolygons")

# Let's also specify an extent within which individuals will be released
ext2 <- extent(ext) - 100
ext2 <- as(ext2, "SpatialPolygons")

# And an extent to which we will expand covariates before extracting covariate
# values
ext3 <- extent(c(-50, 350, -50, 350))
ext3 <- as(ext3, "SpatialPolygons")

# Generate a point of attraction at the center of the study area and compute the
# distance to it
r <- raster(nrows = n, ncols = n, xmn = 0, xmx = n, ymn = 0, ymx = n)
cent <- SpatialPoints(t(c(n / 2, n / 2)))
dist <- distanceFromPoints(r, cent)
dist <- (dist - cellStats(dist, min)) / (cellStats(dist, max) - cellStats(dist, min))

################################################################################
#### Main Simulation Functions
################################################################################
# Function to simulate the different covariates and return them as a stack
simCovars <- function(autocorr_range = 10, proportion_forest = 0.5, seed) {

  # Simulate distance layer
  dist <- simDistance(n = n
    , x = n / 2
    , y = n / 2
  )

  # Simulate elevation layer
  elev <- simElevation(n = n
    , autocorr_range = autocorr_range
    , seed           = seed
  )

  # Simulate forest layer
  forest <- simForest(n = n
    , autocorr_range = autocorr_range
    , prop           = proportion_forest
    , seed           = -seed
  )

  # Put them into a stack
  covars <- stack(dist, elev, forest)
  names(covars) <- c("dist", "elev", "forest")
  return(covars)

}

# Function to simulate movement on a simulated landscape. Returns "observed" GPS
# data
simMove <- function(covars, messages = F) {

  # Simulate movement
  sim <- move(
      xy       = matrix(runif(2, xmin(ext2), xmax(ext2)), ncol = 2)
    , covars   = covars
    , formula  = formula
    , prefs    = prefs
    , sl_dist  = sl_dist
    , ta_dist  = ta_dist
    , ext      = ext
    , n_steps  = n_steps
    , n_rsteps = n_rsteps
    , stop     = stop
    , messages = messages
    , extract  = extract_where
  )

  # Add an artifial timestamp, the ID, as well as a unique step_id to the
  # data
  sim$timestamp <- lubridate::ymd_hms("2000-01-01 00:00:00") + hours(1:nrow(sim))

  # In reality, we don't observe step lengths etc.
  sim$absta <- NULL
  sim$ta    <- NULL
  sim$sl    <- NULL
  return(sim)
}

################################################################################
#### Ensure All Functions Work as Intended
################################################################################
# Simulate covariates and across the layers
cov <- simCovars(autocorr_range = 10, proportion_forest = 0.5, seed = 1)
obs <- simMove(cov, messages = T)

# Visualize
as.data.frame(cov, xy = T) %>%
  gather(key = covariate, value = value, 3:5) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = st_as_sf(ext2), inherit.aes = F, fill = NA, col = "white", lty = 2) +
    geom_path(data = obs, aes(x = x, y = y), inherit.aes = F) +
    scale_fill_viridis_c(option = "viridis") +
    coord_sf() +
    theme_minimal() +
    facet_wrap("covariate") +
    theme(
        axis.title.y    = element_text(angle = 0, vjust = 0.5)
      , legend.position = "none"
    )

################################################################################
#### Simulation
################################################################################
# Specify the different design combinations through which we want to run
dat <- expand_grid(
    AutocorrRange = autocorr_range
  , Replicate     = 1:n_replicates
)

# Prepare filenames for the simulated data
dat <- mutate(dat
  , CovariateFilename  = paste0("03_Data/Simulations/Covariates_A", sprintf("%03d", AutocorrRange), "_R", sprintf("%03d", Replicate), ".grd")
  , MovementFilename   = paste0("03_Data/Simulations/Movement_A", sprintf("%03d", AutocorrRange), "_R", sprintf("%03d", Replicate), ".rds")
)

# Create the respective directory
dir.create("03_Data/Simulations", showWarnings = F, recursive = T)

# Simulate the covariate layers
cat("Simulating covariate layers...\n")
dat$Covariates <- pbmclapply(
    X                  = 1:nrow(dat)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {
    filename <- dat$CovariateFilename[x]
    if (file.exists(filename)) {
      cov <- stack(filename)
    } else {
      cov <- simCovars(autocorr_range = dat$AutocorrRange[x], seed = dat$Replicate[x])
      cov <- writeRaster(cov, filename)
    }
    return(cov)
})

# Now simulate movement across the covariates
cat("Simulating movement...\n")
dat$Movement <- pbmclapply(
    X                  = 1:nrow(dat)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {
    filename <- dat$MovementFilename[x]
    if (file.exists(filename)) {
      sim <- read_rds(filename)
    } else {
      sim         <- simMove(covars = dat$Covariates[[x]], messages = F)
      sim$id      <- x
      sim$step_id <- (1:n_steps) + ((x - 1) * n_steps)
      write_rds(sim, filename)
    }
    return(sim)
})

# Let's expand the covariates slightly so that later we won't have random
# steps leaving the study area
cat("Extending covariate layers...\n")
dat$Covariates <- pbmclapply(
    X                  = 1:nrow(dat)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {
    file <- dat$CovariateFilename[x]
    cov  <- stack(file)
    if (extent(cov) != extent(ext3)) {
      cov <- extendRaster(cov, ext3)
      cov <- writeRaster(cov, file, overwrite = T)
    }
    return(cov)
})

# Ensure that all covariate layers are loaded propery
dat$Covariates <- lapply(dat$CovariateFilename, stack)

# Store the entire simulation to file
write_rds(dat, "03_Data/Simulation.rds")

# Store session information
writeLines(capture.output(sessionInfo()), "02_R-Scripts/01_SimulationAnalysis/SessionInfo/01_SimulationSession.txt")
