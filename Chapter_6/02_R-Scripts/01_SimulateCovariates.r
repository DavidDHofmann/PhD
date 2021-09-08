################################################################################
#### Simulating Covariates for the Movement Simulation
################################################################################
# Description: Generate artificial landscapes through which the simulated
# individuals will move. We will replicate each simulation 100 times

# Clear R's brain
rm(list = ls())

# Load required packages
library(NLMR)          # To simulate covariates
library(raster)        # To handle spatial data
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(sf)            # For plotting spatial features
library(rgeos)         # For manipulating spatial objects
library(ggpubr)        # To put plots together

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_6")

################################################################################
#### Simulate Covariates
################################################################################
# Set seed
set.seed(12345)

# Specify resolution of covariates
n <- 500

# Function to simulate water cover
simWater <- function(n){
  water <- nlm_mosaictess(
      ncol  = n
    , nrow  = n
    , germs = 500
  )
  rcl <- cbind(
      old = unique(values(water))
    , new = rbinom(n = length(unique(values(water))), size = 1, prob = 0.3)
  )
  water <- reclassify(water, rcl)
  return(water)
}

# Function to simulate elevation
simElev <- function(n){
  elev <- nlm_gaussianfield(
      ncol           = n
    , nrow           = n
    , autocorr_range = 25
    , mag_var        = 1
    , nug            = 0
  )
  return(elev)
}

# Try the simulation functions
water <- simWater(n)
elev  <- simElev(n)
plot(stack(water, elev))

# Distribute two source areas
nps <- SpatialPoints(rbind(c(150, 250), c(350, 250)))
nps <- gBuffer(nps, byid = T, width = 20)
nps$ID <- 1:length(nps)

# Calculate distance to points of interest and normalize the distance afterwards
dist <- distanceFromPoints(elev, gCentroid(nps, byid = T))
dist <- (dist - cellStats(dist, min)) / (cellStats(dist, max) - cellStats(dist, min))

# Specify the extent of the study area
ext <- as(extent(dist), "SpatialPolygons")

# Specify the core study area
core <- as(extent(100, 400, 100, 400), "SpatialPolygons")

# Calculate buffer from this
buffer <- ext - core

# Put different polygons together
pols <- rbind(ext, core, buffer, makeUniqueIDs = T)
pols$Name <- c("Extent", "Core", "Buffer")

# Let's replicate the layer generation x times
covars_replicated <- mclapply(1:100, mc.cores = detectCores() - 1, function(x){
  water         <- simWater(n)
  elev          <- simElev(n)
  covars        <- stack(water, elev, dist)
  names(covars) <- c("water", "elev", "dist")
  return(covars)
})

################################################################################
#### Generate Buffer Zone
################################################################################
# Create covariate layers where buffer zones are filled by random values from
# the true layer
covars_replicated_randomized <- mclapply(covars_replicated, mc.cores = detectCores() - 1, function(y){
  randomized <- lapply(1:nlayers(y), function(x){
    values_core <- y[[x]][core][[1]]
    r <- raster(y)
    r <- setValues(r, sample(values_core, size = ncell(r), replace = T))
    r <- mask(r, core, inverse = T)
    r <- cover(r, y[[x]])
    return(r)
  })
  randomized <- stack(randomized)
  names(randomized) <- names(y)
  return(randomized)
})

# Finally, we want to have a set of cropped covariates
covars_replicated_cropped <- mclapply(covars_replicated, mc.cores = detectCores() - 1, function(y){
  crop(y, core)
})

# Plot the real covariates
plot_real <- as.data.frame(covars_replicated[[1]], xy = T) %>%
  gather(key = covariate, value = value, 3:5) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = st_as_sf(nps), inherit.aes = F, fill = NA, col = "white") +
    geom_sf(data = st_as_sf(core), inherit.aes = F, fill = NA, col = "red", lty = 2) +
    geom_sf(data = st_as_sf(buffer), inherit.aes = F, fill = "white", col = NA, alpha = 0.3) +
    scale_fill_viridis_c(option = "viridis") +
    coord_sf() +
    theme_minimal() +
    facet_wrap("covariate") +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

# Plot the covariates with randomized buffer
plot_randomized <- as.data.frame(covars_replicated_randomized[[1]], xy = T) %>%
  gather(key = covariate, value = value, 3:5) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = st_as_sf(nps), inherit.aes = F, fill = NA, col = "white") +
    geom_sf(data = st_as_sf(core), inherit.aes = F, fill = NA, col = "red", lty = 2) +
    geom_sf(data = st_as_sf(buffer), inherit.aes = F, fill = "white", col = NA, alpha = 0.3) +
    scale_fill_viridis_c(option = "viridis") +
    coord_sf() +
    theme_minimal() +
    facet_wrap("covariate") +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

# Plot the cropped covariates
plot_cropped <- as.data.frame(covars_replicated_cropped[[1]], xy = T) %>%
  gather(key = covariate, value = value, 3:5) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = st_as_sf(nps), inherit.aes = F, fill = NA, col = "white") +
    geom_sf(data = st_as_sf(core), inherit.aes = F, fill = NA, col = "red", lty = 2) +
    geom_sf(data = st_as_sf(buffer), inherit.aes = F, fill = "white", col = NA, alpha = 0.3) +
    scale_fill_viridis_c(option = "viridis") +
    coord_sf() +
    theme_minimal() +
    facet_wrap("covariate") +
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

# Put plots together
p <- ggarrange(plot_real, plot_randomized, plot_cropped, ncol = 1)
p

# Store the plot
ggsave(plot = p, "04_Manuscript/99_CovariatePlots.png", width = 8, height = 8)

# Put all covariates together
covars <- tibble(
    Replicate  = 1:length(covars_replicated)
  , Real       = covars_replicated
  , Randomized = covars_replicated_randomized
  , Cropped    = covars_replicated_cropped
)

# Gather them
covars <- gather(covars, key = Type, value = Layers, 2:ncol(covars))

# Look at the final object
print(covars)
format(object.size(covars), "auto")

# Store simulated data to file
save(covars, nps, ext, core, buffer, file = "03_Data/Landscape.Rdata")
