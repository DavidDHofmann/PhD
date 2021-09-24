################################################################################
#### Number of Connections Between Protected Areas
################################################################################
# Description: In this script, we analyze the simulated dispersal trajectories
# and identify the number of connections between protected areas

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(pbmcapply)      # To run on multiple cores with progress bar
library(rgeos)          # For spatial ananylsis
library(davidoff)       # Custom functions
library(rgdal)          # To load spatial data
library(sf)             # To plot spatial stuff
library(igraph)         # To plot networks

################################################################################
#### Rasterize National Parks
################################################################################
# Load protected areas and subset to national parks only
prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
prot <- subset(prot, Desig == "National Park")

# Load reference raster and rasterize all national parks
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")
prot$ID <- 1:nrow(prot)
prot_r <- raster(terra::rasterize(x = vect(prot), y = rast(r), field = "ID"))

# Visualize them
plot(prot_r, main = "National Parks", horizontal = T)
plot(prot, add = T)
text(prot, "Name", cex = 0.5, halo = T)

################################################################################
#### Prepare Simulations
################################################################################
# Load simulations. We'll focus on simulations initiated in the main study area
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")
sims <- subset(sims, Area == "Main")
sims <- dplyr::select(sims, x, y, TrackID, StepNumber)

# Identify the origin (source area) of each trajectory. For this, we create
# SpatialPoints from the first location of each trajectory
first <- sims %>%
  dplyr::select("x", "y", "TrackID") %>%
  group_by(TrackID) %>%
  slice(1) %>%
  SpatialPointsDataFrame(
      coords      = cbind(.[["x"]], .[["y"]])
    , proj4string = CRS("+init=epsg:4326")
  )

# Assess from which protected area each trajectory leaves
first$SourceArea <- raster::extract(prot_r, first)

# Join information to simulated tracks
sims$SourceArea <- first$SourceArea[match(sims$TrackID, first$TrackID)]

# Remove NA's (start points outside national parks)
sims <- subset(sims, !is.na(SourceArea) & !is.nan(SourceArea))

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# Identify through which national parks the dispersers moved
visits <- data.frame(
    TrackID    = sims$TrackID
  , StepNumber = sims$StepNumber
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , Prot       = raster::extract(prot_r, sims)
)

# Add this information to the simulations
sims$CurrentPark <- visits$Prot

# Coerce to regular dataframe again
sims <- as.data.frame(sims, xy = T)
sims$xy <- NULL

# Each national park belongs to a country, let's assign the repsective country
# to the dataframe as well
sims$SourceAreaCountry <- as.character(prot$Country[match(sims$SourceArea, prot$ID)])
sims$CurrentParkCountry <- as.character(prot$Country[match(sims$CurrentPark, prot$ID)])

################################################################################
#### Identify Direct Connections between National Parks
################################################################################
# Note: Here we are going to look at all direct connections between national
# parks. That is, if a trajectory moves from A through B to C, we generate an
# edge between A and B and between A and C, but not between B and C.

# Identify number of trajectories leaving from each area
nsims <- sims %>%
  group_by(TrackID, SourceArea) %>%
  nest() %>%
  ungroup() %>%
  count(SourceArea) %>%
  arrange(SourceArea) %>%
  setNames(c("From", "Simulations"))

# Add the park and country from which they left
nsims$FromPark    <- prot$Name[match(nsims$From, prot$ID)]
nsims$FromCountry <- prot$Country[match(nsims$From, prot$ID)]

# Identify how long it takes on average to reach different areas
visits <- sims %>%
  rename(From = SourceArea, To = CurrentPark) %>%
  group_by(TrackID, From, To) %>%
  summarize(StepNumber = min(StepNumber), .groups = "drop") %>%
  subset(!is.na(From) & !is.nan(To) & !is.na(To)) %>%
  arrange(TrackID, StepNumber)

# Get the source park name, destination park name, source country, destination
# country
visits$FromPark    <- prot$Name[match(visits$From, prot$ID)]
visits$ToPark      <- prot$Name[match(visits$To, prot$ID)]
visits$FromCountry <- prot$Country[match(visits$From, prot$ID)]
visits$ToCountry   <- prot$Country[match(visits$To, prot$ID)]

# Remove factors
visits$FromPark    <- as.character(visits$FromPark)
visits$ToPark      <- as.character(visits$ToPark)
visits$FromCountry <- as.character(visits$FromCountry)
visits$ToCountry   <- as.character(visits$ToCountry)

# Store visits to file
write_rds(visits, "03_Data/03_Results/99_InterpatchConnectivity.rds")
