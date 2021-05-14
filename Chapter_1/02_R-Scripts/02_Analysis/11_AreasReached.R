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

################################################################################
#### Rasterize National Parks
################################################################################
# Load protected areas
prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Load reference raster
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")

# Assign a unique ID to each protected area
prot$ID <- 1:nrow(prot)

# Rasterize IDs
prot_r <- raster(terra::rasterize(x = vect(prot), y = rast(r), field = "ID"))

# Visualize them
plot(prot_r, main = "National Parks")
text(subset(prot, Desig == "National Park"), "Name", cex = 0.5)

################################################################################
#### Identify Connections between Protected Areas
################################################################################
# Load simulations
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")

# Create SpatialPoints for first location of each trajectory
first <- sims %>%
  dplyr::select("x", "y", "TrackID") %>%
  group_by(TrackID) %>%
  slice(1) %>%
  SpatialPointsDataFrame(
      coords      = cbind(.[["x"]], .[["y"]])
    , proj4string = CRS("+init=epsg:4326")
  )

# Assess from which protected area each trajectory leaves
first$From <- raster::extract(prot_r, first)

# Join information to simulated tracks
sims <- left_join(sims, first@data[, c("TrackID", "From")], by = "TrackID")

# Remove potential NA's (start points outside the main study area)
sims <- subset(sims, !is.na(From))

# Identify number of trajectories leaving from each area
nsims <- sims %>%
  group_by(TrackID, From) %>%
  nest() %>%
  ungroup() %>%
  count(From) %>%
  arrange(From) %>%
  setNames(c("From", "Simulations"))

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
sims$To <- visits$Prot

# Coerce to regular dataframe again
sims <- as.data.frame(sims, xy = T)
sims$xy <- NULL

# Identify how long it takes on average to reach different areas
visits <- sims %>%
  group_by(TrackID, From, To) %>%
  summarize(StepNumber = min(StepNumber), .groups = "drop") %>%
  subset(!is.na(From) & !is.nan(To)) %>%
  arrange(TrackID, StepNumber) %>%
  group_by(From, To) %>%
  summarize(
      MeanStepNumber = mean(StepNumber)
    , SDStepNumber   = sd(StepNumber)
    , Frequency      = n()
    , .groups        = "drop"
  ) %>%
  subset(From != To)

# Join the number of dispersers
visits <- left_join(visits, nsims, by = "From")

# Calculate relative frequency
visits$RelFrequency <- visits$Frequency / visits$Simulations

# We also want to know the designation and name of the origin...
visits <- visits %>%
  dplyr::select(From) %>%
  left_join(prot@data[, c("Name", "Desig", "ID")], by = c("From" = "ID")) %>%
  dplyr::select(Name, Desig) %>%
  setNames(c("FromName", "FromDesig")) %>%
  cbind(visits, .)

# And the designation and name of the destination
visits <- visits %>%
  dplyr::select(To) %>%
  left_join(prot@data[, c("Name", "Desig", "ID")], by = c("To" = "ID")) %>%
  dplyr::select(Name, Desig) %>%
  setNames(c("ToName", "ToDesig")) %>%
  cbind(visits, .)

# Sort all nicely
visits <- visits %>%
  dplyr::select(From, To, FromName, ToName, FromDesig, ToDesig, everything())

# Store areas reached to file
write_rds(visits, "03_Data/03_Results/99_AreasReached.rds")
