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
#### Approach I: Identify Direct Connections between National Parks
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

# Identify how long it takes on average to reach different areas
visits_direct <- sims %>%
  rename(From = SourceArea, To = CurrentPark) %>%
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
visits_direct <- left_join(visits_direct, nsims, by = "From")

# Calculate relative frequency
visits_direct$RelFrequency <- visits_direct$Frequency / visits_direct$Simulations

# Store areas reached to file
write_rds(visits_direct, "03_Data/03_Results/99_DirectInterpatchConnectivity.rds")

# Repeat the same analysis on a country level
nsims <- sims %>%
  group_by(TrackID, SourceAreaCountry) %>%
  nest() %>%
  ungroup() %>%
  count(SourceAreaCountry) %>%
  arrange(SourceAreaCountry) %>%
  setNames(c("From", "Simulations"))

# Identify how long it takes on average to reach different areas
visits_direct <- sims %>%
  rename(From = SourceAreaCountry, To = CurrentParkCountry) %>%
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
visits_direct <- left_join(visits_direct, nsims, by = "From")

# Calculate relative frequency
visits_direct$RelFrequency <- visits_direct$Frequency / visits_direct$Simulations

# Store areas reached to file
write_rds(visits_direct, "03_Data/03_Results/99_DirectInterpatchConnectivityCountries.rds")

# ################################################################################
# #### Approach II: Identify Indirect Connections between National Parks
# ################################################################################
# # Note: Here we are going to look at all indirect connections between national
# # parks. That is, if a trajectory moves from A through B to C, we generate an
# # edge between A and B and between B and C, but not between A and C.
#
# # Replace NaNs
# visits$Prot[is.nan(visits$Prot)] <- NA
#
# # Function to retrieve the visitation history from a sequence of values. Note
# # that we need to preserve self-loops so that we can correctly calculate the
# # number of outgoing trajectories from each national park
# visitHist <- function(x, singlecount = F){
#
#   # Indicate the step number
#   StepNumber <- 0:(length(x) - 1)
#
#   # Identify transition of each step (from A to B)
#   transitions <- data.frame(From = lag(x), To = x, StepNumber) %>% na.omit()
#
#   # Logical check if the current trajectory ever left its source area. If this
#   # is the case, we want to know the number of steps it took to move from A to
#   # B, from B to C and so on.
#   left <- length(unique(c(transitions$From, transitions$To))) > 1
#
#   # In case the trajectory ever left the source area, we can calculate the
#   # duration it took to move between national parks
#   if (left){
#     transitions_out <- transitions %>%
#       subset(From != To) %>%
#       group_by(From, To) %>%
#       summarize(Frequency = n(), StepNumber = min(StepNumber), .groups = "drop") %>%
#       arrange(StepNumber) %>%
#       mutate(StepNumber = StepNumber - lag(StepNumber, default = 0))
#   }
#
#   # In case the trajectory never left the source area, we simply add a duration
#   # of 0
#   transitions <- transitions %>%
#     subset(From == To) %>%
#     group_by(From, To) %>%
#     summarize(Frequency = n(), StepNumber = min(StepNumber), .groups = "drop") %>%
#     arrange(StepNumber) %>%
#     mutate(StepNumber = 0)
#   if (left){
#     transitions <- rbind(transitions, transitions_out)
#   }
#   if (singlecount){
#     transitions$Frequency = 1
#   }
#   return(transitions)
# }
#
# # Test it. Note that this version also returns the duration required for the
# # links
# visitHist(c(1, 2, 2, 4))
# visitHist(c(1, 1, 1, 1))
#
# # Nest by trajectory
# visits_indirect <- nest(visits, data = -TrackID)
#
# # Identify visitation history of each track
# visits_indirect$hist <- pbmclapply(
#     X                  = visits_indirect$data
#   , mc.cores           = detectCores() - 1
#   , ignore.interactive = T
#   , FUN                = function(x){
#
#     # Check if last visited area is an NA, if so, we need to keep that so that
#     # we can correctly calculate the number of outgoing trajectories at each
#     # national park
#     last_na <- is.na(x$Prot[length(x$Prot)])
#
#     # Identify last visit before entering NA (i.e. where the trajectory comes
#     # from before ending up in NA)
#     last_visit <- na.omit(x$Prot)[length(na.omit(x$Prot))]
#
#     # Get the visitation history
#     hist <- visitHist(na.omit(x$Prot), singlecount = T)
#
#     # If the last entry was an NA, add it back to the visits and indicate where
#     # the trajectory came from
#     if (last_na){
#       last <- data.frame(From = last_visit, To = NA, Frequency = 1, StepNumber = NA)
#       hist <- rbind(hist, last)
#     }
#
#     # Return the visits
#     return(hist)
# })
#
# # Unnest all visits and put them into a single dataframe
# hist <- visits_indirect %>%
#   dplyr::select(TrackID, hist) %>%
#   unnest(hist) %>%
#   dplyr::select(-Frequency) %>%
#   mutate(To = replace_na(To, 9999)) %>%
#   subset(From != To)
#
# # Count the number of outgoing trajectories at each national park. This will
# # include trajectories that originate form another national park but pass
# # through the respective national park!!!
# nsims <- hist %>%
#   count(From) %>%
#   setNames(c("From", "Simulations"))
#
# # Now we remove all NAs (replaced with 9999)
# hist <- subset(hist, To != 9999)
#
# # Summarize connections
# hist <- hist %>%
#   group_by(From, To) %>%
#   summarize(
#       MeanStepNumber = mean(StepNumber)
#     , SDStepNumber   = sd(StepNumber)
#     , Frequency      = n()
#     , .groups         = "drop"
#   )
#
# # Join the number of dispersers
# visits_indirect <- left_join(hist, nsims, by = "From")
#
# # Calculate relative frequency
# visits_indirect$RelFrequency <- visits_indirect$Frequency / visits_indirect$Simulations
#
# # Store to file
# write_rds(visits_indirect, "03_Data/03_Results/99_IndirectInterpatchConnectivity.rds")
