################################################################################
#### Combine Simulations into Single Object
################################################################################
# Description: Take all simulations and put them into a single r file

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/media/david/My Passport/Backups/WildDogs/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(tidyverse)    # For data wrangling
library(raster)       # For handling spatial data
library(rgdal)        # For handling spatial data

################################################################################
#### Putting Data Together
################################################################################
# Identify all simulation files
dir1 <- dir(
    "03_Data/03_Results/99_Simulations/BigComputer/Static"
  , pattern     = ".rds$"
  , full.names  = T
)
dir2 <- dir(
    "03_Data/03_Results/99_Simulations/BigComputer/Random"
  , pattern     = ".rds$"
  , full.names  = T
)

# Put them together
files <- rbind(
    data.frame(filename = dir1, StartPoints = rep("Static", length(dir1)))
  , data.frame(filename = dir2, StartPoints = rep("Random", length(dir2)))
)

# Let's check the files
files

# Load the files and bind their rows together
sims <- lapply(1:nrow(files), function(x){

  # Read the respective file
  data <- read_rds(as.character(files$filename[x]))

  # Add an indicator of the simID
  data$SimID = x

  # Indicate if start points where sampled statically or randomly
  data$PointSampling = files$StartPoints[x]

  # Return the respective data
  return(data)

  # Bind everything together
}) %>% do.call(rbind, .)

# Free some memory
gc()

# Take a look at the data
head(sims)

# Count the number of simulated steps (in Mio)
nrow(sims) / 1e6

# Because in each simulation we start off with new IDs for trajectories, they
# are not unique across simulations. We thus combine ID and SimID to create an
# ID that is unqiue to each simulated path, across all simulations
sims <- sims %>% mutate(ID = group_indices(., SimID, ID))

# Make sure it worked
table(table(sims$ID))

# Collect garbage
gc()

# Let's also create a step counter, indicating the number of the step in its
# respective trajectory
sims <- sims %>% group_by(ID) %>% mutate(StepNumber = (row_number() - 1))

# Check object size
format(object.size(sims), units = "Gb")

# Write to an rds
write_rds(sims, "03_Data/03_Results/99_DispersalSimulationSub.rds")

# ################################################################################
# #### Create Shapefile
# ################################################################################
# # Maybe we want to visualize some of the trajectories in QGIS. Let's subset the
# # data and store it as shapefile
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationTest.rds")
#
# # Subset data
# sims <- sims %>%
#   group_by(StartPoint, StartPoints, ID) %>%
#   nest() %>%
#   group_by(StartPoints, StartPoint) %>%
#   sample_n(10) %>%
#   unnest()
#
# # Create trajectories
# sub <- sims %>% group_by(ID) %>% nest()
# sub_traj <- lapply(1:length(sub$data), function(x){
#   p <- SpatialPoints(sub$data[[x]][, c("x", "y")])
#   l <- spLines(p)
#   l$StartPoint <- sub$data[[x]]$StartPoint[1]
#   l$StartPoints <- sub$data[[x]]$StartPoints[1]
#   return(l)
# }) %>% do.call(rbind, .)
#
# # Visualize
# plot(sub_traj, col = sub_traj$StartPoints)
#
# # Store them to file
# writeOGR(sub_traj, ".", "test.shp", driver = "ESRI Shapefile")
