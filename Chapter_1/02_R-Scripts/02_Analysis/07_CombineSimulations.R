################################################################################
#### Combine Simulations into Single Object
################################################################################
# Description: Take all simulations and put them into a single r file

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)    # For data wrangling
library(raster)       # For handling spatial data
library(rgdal)        # For handling spatial data
library(pbmcapply)    # For multicore processes with progress bar
library(davidoff)     # Access to custom functions
library(rgeos)        # For spatial data manipulation
library(lemon)        # For capped coordinates

################################################################################
#### Putting Data Together
################################################################################
# Identify all simulation files
files <- dir(
    "03_Data/03_Results/99_Simulations/Main"
  , pattern     = ".rds$"
  , full.names  = T
)

# How many are there?
length(files)

# Load the files and bind their rows together
sims <- lapply(1:length(files), function(x){
  dat <- read_rds(files[x])
  dat$SimID <- x
  return(dat)
}) %>% do.call(rbind, .)

# Take a look at the data
head(sims)

# Count the number of simulated steps (in Mio)
nrow(sims) / 1e6

# Because in each simulation we start off with new IDs for trajectories, they
# are not unique across simulations. We thus combine ID and SimID to create an
# ID that is unqiue to each simulated path, across all simulations
sims <- sims %>% mutate(TrackID = group_indices(., SimID, TrackID))

# Make sure it worked
table(table(sims$TrackID))

# Collect garbage
gc()

# Let's also create a step counter, indicating the number of the step in its
# respective trajectory
sims <- sims %>% group_by(TrackID) %>% mutate(StepNumber = (row_number()))

# Check object size
format(object.size(sims), units = "Gb")

# Write to an rds
write_rds(sims, "03_Data/03_Results/99_DispersalSimulationSub.rds")

################################################################################
#### Test Convergence
################################################################################
# Coerce simulations to proper trajectories
tracks <- sims2tracks(sims, keep.data = T, steps = 2000)

# Load reference shapefile
s <- readOGR("03_Data/02_CleanData/00_General_Shapefile.shp")

# Distribute 1000 "chek points" in our study area. For this, generate 5km2
# squares in our study area
check <- raster(s, res = metersToDegrees(5000))
values(check) <- 0
values(check)[sample(1:ncell(check), 1000)] <- 1
check <- rasterToPolygons(check, fun = function(x){x == 1})

# Assign a checkpoint ID to each polygon
check$CheckID <- 1:length(check)

# Extract all trajectory ids
ids <- unique(tracks$TrackID)
length(ids)

# Load source areas
source_areas <- readOGR("03_Data/03_Results/99_SourceAreas.shp")
buffer_areas <- readOGR("03_Data/03_Results/99_BufferArea.shp")

# Visualize them
plot(buffer_areas
  , axes   = T
  , las    = 1
  , col    = "gray40"
  , border = NA
  , main   = "Checkpoints in Study Area"
)
plot(source_areas
  , col    = colTrans("darkgreen")
  , border = NA
  , main   = "Checkpoints in Study Area"
  , add    = T
)
plot(check
  , cex    = 0.3
  , pch    = 20
  , col    = "blue"
  , border = NA
  , add    = T
)
plot(tracks[sample(nrow(tracks), 10), ]
  , add = T
  , col = "gold"
  , lwd = 0.6
)
legend("bottomright"
  , legend = c("Checkpoints", "Source Areas", "Buffer Zone", "Simulated Dispersers")
  , pch    = c(15, 15, 15, NA)
  , lty    = c(NA, NA, NA, 1)
  , col    = c("blue", colTrans("darkgreen"), "gray", "gold")
)

# Prepare a design matrix through which we will loop
design <- expand_grid(
    NTracks    = seq(0, length(unique(tracks$TrackID)), by = 500)
  , Replicate  = c(1:100)
)

# Let's shuffle the matrix (makes prediction of calculation time in pbmclapply
# more accurate)
design <- design[sample(nrow(design), replace = F), ]

# Go through the design matrix and assess the traversal frequency at the
# checkpoints
convergence <- mutate(design, Frequency = pbmclapply(
    X                  = 1:nrow(design)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x){

  # Sample simulated trajectories. Some lines may be sampled twice, which leads
  # to an issue with non-unique lineIDs. Hence, I'll change them if needed
  sub <- tracks[sample(length(tracks), size = design$NTracks[x], replace = T), ]
  if (length(sub) > 0){
    quiet <- lapply(1:length(sub), function(z){
      sub@lines[[z]]@ID <<- as.character(z)
    })
  }

  # Check number of intersections at each checkpoint
  ints <- suppressWarnings(gIntersects(sub, check, byid = T, prepared = T))
  ints <- rowSums(ints)
  ints <- enframe(ints, name = "CheckID", value = "Traversals")

  # Remove garbage
  gc()

  # Return the results
  return(ints)

}))

# Check object size
format(object.size(convergence), units = "Gb")

# Store results
write_rds(convergence, "03_Data/03_Results/99_Convergence.rds")
convergence <- read_rds("03_Data/03_Results/99_Convergence.rds")

# Clean results and calculate relative traversal frequency
convergence <- convergence %>%
  unnest(cols = Frequency) %>%
  mutate(RelativeTraversals = Traversals / NTracks) %>%
  mutate(RelativeTraversals = ifelse(is.nan(RelativeTraversals)
  , 0
  , RelativeTraversals)
)

# Look for local convergence
convergence %>%
  subset(CheckID %in% sample(unique(convergence$CheckID), 8)) %>%
  group_by(NTracks, CheckID) %>%
  summarize(
      MeanRelativeTraversals = mean(RelativeTraversals)
    , Upper = quantile(RelativeTraversals, 0.975)
    , Lower = quantile(RelativeTraversals, 0.025)
  ) %>%
  # subset(MeanRelativeTraversals > 0) %>%
  ggplot(aes(x = NTracks, y = MeanRelativeTraversals)) +
  geom_ribbon(aes(
      ymin = Lower
    , ymax = Upper
  ), alpha = 0.5, fill = "orange", color = "orange") +
  geom_line() +
  geom_point(size = 0.2) +
  facet_wrap("CheckID", scales = "free", nrow = 4) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  )

# Look for global convergence
convergence %>%
  group_by(NTracks, Replicate) %>%
  summarize(
      RelativeTraversals = mean(RelativeTraversals)
  ) %>%
  group_by(NTracks) %>%
  summarize(
      MeanRelativeTraversals = mean(RelativeTraversals)
    , Upper = quantile(RelativeTraversals, 0.975)
    , Lower = quantile(RelativeTraversals, 0.025)
  ) %>%
  subset(MeanRelativeTraversals > 0) %>%
  ggplot(aes(x = NTracks, y = MeanRelativeTraversals)) +
  geom_ribbon(aes(
      ymin = Lower
    , ymax = Upper
  ), alpha = 0.5, fill = "orange", color = "orange") +
  geom_line() +
  geom_point() +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  )
