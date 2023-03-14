################################################################################
#### Interpatch Connectivity
################################################################################
# Description: Computing interpatch connectivity from simulated dispersal paths

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(pbmcapply)      # To run stuff in parallel
library(igraph)         # For network analysis
library(rgeos)          # To manipulate spatial objects
library(ggnetwork)      # To plot network using ggplot
library(sf)             # To plot spatial features
library(ggspatial)      # To add scale bars etc to plots

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load reference raster
r <- rast("03_Data/02_CleanData/ReferenceRaster.tif")

# Load source areas
area <- vect("03_Data/02_CleanData/SourceAreas.shp")

# Rasterize them to the reference raster
area_r <- terra::rasterize(area, y = r, field = "ID")

################################################################################
#### Inter-Patch Connectivity between Source Areas
################################################################################
# Load dispersal simulations
sims <- read_rds("03_Data/03_Results/DispersalSimulation.rds")

# Keep only desired columns
sims <- sims[, c("x", "y", "TrackID", "StepNumber", "SourceArea", "FloodLevel")]

# Count number of simulations per source area (should be 1000 each)
n_sample <- sims %>%
  select(TrackID, SourceArea, FloodLevel) %>%
  distinct() %>%
  count(SourceArea, FloodLevel) %>%
  pull(n) %>%
  unique()
print(n_sample)

# Make coordinates of simulated trajectories spatial
coordinates(sims) <- c("x", "y")
crs(sims) <- CRS("+init=epsg:4326")

# Identify through which designated areas the dispersers moved
visits <- data.frame(
    TrackID    = sims$TrackID
  , StepNumber = sims$StepNumber
  , FloodLevel = sims$FloodLevel
  , SourceArea = sims$SourceArea
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , Area       = raster::extract(raster(area_r), sims)
)

# Ignore any step in an "na" or "nan" area
visits <- subset(visits, !is.na(Area) & !is.nan(Area))

# Calculate, for each step, the distance to the first coordinate. We'll use this
# to determine how far (Euclidean Distance!!!) an individual had to disperse
# before reaching another area
visits <- visits %>%
  nest(data = -TrackID) %>%
  mutate(data = pbmclapply(data
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(x) {

      # Project coordinates
      coords <- reprojCoords(
          xy   = x[, c("x", "y")]
        , from = CRS("+init=epsg:4326")
        , to   = CRS("+init=epsg:32734")
      )

      # Compute distance to first coordinate
      first <- coords[1, ]
      distance <- sqrt((coords[, 1] - first[1]) ** 2 + (coords[, 2] - first[2]) ** 2)
      x$DistanceFromFirst <- distance

      # Return the resulting object
      return(x)
  })) %>%
  unnest(data)

# Identify how long it takes to reach the different areas
visits <- visits %>%
  rename(SourceArea = SourceArea, CurrentArea = Area) %>%
  group_by(TrackID, FloodLevel, SourceArea, CurrentArea) %>%
  summarize(
      StepNumber        = min(StepNumber)
    , DistanceFromFirst = min(DistanceFromFirst)
    , .groups           = "drop"
  ) %>%
  arrange(TrackID, StepNumber)

# Compute summary statistics by source area
summarizeVisits <- function(visits) {
  visits %>%
    group_by(FloodLevel, SourceArea, CurrentArea) %>%
    summarize(
        MeanStepNumber = mean(StepNumber)
      , SDStepNumber   = sd(StepNumber)
      , Frequency      = n()
      , .groups        = "drop"
    )
}

# Try it
summarizeVisits(visits)

# Generate bootstrap samples to compute standard deviation for the number of
# successful dispersals between source areas and the number of steps required
visits_nested <- nest(visits, Data = -c(TrackID, FloodLevel, SourceArea))
bootstrapped <- pbmclapply(1:1000, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  boot <- visits_nested %>%
    group_by(FloodLevel, SourceArea) %>%
    slice_sample(replace = T, n = n_sample) %>%
    unnest(Data) %>%
    summarizeVisits()
  boot$Bootstrap <- x
  return(boot)
})
bootstrapped <- do.call(rbind, bootstrapped)

# Bind and compute confidence intervals
visits_bootstrapped <- bootstrapped %>%
  group_by(FloodLevel, SourceArea, CurrentArea) %>%
  summarize(
      StepNumber   = mean(MeanStepNumber)
    , StepNumberSE = sd(MeanStepNumber)   # Bootstrapped sd is the metric's se
    , Freq         = mean(Frequency)
    , FreqSE       = sd(Frequency)        # Bootstrapped sd is the metric's se
    , .groups      = "drop"
  )

# Store visits to file
write_rds(visits, "03_Data/03_Results/InterpatchConnectivity.rds")
write_rds(visits_bootstrapped, "03_Data/03_Results/BootstrappedInterpatchConnectivity.rds")

# Reload them
visits <- read_rds("03_Data/03_Results/InterpatchConnectivity.rds")
visits_bootstrapped <- read_rds("03_Data/03_Results/BootstrappedInterpatchConnectivity.rds")

################################################################################
#### Additional Metrics: Number of Individuals Reaching any other Area
################################################################################
# Let's prepare a list of additional metrics we may want to report
metrics <- list()

# Get an idea of the number of individuals dispersing into any other area
visits %>%
  subset(SourceArea != CurrentArea & CurrentArea <= 6) %>%
  select(TrackID, FloodLevel) %>%
  distinct(TrackID, FloodLevel) %>%
  count(FloodLevel) %>%
  mutate(Percent = n / 6000)

# Bootstrap number of dispersers reaching another source area
visits_nested <- nest(visits, Data = -c(TrackID, FloodLevel))
n_sample <- unique(count(visits_nested, FloodLevel)$n)
bootstrapped <- pbmclapply(1:1000, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  boot <- visits_nested %>%
    group_by(FloodLevel) %>%
    slice_sample(replace = T, n = n_sample) %>%
    ungroup() %>%
    mutate(TrackID = 1:n()) %>%
    unnest(Data) %>%
    subset(SourceArea != CurrentArea & CurrentArea <= 6) %>%
    select(TrackID, FloodLevel) %>%
    distinct() %>%
    count(FloodLevel)
  boot$Bootstrap <- x
  return(boot)
})
metrics[[1]] <- bootstrapped %>%
  do.call(rbind, .) %>%
  group_by(FloodLevel) %>%
  summarize(
      Number   = mean(n)
    , NumberSE = sd(n)
  ) %>%
  mutate(combined = paste0(format(round(Number), big.mark = "'"), "$ \\pm $", sprintf("%.2f", round(NumberSE, 2)))) %>%
  mutate(NumberRev = rev(Number)) %>%
  mutate(Percent = round(Number / NumberRev * 100 - 100, 0))
names(metrics)[1] <- "Number individuals reaching other source areas at maximum and minimum extent"# Bootstrap number of dispersers reaching another source area

################################################################################
#### Additional Metrics: Dispersal Duration before Reaching any other Area
################################################################################
# Get an idea of the dispersal duration before reaching any other area
visits %>%
  subset(SourceArea != CurrentArea & CurrentArea <= 6) %>%
  select(TrackID, StepNumber, FloodLevel) %>%
  group_by(FloodLevel, TrackID) %>%
  summarize(StepNumber = min(StepNumber), .groups = "drop") %>%
  group_by(FloodLevel) %>%
  summarize(StepNumber = mean(StepNumber))

# Bootstrap the dispersal duration before reaching any other area
visits_nested <- nest(visits, Data = -c(TrackID, FloodLevel))
n_sample <- unique(count(visits_nested, FloodLevel)$n)
bootstrapped <- pbmclapply(1:1000, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  boot <- visits_nested %>%
    group_by(FloodLevel) %>%
    slice_sample(replace = T, n = n_sample) %>%
    ungroup() %>%
    mutate(TrackID = 1:n()) %>%
    unnest() %>%
    subset(SourceArea != CurrentArea & CurrentArea <= 6) %>%
    select(TrackID, StepNumber, FloodLevel) %>%
    group_by(FloodLevel, TrackID) %>%
    summarize(StepNumber = min(StepNumber), .groups = "drop") %>%
    group_by(FloodLevel) %>%
    summarize(StepNumber = mean(StepNumber))
  boot$Bootstrap <- x
  return(boot)
})
metrics[[2]] <- bootstrapped %>%
  do.call(rbind, .) %>%
  group_by(FloodLevel) %>%
  summarize(
      Number   = mean(StepNumber)
    , NumberSE = sd(StepNumber)
  ) %>%
  mutate(combined = paste0(format(round(Number), big.mark = "'"), "$ \\pm $", sprintf("%.2f", round(NumberSE, 2)))) %>%
  mutate(NumberRev = rev(Number)) %>%
  mutate(Percent = round(Number / NumberRev * 100 - 100, 0))
names(metrics)[2] <- "Steps until reaching any other area at maximum and minimum extent"

################################################################################
#### Additional Metrics: Number of Individuals Moving into area 6
################################################################################
# Get an idea of the number of individuals moving into source area 6
visits %>%
  subset(SourceArea != 6 & CurrentArea == 6) %>%
  select(TrackID, FloodLevel) %>%
  distinct() %>%
  count(FloodLevel)

# Bootstrap the number of individuals reaching area 6
visits_nested <- nest(visits, Data = -c(TrackID, FloodLevel))
n_sample <- unique(count(visits_nested, FloodLevel)$n)
bootstrapped <- pbmclapply(1:1000, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  boot <- visits_nested %>%
    group_by(FloodLevel) %>%
    slice_sample(replace = T, n = n_sample) %>%
    ungroup() %>%
    mutate(TrackID = 1:n()) %>%
    unnest(Data) %>%
    subset(SourceArea != 6 & CurrentArea == 6) %>%
    select(TrackID, FloodLevel) %>%
    distinct() %>%
    count(FloodLevel)
  boot$Bootstrap <- x
  return(boot)
})
metrics[[3]] <- bootstrapped %>%
  do.call(rbind, .) %>%
  group_by(FloodLevel) %>%
  summarize(
      Number   = mean(n)
    , NumberSE = sd(n)
  ) %>%
  mutate(combined = paste0(format(round(Number), big.mark = "'"), "$ \\pm $", sprintf("%.2f", round(NumberSE, 2)))) %>%
  mutate(NumberRev = rev(Number)) %>%
  mutate(Percent = round(Number / NumberRev * 100 - 100, 0))
names(metrics)[3] <- "Number individuals reaching source area 6 at maximum and minimum extent"

################################################################################
#### Additional Metrics: Number of Steps before Moving into area 6
################################################################################
# Get an idea of the dispersal duration until individuals reach area 6
visits %>%
  subset(SourceArea != 6 & CurrentArea == 6) %>%
  group_by(TrackID, FloodLevel) %>%
  summarize(StepNumber = min(StepNumber), .groups = "drop") %>%
  group_by(FloodLevel) %>%
  summarize(StepNumber = mean(StepNumber))

# Bootstrap the dispersal duration unti individuals reach area 6
visits_nested <- nest(visits, Data = -c(TrackID, FloodLevel))
n_sample <- unique(count(visits_nested, FloodLevel)$n)
bootstrapped <- pbmclapply(1:1000, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  boot <- visits_nested %>%
    group_by(FloodLevel) %>%
    slice_sample(replace = T, n = n_sample) %>%
    ungroup() %>%
    mutate(TrackID = 1:n()) %>%
    unnest(Data) %>%
    subset(SourceArea != 6 & CurrentArea == 6) %>%
    select(TrackID, FloodLevel, StepNumber) %>%
    group_by(TrackID, FloodLevel) %>%
    summarize(StepNumber = min(StepNumber), .groups = "drop") %>%
    group_by(FloodLevel) %>%
    summarize(StepNumber = mean(StepNumber))
  boot$Bootstrap <- x
  return(boot)
})
metrics[[4]] <- bootstrapped %>%
  do.call(rbind, .) %>%
  group_by(FloodLevel) %>%
  summarize(
      Number   = mean(StepNumber)
    , NumberSE = sd(StepNumber)
  ) %>%
  mutate(combined = paste0(format(round(Number), big.mark = "'"), "$ \\pm $", sprintf("%.2f", round(NumberSE, 2)))) %>%
  mutate(NumberRev = rev(Number)) %>%
  mutate(Percent = round(Number / NumberRev * 100 - 100, 0))
names(metrics)[4] <- "Steps until reaching source area 6 at maximum and minimum extent"

################################################################################
#### Additional Metrics: Metrics for improved Connectivity
################################################################################
# For a few selected source areas we get an increase in connectivity / dispersal
# success
better <- visits_bootstrapped %>%
  subset(CurrentArea <= 6 & SourceArea != CurrentArea) %>%
  select(SourceArea, CurrentArea, FloodLevel, Freq) %>%
  pivot_wider(names_from = FloodLevel, values_from = Freq) %>%
  mutate(Percent = Max / Min * 100 - 100) %>%
  mutate(Better = Percent > 0) %>%
  subset(Better) %>%
  arrange(desc(Percent))

# Let's create a combined dataframe (including the SEs)
metrics[[5]] <- visits_bootstrapped %>%
  select(FloodLevel, SourceArea, CurrentArea, Freq, FreqSE) %>%
  left_join(., better, by = c("SourceArea", "CurrentArea")) %>%
  subset(!is.na(Better)) %>%
  subset(Percent == max(Percent)) %>%
  mutate(Percent = round(Percent)) %>%
  mutate(combined = paste0(format(round(Freq), big.mark = "'"), "$ \\pm $", sprintf("%.2f", round(FreqSE, 2))))
names(metrics)[5] <- "Improved connectivity from source area 1 to 5"

################################################################################
#### Additional Metrics: Temporary Emigration
################################################################################
# How many individuals emigrate, at least temporarily? Also record through
# which area they emigrate
visits %>%
  subset(CurrentArea > 6) %>%
  group_by(TrackID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  count(FloodLevel, CurrentArea)

# Bootstrap that number
visits_nested <- nest(visits, Data = -c(TrackID, FloodLevel))
n_sample <- unique(count(visits_nested, FloodLevel)$n)
bootstrapped <- pbmclapply(1:1000, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  boot <- visits_nested %>%
    group_by(FloodLevel) %>%
    slice_sample(replace = T, n = n_sample) %>%
    ungroup() %>%
    mutate(TrackID = 1:n()) %>%
    unnest(Data) %>%
    subset(CurrentArea > 6) %>%
    group_by(TrackID) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    count(FloodLevel, CurrentArea)
  boot$Bootstrap <- x
  return(boot)
})
metrics[[6]] <- bootstrapped %>%
  do.call(rbind, .) %>%
  group_by(FloodLevel, Bootstrap) %>%
  summarize(
    n = sum(n)
  ) %>%
  ungroup() %>%
  group_by(FloodLevel) %>%
  summarize(
      Number   = mean(n)
    , NumberSE = sd(n)
    , .groups  = "drop"
  ) %>%
  mutate(combined = paste0(format(round(Number), big.mark = "'"), "$ \\pm $", sprintf("%.2f", round(NumberSE, 2)))) %>%
  mutate(NumberRev = rev(Number)) %>%
  mutate(Percent = round(Number / NumberRev * 100 - 100, 0))
names(metrics)[6] <- "Temporary Emigration"

# Let's also store the metrics per source area separately
emigration_temp <- bootstrapped %>%
  do.call(rbind, .) %>%
  group_by(FloodLevel, CurrentArea) %>%
  summarize(
      Number   = mean(n)
    , NumberSE = sd(n)
    , .groups  = "drop"
  ) %>%
  mutate(Type = "Temporary")

################################################################################
#### Additional Metrics: Permanent Emigration
################################################################################
# How many individuals emigrate, at least temporarily? Also record through
# which area they emigrate
# How many dispersers leave the main study area for good?
visits %>%
  group_by(TrackID) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  subset(CurrentArea > 6) %>%
  count(FloodLevel)

# Bootstrap that number
visits_nested <- nest(visits, Data = -c(TrackID, FloodLevel))
n_sample <- unique(count(visits_nested, FloodLevel)$n)
bootstrapped <- pbmclapply(1:1000, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  boot <- visits_nested %>%
    group_by(FloodLevel) %>%
    slice_sample(replace = T, n = n_sample) %>%
    ungroup() %>%
    mutate(TrackID = 1:n()) %>%
    unnest(Data) %>%
    group_by(TrackID) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    subset(CurrentArea > 6) %>%
    count(FloodLevel, CurrentArea)
  boot$Bootstrap <- x
  return(boot)
})
metrics[[7]] <- bootstrapped %>%
  do.call(rbind, .) %>%
  group_by(FloodLevel, Bootstrap) %>%
  summarize(
    n = sum(n)
  ) %>%
  ungroup() %>%
  group_by(FloodLevel) %>%
  summarize(
      Number   = mean(n)
    , NumberSE = sd(n)
    , .groups  = "drop"
  ) %>%
  mutate(combined = paste0(format(round(Number), big.mark = "'"), "$ \\pm $", sprintf("%.2f", round(NumberSE, 2)))) %>%
  mutate(NumberRev = rev(Number)) %>%
  mutate(Percent = round(Number / NumberRev * 100 - 100, 0))
names(metrics)[7] <- "Permanent Emigration"

# Let's also store the metrics per source area separately
emigration_perm <- bootstrapped %>%
  do.call(rbind, .) %>%
  group_by(FloodLevel, CurrentArea) %>%
  summarize(
      Number   = mean(n)
    , NumberSE = sd(n)
    , .groups  = "drop"
  ) %>%
  mutate(Type = "Permanent")

# Put source-areas specific metrics together and store them. We can use them
# for a nice plot
emigration <- rbind(emigration_temp, emigration_perm)

# Store this to file
write_rds(emigration, "03_Data/03_Results/Emigration.rds")

################################################################################
#### Store Metrics
################################################################################
# Store all the metrics to file
writeLines(metrics[[1]]$combined[1], "04_Manuscript/99_NumberReachingOthersMaxFlood.tex")
writeLines(metrics[[1]]$combined[2], "04_Manuscript/99_NumberReachingOthersMinFlood.tex")
writeLines(metrics[[2]]$combined[1], "04_Manuscript/99_StepsToReachingOthersMaxFlood.tex")
writeLines(metrics[[2]]$combined[2], "04_Manuscript/99_StepsToReachingOthersMinFlood.tex")
writeLines(metrics[[3]]$combined[1], "04_Manuscript/99_NumberReachingArea6MaxFlood.tex")
writeLines(metrics[[3]]$combined[2], "04_Manuscript/99_NumberReachingArea6MinFlood.tex")
writeLines(metrics[[4]]$combined[1], "04_Manuscript/99_StepsToReachingArea6MaxFlood.tex")
writeLines(metrics[[4]]$combined[2], "04_Manuscript/99_StepsToReachingArea6MinFlood.tex")
writeLines(metrics[[5]]$combined[1], "04_Manuscript/99_ImprovedConnectivity1to5MaxFlood.tex")
writeLines(metrics[[5]]$combined[2], "04_Manuscript/99_ImprovedConnectivity1to5MinFlood.tex")
writeLines(metrics[[6]]$combined[1], "04_Manuscript/99_TemporaryEmigrationMaxFlood.tex")
writeLines(metrics[[6]]$combined[2], "04_Manuscript/99_TemporaryEmigrationMinFlood.tex")
writeLines(metrics[[7]]$combined[1], "04_Manuscript/99_PermanentEmigrationMaxFlood.tex")
writeLines(metrics[[7]]$combined[2], "04_Manuscript/99_PermanentEmigrationMinFlood.tex")

# Let's also store the percentage differences between min and max flood
writeLines(paste0(abs(metrics[[1]]$Percent[[1]])), "04_Manuscript/99_NumberReachingOthersPercentage.tex")
writeLines(paste0(abs(metrics[[2]]$Percent[[1]])), "04_Manuscript/99_StepsToReachingOthersPercentage.tex")
writeLines(paste0(abs(metrics[[3]]$Percent[[1]])), "04_Manuscript/99_NumberReachingArea6Percentage.tex")
writeLines(paste0(abs(metrics[[4]]$Percent[[1]])), "04_Manuscript/99_StepsToReachingArea6Percentage.tex")
writeLines(paste0(abs(metrics[[5]]$Percent[1])), "04_Manuscript/99_ImprovedConnectivity1to5Percentage.tex")
writeLines(paste0(abs(metrics[[6]]$Percent[1])), "04_Manuscript/99_TemporaryEmigrationPercentage.tex")
writeLines(paste0(abs(metrics[[7]]$Percent[1])), "04_Manuscript/99_PermanentEmigrationPercentage.tex")

# ################################################################################
# #### Number Connections in Relation to Distance
# ################################################################################
# # Function to determine how many of the simulated individuals reached another
# # national park, as well as the average dispersal duration required for this
# getConnections <- function(min_distance, grouping = "SourceArea") {
#   reached <- visits %>%
#       subset(DistanceFromFirst >= min_distance & CurrentArea != SourceArea) %>%
#       group_by(FloodLevel, TrackID, SourceArea, CurrentArea) %>%
#       summarize(
#           StepNumber = min(StepNumber)
#         , .groups    = "drop"
#       ) %>%
#       group_by_("FloodLevel", grouping) %>%
#       summarize(
#           MeanStepNumber = mean(StepNumber)
#         , SDStepNumber   = sd(StepNumber)
#         , Frequency      = length(unique(TrackID))
#         , .groups        = "drop"
#       )
#   return(reached)
# }
# 
# # Try it
# getConnections(0, grouping = "SourceArea")
# getConnections(0, grouping = "CurrentArea")
# 
# # Run the function for different distances to see how the success rate
# # deacreases if one only considers more distant areas
# results <- tibble(MinDistance = seq(0, 300 * 1000, length.out = 25))
# results$Results <- pbmclapply(results$MinDistance
#     , ignore.interactive = T
#     , mc.cores           = detectCores() - 1
#     , FUN                = function(x) {
#       conns1 <- getConnections(x, grouping = "SourceArea")
#       conns2 <- getConnections(x, grouping = "CurrentArea")
#       conns1$Type <- "IntoOther"
#       conns2$Type <- "FromOther"
#       names(conns1)[2] <- "Area"
#       names(conns2)[2] <- "Area"
#       conns <- rbind(conns1, conns2)
#       conns$Type <- factor(conns$Type, levels = c("IntoOther", "FromOther"))
#     return(conns)
# })
# 
# # Unnest and convert fractions to percentages
# reached <- unnest(results, Results)
# reached$MinDistance <- reached$MinDistance / 1000
# reached <- subset(reached, Area <= 6)
# 
# # Plot
# p <- ggplot(reached, aes(x = MinDistance, y = Frequency, col = FloodLevel)) +
#   geom_hline(yintercept = 1000, lty = 2, col = "gray") +
#   geom_line(size = 0.3) +
#   geom_point(size = 1) +
#   theme_classic() +
#   xlab("Minimum Distance Considered (km)") +
#   ylab("% Trajectories Reaching\nanother National Park") +
#   theme(
#       panel.grid.major = element_line(colour = "gray90", size = 0.1)
#     , panel.grid.minor = element_line(colour = "gray90", size = 0.1)
#     , legend.position = c(0.85, 0.7)
#   ) +
#   facet_grid(Area ~ Type) +
#   scale_color_manual(values = c("cornflowerblue", "orange"))
# p
