################################################################################
#### Lookup Table
################################################################################
# Description: Prepare a lookup table that we can use to generate the landscape
# for a desired date

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_3")

# Load required packages
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(terra)          # To handle spatial data

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load covariate layers
covs <- "03_Data/02_CleanData/00_SeasonalAggregates/Covariates.rds" %>%
  read_rds() %>%
  mutate(Raster = map(Filename, rast))

# We need to apply a little trick here. The covariates should "wrap around" at
# the end of the year. We therefore need at least one timestamp before the first
# date in 2000 and one after the last date in 2000. We'll simply reuse the
# layers that match the same date, but update the year.
covs <- covs %>%
  select(Covariate, Layerdates, Raster) %>%
  mutate(Layerdates = map(Layerdates, function(x) {

    # Assign layer numbers
    new <- data.frame(Layerdate = x, Layer = 1:length(x))

    # First entry is going to be the same timestamp as the last, but mapped to
    # an earlier year. Last entry is goign to be the same as the first, but
    # mapped to a later year.
    wrapped_first <- last(new)
    wrapped_last  <- first(new)
    wrapped_first$Layerdate <- update(wrapped_first$Layerdate, year = 1999)
    wrapped_last$Layerdate  <- update(wrapped_last$Layerdate, year = 2001)

    # Put all together
    new <- rbind(wrapped_first, new, wrapped_last)
    return(new)
  }))

# Check it
covs$Layerdates[[1]]

# Visualize the covariate dates
covs %>%
  select(Covariate, Layerdates) %>%
  subset(!(Covariate %in% c("Humans", "ForestCover"))) %>%
  unnest(Layerdates) %>%
  ggplot(aes(x = Layerdate, y = Covariate)) +
    geom_point() +
    theme_minimal()

# Span a sequence of timestamps the simulation. These are going to be the
# timestamps that we will simulate dispersal for. Note that the year has no real
# meaning as we'll just loop from year 2000 to year 2000 forever.
simul_dates <- tibble(
  Timestamp = seq(
        from = ymd_hm("2000-01-01 03:00")
      , to   = ymd_hm("2000-12-31 23:00")
      , by   = "4 hours"
    )
  ) %>% subset(hour(Timestamp) %in% c(3, 7, 15, 19, 23))

# Prepare a lookup table that gives the index of the layer that is closest in
# date
lookup <- expand_grid(
    Timestamp = simul_dates$Timestamp
  , Covariate = unique(covs$Covariate)
)

# Find the closest dates
lookup <- covs %>%
  left_join(lookup, ., by = "Covariate") %>%
  unnest(Layerdates) %>%
  mutate(TimeDifference = abs(difftime(Timestamp, Layerdate, units = "hours"))) %>%
  group_by(Timestamp, Covariate) %>%
  filter(TimeDifference == min(TimeDifference)) %>%
  ungroup()

# Get a feeling for the maximum timelags
lookup %>%
  group_by(Covariate) %>%
  summarize(
      min    = min(TimeDifference)
    , mean   = mean(TimeDifference)
    , max    = max(TimeDifference)
    , perc25 = quantile(TimeDifference, 0.25)
    , perc50 = quantile(TimeDifference, 0.50)
    , perc75 = quantile(TimeDifference, 0.75)
  )

# Simplify the table a bit
lookup <- lookup %>%
  select(Timestamp, Covariate, Layer) %>%
  pivot_wider(names_from = Covariate, values_from = Layer)

# Store the table to file
write_rds(lookup, "03_Data/02_CleanData/00_SeasonalAggregates/LookupTable.rds")

# # Run through the lookup table and prepare the associated landscapes
# backup <- lookup
# lookup <- lookup[1:100, ]
#
# cat("Generating landscapes for the different simulation timestamps...\n")
# pb <- txtProgressBar(min = 0, max = nrow(lookup), style = 3)
# lookup$SpatialRasterCollection <- lapply(1:nrow(lookup), function(i) {
#   land <- generateLandscape(lookup$Timestamp[i], covs, lookup)
#   land <- wrap(land)
#   setTxtProgressBar(pb, i)
#   return(land)
# }) %>% setNames(lookup$Timestamp)
#
# # Store the table to file
# write_rds(lookup, "03_Data/02_CleanData/00_SeasonalAggregates/LookupTable.rds")
#
# # Visualize one of the landscapes
# plot(lookup$SpatialRasterCollection[[1]], col = hcl.colors(100))
# plot(lookup$SpatialRasterCollection[[1000]], col = hcl.colors(100))
