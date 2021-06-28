################################################################################
#### Define Study Areas for the Different Species
################################################################################
# Description: Use the Movement data downloaded from Movebank to define study
# areas for each of the considered species
# Author: David Hofmann
# Date: June 2021

# Clear R's brain
rm(list = ls())

# # Load required packages
library(raster)       # To handle spatial data
library(rgeos)        # To manipulate spatial data
library(tidyverse)    # For data wrangling
library(rworldmap)    # To plot a map of the world

# Set the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_5")

# Load movement data
dat <- dir(pattern = "csv$", path = "03_Data/01_RawData/MOVEBANK", full.names = T)
dat <- lapply(dat, function(x){
  data <- read.csv(x, stringsAsFactors = F)
  data <- dplyr::select(data, c(
      Timestamp = timestamp
    , x         = location.long
    , y         = location.lat
    , Species   = individual.taxon.canonical.name
    , ID        = individual.local.identifier
    , Source    = study.name
  ))
  return(data)
})
dat <- do.call(rbind, dat)

# Let's make sure that each individual gets a unique ID
dat$ID <- dat %>% group_indices(Species, ID)

# Make proper timestamps
dat$Timestamp <- as.POSIXct(as.character(dat$Timestamp), tz = "UTC")

# Remove NA coordinates
dat <- subset(dat, !is.na(x) & !is.na(y))

# Determine the extent for each species
dat <- dat %>%
  nest(data = -Species) %>%
  mutate(extent = map(data, function(x){
    xmin <- min(x$x)
    xmax <- max(x$x)
    ymin <- min(x$y)
    ymax <- max(x$y)
    xlen <- xmax - xmin
    ylen <- ymax - ymin
    fac <- 0.1
    ext <- extent(c(xmin - fac * xlen, xmax + fac * xlen, ymin - fac * ylen, ymax + fac * ylen))
    ext <- as(ext, "SpatialPolygons")
    return(ext)
  }))

# Put extents into a single object
ext <- rbind(dat$extent[[1]], dat$extent[[2]], makeUniqueIDs = T)
ext <- as(ext, "SpatialPolygonsDataFrame")
crs(ext) <- CRS("+init=epsg:4326")

# Plot the extents on a world map
world <- getMap(resolution = "coarse")
plot(ext, border = "red")
plot(world, add = T)

# Try to download some elevation data
library(elevatr)
test <- get_elev_raster(ext[2, ], z = 9)
test <- crop(test, ext[2, ])
plot(test)
points(dat$data[[2]][, c("x", "y")], pch = 16, cex = 0.3)
test <- getData("GADM", country = "CHE", level = 0)
plot(ext)
