################################################################################
#### Step Selection Function - Extraction of Covariates
################################################################################
# Description: In this script we will extract all covariates underlying the
# generated steps

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load packages
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(lubridate)    # To handle dates
library(amt)          # To coerce gps fixes to steps
library(raster)       # To handle spatial data
library(rgdal)        # To handle spatial data
library(pbmcapply)    # For multicore abilities
library(spatstat)     # For point patter analysis (calculate distance)
library(maptools)     # For point patter analysis (calculate distance)
library(rgeos)        # For spatial manipulation

# Load the generated steps
dat1 <- readOGR("03_Data/02_CleanData/00_General_Dispersers_POPECOL(iSSF).shp")
dat2 <- readOGR("03_Data/02_CleanData/00_General_Dispersers_POPECOL(TiSSF).shp")
lines <- rbind(dat1, dat2)

################################################################################
#### Land Cover - Globeland
################################################################################
# Load the copernicus land cover map
dat <- raster("03_Data/02_CleanData/01_LandCover_LandCover_GLOBELAND.tif")
inf <- read_csv("03_Data/02_CleanData/01_LandCover_LandCover_GLOBELAND.csv")

# Split categories into different layers
dat <- layerize(dat)

# Put nice layernames
names(dat) <- inf$Class

# Extract land cover along each line
extracted <- extrCov(dat, lines)

# Add nice column names
names(extracted) <- paste0("Globeland_", names(dat))

# Let's look at the result
head(extracted)

# Add the data to the lines
lines@data <- cbind(lines@data, extracted)

################################################################################
#### Land Cover - Copernicus
################################################################################
# Load the copernicus land cover map
dat <- raster("03_Data/02_CleanData/01_LandCover_LandCover_COPERNICUS.tif")
inf <- read_csv("03_Data/02_CleanData/01_LandCover_LandCover_COPERNICUS.csv")

# Split categories into different layers
dat <- layerize(dat)

# Put nice layernames
names(dat) <- inf$Class

# Extract land cover along each line
extracted <- extrCov(dat, lines)

# Add nice column names
names(extracted) <- paste0("Copernicus_", names(dat))

# Let's look at the result
head(extracted)

# Add the data to the lines
lines@data <- cbind(lines@data, extracted)

################################################################################
#### LandCover - DistanceToWater
################################################################################
# Transform the lines to utm
lines <- spTransform(lines, CRS("+init=epsg:32734"))

# To extract distances we need to coerce the lines to a psp object
linesppp <- lapply(lines@lines, function(z){lapply(z@Lines, as.psp)})
linesppp <- do.call("c", linesppp)

# We also need to create points on the lines to extract average distances. We
# can do so by setting a regular distance (100 meters in this case)
linesppp <- lapply(linesppp, pointsOnLines, eps = 100)

# Load the merged water cover dataset
dat <- stack("03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd")

# Extract dates from layernames
dates <- names(dat) %>%
  substr(start = 2, stop = 11) %>%
  as.Date(format = "%Y.%m.%d")

# Prepare a list that stores a ppp layer for water (Code 1) for each floodmap
datppp <- suppressMessages(
  pbmclapply(
      X                   = 1:nlayers(dat)
    , mc.cores            = detectCores() - 1
    , ignore.interactive  = T
    , FUN                 = function(x){
      points <- rasterToPoints(dat[[x]], fun = function(z){z == 1}, spatial = T)
      points <- spTransform(points, CRS("+init=epsg:32734"))
      points <- as(points, "ppp")
      return(points)
  })
)

# Calculate the average distance to water on the ppp object that is closest in
# date to the actual step
lines$DistanceToWater <- suppressMessages(
  pbmclapply(
      X                  = 1:nrow(lines)
    , mc.cores           = detectCores() / 2
    , ignore.interactive = T
    , FUN                = function(x){
    index <- which.min(abs(as.Date(lines$t1_[x]) - dates))[1]
    distance <- nncross(linesppp[[x]], datppp[[index]])
    distance <- mean(distance$dist)
    return(distance)
    gc()
  }) %>% do.call(rbind, .)
)

# Transform lines back to WGS84
lines <- spTransform(lines, CRS("+init=epsg:4326"))

################################################################################
#### LandCover - Water
################################################################################
# Load the merged water cover dataset
dat <- stack("03_Data/02_CleanData/01_LandCover_WaterCover_MERGED.grd")

# Extract dates from layernames
dates <- names(dat) %>%
  substr(start = 2, stop = 11) %>%
  as.Date(format = "%Y.%m.%d")

# Now we can extract the percentage cover of Water along each step
extracted <- extrCov(dat, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(dates)

# Let's look at the result
head(extracted)

# We only want to keep the values from the dates that are closest in time to the
# steps
lines$Water <- pbmclapply(1:nrow(lines)
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(x){
    index <- which.min(abs(as.Date(lines$t1_[x]) - dates))[1]
    value <- extracted[x, index]
    return(value)
}) %>% do.call(c, .)

################################################################################
#### LandCover - Trees
################################################################################
# Load the treecover map
dat <- stack("03_Data/02_CleanData/01_LandCover_TreeCover_MODIS.grd")

# Extract dates from layernames
dates <- names(dat) %>%
  substr(start = 2, stop = 11) %>%
  as.Date(format = "%Y.%m.%d")

# Now we can extract the percentage cover of water along each step
extracted <- extrCov(dat, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(dates)

# Keep only values closest in date
lines$Trees <- pbmclapply(1:nrow(lines)
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(x){
    index <- which.min(abs(as.Date(lines$t1_[x]) - dates))[1]
    value <- extracted[x, index]
    return(value)
}) %>% do.call(c, .)

################################################################################
#### LandCover - Shrubs
################################################################################
# Load the shrubcover map
dat <- stack("03_Data/02_CleanData/01_LandCover_NonTreeVegetation_MODIS.grd")

# Extract dates from layernames
dates <- names(dat) %>%
  substr(start = 2, stop = 11) %>%
  as.Date(format = "%Y.%m.%d")

# Extract the average shrub cover along each line
extracted <- extrCov(dat, lines)

# For completeness we might want to add the dates into the dataframe
names(extracted) <- as.character(dates)

# Let's look at the result
head(extracted)

# Keep only values closest in date
lines$Shrubs <- pbmclapply(1:nrow(lines)
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(x){
    index <- which.min(abs(as.Date(lines$t1_[x]) - dates))[1]
    value <- extracted[x, index]
    return(value)
}) %>% do.call(c, .)

################################################################################
#### Land Use - Protection
################################################################################
# Load protection zones
dat <- raster("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS(3Classes).tif")
inf <- read_csv("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS(3Classes).csv")

# Separate layers
dat <- layerize(dat)

# Extract the percentage coverage
extracted <- extrCov(dat, lines)

# Add nice column names
names(extracted) <- inf$Class

# Put the extracted values into the dataframe
lines@data <- cbind(lines@data, extracted)

################################################################################
#### Anthropogenic - Road Crossings
################################################################################
# Load data
dat <- readOGR("03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.shp")

# Simplify to single geometry
dat <- gLineMerge(dat)

# Identify intersections with roads
lines$RoadCrossing <- as.vector(gIntersects(dat, lines, byid = T))

################################################################################
#### Anthropogenic - Various
################################################################################
# Load different anthropogenic layers
dat <- stack(c(
    "03_Data/02_CleanData/04_AnthropogenicFeatures_Villages_FACEBOOK.tif"
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_Villages_WORLDPOP.tif"
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToVillages_FACEBOOK.tif"
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToVillages_WORLDPOP.tif"
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToRoads_GEOFABRIK.tif"
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToHumans_FACEBOOK.tif"
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToHumans_WORLDPOP.tif"
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanDensity_FACEBOOK.tif"
  , "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanDensity_WORLDPOP.tif"
))

# Assign nice layernames
names(dat) <- c(
    "Facebook_Villages"
  , "Worldpop_Villages"
  , "Facebook_DistanceToVillages"
  , "Worldpop_DistanceToVillages"
  , "DistanceToRoads"
  , "Facebook_DistanceToHumans"
  , "Worldpop_DistanceToHumans"
  , "Facebook_HumanDensity"
  , "Worldpop_HumanDensity"
)

# Extract values
extracted <- extrCov(dat, lines)

# Assign nice names
names(extracted) <- names(dat)

# Put the extracted values into the dataframe
lines@data <- cbind(lines@data, extracted)

################################################################################
#### Anthropogenic - Human Influence (Facebook)
################################################################################
# Load human influence data
dat <- stack("03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluenceBuff_FACEBOOK.grd")

# Extract values
extracted <- extrCov(dat, lines)

# Assign names
names(extracted) <- paste0("Facebook_HumanInfluence", names(dat))

# Put the extracted values into the dataframe
lines@data <- cbind(lines@data, extracted)

################################################################################
#### Anthropogenic - Human Influence (Worldpop)
################################################################################
# Load human influence data
dat <- stack("03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluenceBuff_WORLDPOP.grd")

# Extract values
extracted <- extrCov(dat, lines)

# Assign names
names(extracted) <- paste0("Worldpop_HumanInfluence", names(dat))

# Put the extracted values into the dataframe
lines@data <- cbind(lines@data, extracted)

################################################################################
#### Storing
################################################################################
names(lines)

# Reorder the columns
lines@data <- dplyr::select(lines@data, c(
  , dog   = DogName
  , burst = id
  , step_id_ = step_d_
  , State
  , case_
  , x1_
  , x2_
  , y1_
  , y2_
  , t1_
  , t2_
  , dt_
  , sl_
  , ta_
  , absta_
  , everything()
  )
)

# Rename dog
names(lines)[1] <- "id"

# To store the files we need to coerce the duration column to a numeric
lines$DistanceToWater <- as.vector(lines$DistanceToWater)

# Prepare filenames
filename1 <- "00_General_Dispersers_POPECOL(iSSF_Extracted)"
filename2 <- "00_General_Dispersers_POPECOL(TiSSF_Extracted)"

# Split the data
lines1 <- subset(lines, method == "iSSF")
lines2 <- subset(lines, method == "TiSSF")

# Save the lines to a spatial lines dataframe
writeOGR(lines1
  , "03_Data/02_CleanData"
  , filename1
  , driver = "ESRI Shapefile"
  , overwrite = TRUE
)
writeOGR(lines2
  , "03_Data/02_CleanData"
  , filename2
  , driver = "ESRI Shapefile"
  , overwrite = TRUE
)

# Let's also store the data to a regular csv. We can use this file to restore
# the original column names since the ESRI shapefiles will store abbreviated
# names
write.csv(lines1@data, paste0("03_Data/02_CleanData/", filename1, ".csv"))
write.csv(lines2@data, paste0("03_Data/02_CleanData/", filename2, ".csv"))
