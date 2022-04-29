################################################################################
#### Covariate Extraction
################################################################################
# Description: Use the step selection data generated in the previous script and
# extract covariate values below each observed and random step. To speed up
# extraction speeds, we might extract covariate values along interpolated points
# instead of along lines. For dynamic covariates, we'll extract data from the
# layer that is closest in date to the actual step.

# Clear R's brain
rm(list = ls())

# Load packages
library(tidyverse)  # To wrangle data
library(terra)      # To handle spatial data
library(raster)     # To handle spatial data
library(lubridate)  # To handle dates
library(pbmcapply)  # To run stuff on multiple cores
library(velox)      # To extract data quickly

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

# Load step selection data
ssf <- read_csv("03_Data/02_CleanData/00_General_SSF.csv", show_col_types = F)
ssf <- subset(ssf, ID == "Abel")

# Let's generate a set of interpolated points along each of the steps
ssf$Points <- pbmclapply(
    X                  = 1:nrow(ssf)
  , ignore.interactive = T
  , FUN                = function(i) {
  p <- interpolatePoints(
      x1 = ssf$x[i]
    , x2 = ssf$x_to[i]
    , y1 = ssf$y[i]
    , y2 = ssf$y_to[i]
    , by = 250 / 111000
  )
  return(p)
})

################################################################################
#### Precipitation
################################################################################
# Load precipitation data
covariate <- stack("03_Data/02_CleanData/05_Climate_Precipitation.grd")

# Generate timestamps from the layernames
covariate_dates <- covariate %>%
  names() %>%
  substr(start = 2, stop = 21) %>%
  ymd_hms()

# Generate velox rasters
covariate <- velox(covariate)

# Extract precipitation data
ssf$Precipitation <- pbmclapply(
    X                  = 1:nrow(ssf)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(i) {
    index <- which.min(abs(ssf$Timestamp[i] - covariate_dates))[1]
    pts <- SpatialPoints(ssf$Points[[i]])
    extr <- covariate$extract_points(pts)
    extr <- colMeans(extr)
    extr <- extr[index]
    return(extr)
}) %>% do.call(c, .)

################################################################################
#### Temperature
################################################################################
# Load temperature data
covariate <- stack("03_Data/02_CleanData/05_Climate_Temperature.grd")

# Generate timestamps from the layernames
covariate_dates <- covariate %>%
  names() %>%
  substr(start = 2, stop = 21) %>%
  ymd_hms()

# Generate velox rasters
covariate <- velox(covariate)

# Extract precipitation data
ssf$Temperature <- pbmclapply(
    X                  = 1:nrow(ssf)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(i) {
    index <- which.min(abs(ssf$Timestamp[i] - covariate_dates))[1]
    pts <- SpatialPoints(ssf$Points[[i]])
    extr <- covariate$extract_points(pts)
    extr <- colMeans(extr)
    extr <- extr[index]
    return(extr)
}) %>% do.call(c, .)

################################################################################
#### Water (Flood, Pans, Distance)
################################################################################
# Load floodmaps
covariate <- stack("03_Data/02_CleanData/01_LandCover_WaterCoverDynamic.grd")

# Extract dates from layernames
covariate_dates <- covariate %>%
  names() %>%
  substr(start = 2, stop = 11) %>%
  ymd(tz = "UTC")

# Generate velox rasters
covariate <- velox(covariate)

# Extract precipitation data
ssf$Water <- pbmclapply(
    X                  = 1:nrow(ssf)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(i) {
    index <- which.min(abs(ssf$Timestamp[i] - covariate_dates))[1]
    pts <- SpatialPoints(ssf$Points[[i]])
    extr <- covariate$extract_points(pts)
    extr <- colMeans(extr)
    extr <- extr[index]
    return(extr)
}) %>% do.call(c, .)

################################################################################
#### Shrubs / Grassland
################################################################################
# Load shrubs / grassland data
covariate <- stack("03_Data/02_CleanData/01_LandCover_ShrubCoverDynamic.grd")

# Extract dates from layernames
covariate_dates <- covariate %>%
  names() %>%
  substr(start = 9, stop = 18) %>%
  ymd(tz = "UTC")

# Generate velox rasters
covariate <- velox(covariate)

# Extract precipitation data
ssf$Shrubs <- pbmclapply(
    X                  = 1:nrow(ssf)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(i) {
    index <- which.min(abs(ssf$Timestamp[i] - covariate_dates))[1]
    pts <- SpatialPoints(ssf$Points[[i]])
    extr <- covariate$extract_points(pts)
    extr <- colMeans(extr)
    extr <- extr[index]
    return(extr)
}) %>% do.call(c, .)

################################################################################
#### Trees
################################################################################
# Load trees data
covariate <- stack("03_Data/02_CleanData/01_LandCover_TreeCoverDynamic.grd")

# Extract dates from layernames
covariate_dates <- covariate %>%
  names() %>%
  substr(start = 8, stop = 17) %>%
  ymd(tz = "UTC")

# Generate velox rasters
covariate <- velox(covariate)

# Extract precipitation data
ssf$Trees <- pbmclapply(
    X                  = 1:nrow(ssf)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(i) {
    index <- which.min(abs(ssf$Timestamp[i] - covariate_dates))[1]
    pts <- SpatialPoints(ssf$Points[[i]])
    extr <- covariate$extract_points(pts)
    extr <- colMeans(extr)
    extr <- extr[index]
    return(extr)
}) %>% do.call(c, .)

################################################################################
#### Human Influence
################################################################################
# Load human influence data
covariate <- raster("03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence.tif")


################################################################################
#### CONTINUE HERE
################################################################################

# Extract dates from layernames
covariate_dates <- covariate %>%
  names() %>%
  substr(start = 8, stop = 17) %>%
  ymd(tz = "UTC")

# Generate velox rasters
covariate <- velox(covariate)

# Extract precipitation data
ssf$Trees <- pbmclapply(
    X                  = 1:nrow(ssf)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(i) {
    index <- which.min(abs(ssf$Timestamp[i] - covariate_dates))[1]
    pts <- SpatialPoints(ssf$Points[[i]])
    extr <- covariate$extract_points(pts)
    extr <- colMeans(extr)
    extr <- extr[index]
    return(extr)
}) %>% do.call(c, .)

################################################################################
#### NDVI
################################################################################

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/02_Analysis/01_CovariateExtraction.rds")
cat("Done :)\n")
