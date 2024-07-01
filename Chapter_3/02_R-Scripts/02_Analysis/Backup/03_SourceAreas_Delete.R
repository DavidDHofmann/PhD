################################################################################
#### Source Areas
################################################################################
# Description: Generation of source area polygons

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_3")

# Load required packages
library(terra)      # To handle spatial data
library(tidyterra)  # To get tidy terra functions
library(tidyverse)  # To wrangle data

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load some landmarks according to which we will delineate the source areas
wat <- vect("03_Data/02_CleanData/MajorWaters.gpkg")
roa <- vect("03_Data/02_CleanData/Roads.gpkg")
riv <- vect("03_Data/02_CleanData/MajorRivers.gpkg")
pro <- vect("03_Data/02_CleanData/Protected.gpkg") %>%
  filter(Name %in% c("Moremi", "Hwange"))

# Specify centroids for the different source areas
src <- data.frame(
    x  = c(22.7, 23.75, 26.25)
  , y  = c(-20.1, -19.4, -19.1)
  , ID = 1:3
)
src <- vect(as.matrix(src[, c("x", "y")]), crs = crs(wat), atts = src)
src <- buffer(src, width = 20000)

# Visualize all
plot(wat, col = "cornflowerblue", border = NA)
plot(riv, col = "cornflowerblue", border = NA, add = T)
plot(roa, col = "gray50", add = T)
plot(pro, add = T, col = adjustcolor("#039d09", alpha.f = 0.5), lty = 2)
plot(src, add = T)
text(src, label = "ID")

# Store source areas
writeVector(src, "03_Data/02_CleanData/Sources.gpkg", overwrite = T)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/03_SourceAreas.rds")

# Print to terminal
cat("Done :) \n")
