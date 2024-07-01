################################################################################
#### Source Areas from which to Simulate Dispersal
################################################################################
# Description: Generate source areas from which dispersers will be simulated

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_3")

# Load required packages
library(tidyverse)      # To wrangle data
library(tidyterra)      # To wrangle spatial data
library(terra)          # To handle spatial data

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load protected areas and water shapefiles
prot <- vect("03_Data/02_CleanData/Protected.gpkg")
water <- vect("03_Data/02_CleanData/MajorWaters.gpkg")

# Generate polygons for the three areas from which we want to simulate dispersal
# from
source <- prot %>%
  filter(Name %in% c("Nxai Pan", "Moremi", "NG/26", "NG/29", "NG/33", "NG/34")) %>%
  filter(Desig != "Forest Reserve") %>%
  mutate(Group = case_when(
      Name %in% c("NG/26", "NG/29") ~ "West Delta"
    , Name %in% c("Moremi", "NG/33", "NG/34") ~ "East Delta"
    , .default = "Nxai Pan"
  )) %>%
  aggregate(by = "Group") %>%
  erase(water) %>%
  disagg() %>%
  mutate(Area = expanse(., unit = "km")) %>%
  filter(Area %in% sort(Area, decreasing = T)[1:3]) %>%
  tidyterra::select(Name = Group)

# Visualize them
# plot(source, col = adjustcolor("darkgreen", alpha.f = 0.5), border = NA)
# plot(water, add = T, col = "cornflowerblue", border = NA)
# text(source, label = "Name", col = "darkgreen")

# Store the source areas to file
writeVector(source, "03_Data/02_CleanData/Sources.gpkg")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/21_SourceAreas.rds")
cat("Done :)\n")
