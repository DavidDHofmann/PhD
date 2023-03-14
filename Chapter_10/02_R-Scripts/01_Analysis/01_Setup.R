################################################################################
#### Setup
################################################################################
# Description: Set up folder structure for this project and copy all the
# relevant data

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_10")

# Load required packages
library(terra)          # To handle spatial data
library(raster)         # To handle spatial data
library(sf)             # To handle spatial data
library(glmmTMB)        # To handle glmmTMB models
library(tidyverse)      # To wrangle data

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Generate an extent that represents our study area so that we can crop the data
# accordingly
# ext <- ext(21.5, 27.5, -21, -17)
ext <- ext(20.5, 26.5, -21.5, -17.5)

# Copy all files that are needed from the first PhD Chapter
old_chap <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
new_chap <- "/home/david/ownCloud/University/15. PhD/Chapter_10"

################################################################################
#### Setup Directories
################################################################################
# Let's prepare the main folders that we need
dir.create("03_Data/01_RawData", showWarnings = F)
dir.create("03_Data/02_CleanData", showWarnings = F)
dir.create("03_Data/03_Results", showWarnings = F)
dir.create("04_Manuscript", showWarnings = F)
dir.create("05_Presentation", showWarnings = F)

################################################################################
#### Copy Data from Chapter 1
################################################################################
# Print to terminal
cat("Copying files to project directory...\n")

# Copy the most parsimonious movement model from chapter one
mod <- read_rds(file.path(old_chap, "03_Data/03_Results/99_MovementModel.rds"))
mod <- mod$Model[[1]]
write_rds(mod, file.path(new_chap, "03_Data/02_CleanData/MovementModel.rds"))

# We also need to copy the scaling parameters used in that model
file.copy(
    from = file.path(old_chap, "03_Data/03_Results/99_Scaling.rds")
  , to   = file.path(new_chap, "03_Data/02_CleanData/Scaling.rds")
)

# Copy the step-length distribution
file.copy(
    from = file.path(old_chap, "03_Data/03_Results/99_GammaDistribution.rds")
  , to   = file.path(new_chap, "03_Data/02_CleanData/GammaDistribution.rds")
)

# Copy map of africa
afr <- vect(file.path(old_chap, "03_Data/02_CleanData/00_General_Africa_ESRI.shp"))
writeVector(afr, file.path(new_chap, "03_Data/02_CleanData/Africa.shp"), overwrite = T)

# Copy map of kaza
kaz <- vect(file.path(old_chap, "03_Data/02_CleanData/00_General_KAZA_KAZA.shp"))
writeVector(kaz, file.path(new_chap, "03_Data/02_CleanData/KAZA.shp"), overwrite = T)

# Copy roads
roa <- vect(file.path(old_chap, "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.shp"))
roa <- crop(roa, ext)
writeVector(roa, file.path(new_chap, "03_Data/02_CleanData/Roads.shp"), overwrite = T)

# Copy major waters
maj <- vect(file.path(old_chap, "03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters_GEOFABRIK.shp"))
maj <- crop(maj, ext)
writeVector(maj, file.path(new_chap, "03_Data/02_CleanData/MajorWaters.shp"), overwrite = T)

# Copy rivers (those that Gabs provided)
riv <- vect(file.path(old_chap, "03_Data/01_RawData/GABRIELE/Rivers.shp"))
values(riv) <- data.frame(Name = values(riv)[, c("Name")])
riv <- riv[riv$Name %in% c("Tamalakane", "Mababe-Dombo", "Boteti", "Selinda", "Selinda2", "Savuti")]
riv <- crop(riv, ext)
writeVector(riv, file.path(new_chap, "03_Data/02_CleanData/MajorRivers.shp"), overwrite = T)

# Copy faults
fau <- vect(file.path(old_chap, "03_Data/01_RawData/DAVID/Faults.shp"))
writeVector(fau, file.path(new_chap, "03_Data/02_CleanData/Faults.shp"), overwrite = T)

# Copy protected areas
pro <- vect(file.path(old_chap, "03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp"))
pro <- crop(pro, ext)
writeVector(pro, file.path(new_chap, "03_Data/02_CleanData/Protected.shp"), overwrite = T)

# Copy floodmaps
dir.create("03_Data/01_RawData/FLOODMAPS", showWarnings = F)
file.copy(
    from      = dir(path = file.path(old_chap, "03_Data/02_CleanData/00_Floodmaps/01_Original"), full.names = T)
  , to        = file.path(new_chap, "03_Data/01_RawData/FLOODMAPS")
  , recursive = T
)

# Copy river layer (Merit), but cropped to a smaller extent
dir.create("03_Data/01_RawData/MERIT", showWarnings = F)
riv <- rast(file.path(old_chap, "03_Data/02_CleanData/03_LandscapeFeatures_Rivers_MERIT.tif"))
riv <- crop(riv, ext)
names(riv) <- "Rivers"
writeRaster(riv, file.path(new_chap, "03_Data/01_RawData/MERIT/Rivers.tif"), overwrite = T)

# Copy river layer (Merit), but cropped to a smaller extent
dir.create("03_Data/01_RawData/GLOBELAND", showWarnings = F)
wat <- rast(file.path(old_chap, "03_Data/02_CleanData/01_LandCover_WaterCover_GLOBELAND.tif"))
wat <- crop(wat, ext)
names(wat) <- "Water"
writeRaster(wat, file.path(new_chap, "03_Data/01_RawData/GLOBELAND/Water.tif"), overwrite = T)

# Copy human influence layer, but cropped to a smaller extent
hum <- rast(file.path(old_chap, "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluenceBuff_FACEBOOK.grd"))
hum <- hum[[nlyr(hum)]]
hum <- crop(hum, ext)
names(hum) <- "HumanInfluence"
writeRaster(hum, file.path(new_chap, "03_Data/02_CleanData/HumanInfluence.tif"), overwrite = T)

# Copy Vegetation layers, also cropped to a smaller extent
dir.create("03_Data/01_RawData/MODIS", showWarnings = F)
tre <- rast(file.path(old_chap, "03_Data/01_RawData/MODIS/MOD44B/Stitched/01_LandCover_TreeCover_MODIS.tif"))
gra <- rast(file.path(old_chap, "03_Data/01_RawData/MODIS/MOD44B/Stitched/01_LandCover_NonTreeVegetation_MODIS.tif"))
tre <- crop(tre, ext)
gra <- crop(gra, ext)
names(tre) <- "Trees"
names(gra) <- "Shrubs"
writeRaster(tre, file.path(new_chap, "03_Data/01_RawData/MODIS/Trees.tif"), overwrite = T)
writeRaster(gra, file.path(new_chap, "03_Data/01_RawData/MODIS/Shrubs.tif"), overwrite = T)

################################################################################
#### Dispersal Data
################################################################################
# Copy data of the dispersers
file.copy(
    from = file.path(old_chap, "03_Data/02_CleanData/00_General_Dispersers_POPECOL.csv")
  , to   = file.path(new_chap, "03_Data/02_CleanData/Dispersers.csv")
)

################################################################################
#### Reference Layers
################################################################################
# Print to terminal
cat("Generating reference raster and shapefile...\n")

# Let's store the extent as a reference shapefile
s <- as.polygons(ext)
crs(s) <- "+init=epsg:4326"

# Also generate a reference raster with 250m resolution
r <- rast(ext, resolution = 250 / 111000)
values(r) <- rbinom(ncell(r), size = 1, prob = 0.5)

# Store both to file
writeVector(s, "03_Data/02_CleanData/ReferenceShape.shp", overwrite = T)
writeRaster(r, "03_Data/02_CleanData/ReferenceRaster.tif", overwrite = T)

################################################################################
#### Session Information
################################################################################
# Create folder
dir.create("02_R-Scripts/99_SessionInformation", showWarnings = F)

# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_Setup.rds")

# Print to terminal
cat("Done :)\n")
