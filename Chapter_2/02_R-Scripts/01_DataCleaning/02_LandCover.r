################################################################################
#### Preparation of Globeland Land Cover Map
################################################################################
# Description: Simplify and resample of the Globeland Land Cover Dataset

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load packages
library(terra)       # To handle raster data
library(tidyverse)   # For data wrangling

# Load the globeland land cover map
glob <- rast("03_Data/01_RawData/GLOBELAND/Globeland.tif")

# Load the reference shapefile and raster
s <- vect("03_Data/02_CleanData/00_General_Shapefile.shp")
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")

# Crop the merged tiles to our reference shapefile
glob <- crop(glob, s)

# Aggregate to coarser resolution
fact <- res(r)[1] / res(glob)[1]
glob <- terra::aggregate(glob, fact = round(fact), fun = "modal")

# Check out the distribution of values
freq(glob)

# Prepare reclassification table
info <- data.frame(
    Code = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 255)
  , Class = c(
      "CultivatedLand"
    , "Forest"
    , "Grassland"
    , "Shrubland"
    , "Wetland"
    , "WaterBodies"
    , "Tundra"
    , "ArtificialSurfaces"
    , "Bareland"
    , "PermanentIceSnow"
    , "NA"
  )
  , CodeNew = c(3, 4, 6, 5, 1, 1, 6, 2, 7, 7, 6)
  , ClassNew = c(
      "Cropland"
    , "Forest"
    , "Grassland"
    , "Shrubs"
    , "Water"
    , "Water"
    , "Grassland"
    , "Urban"
    , "Bare"
    , "Bare"
    , "Grassland"
  )
)

# Arrange
info <- arrange(info, CodeNew)

# Assign a color to the new classes
info$Color[info$CodeNew == 1] <- "#0e16ed"
info$Color[info$CodeNew == 2] <- "#e6494b"
info$Color[info$CodeNew == 3] <- "#f8f261"
info$Color[info$CodeNew == 4] <- "#1e7b40"
info$Color[info$CodeNew == 5] <- "#cca12b"
info$Color[info$CodeNew == 6] <- "#85cc7a"
info$Color[info$CodeNew == 7] <- "#cccccc"

# Reclassify raster
rcl <- dplyr::select(info, c(Code, CodeNew))
glob <- classify(glob, rcl)

# Resample to reference raster
glob <- resample(glob, r, method = "near")

# Visualize results
ggplot(as.data.frame(glob, xy = T), aes(x = x, y = y, fill = as.factor(Globeland))) +
  geom_raster() +
  coord_sf() +
  scale_fill_manual(
      values = info$Color
    , breaks = info$CodeNew
    , labels = info$ClassNew
    , name   = "Cover Class"
  ) +
  theme_minimal()

# Store the raster
writeRaster(
    x         = glob
  , filename  = "03_Data/02_CleanData/01_LandCover_LandCover.tif"
  , overwrite = TRUE
)

# Also store the water-cover layer seperately
water <- glob == 1
writeRaster(
    x         = water
  , filename  = "03_Data/02_CleanData/01_LandCover_WaterCover.tif"
  , overwrite = TRUE
)

# Store the information table
info %>%
  dplyr::select(Class = ClassNew, Code = CodeNew, Color) %>%
  distinct() %>%
  write_csv("03_Data/02_CleanData/01_LandCover_LandCover.csv")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
write_rds(session, file = "02_R-Scripts/99_SessionInformation/02_LandCover_SessionInfo.rds")
