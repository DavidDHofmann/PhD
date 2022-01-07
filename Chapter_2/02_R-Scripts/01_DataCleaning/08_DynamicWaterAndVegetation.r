################################################################################
#### Dynamic Vegetation Maps
################################################################################
# Description: Use floodmaps to mask out anything that is covered by water in
# the vegetation maps. Specifically, we will for each floodmap create a
# corresponding vegetation map where water has been masked out. This requires us
# to match the date of the floodmaps with the date of the vegetation maps. Note:
# vegetation maps are based on data from day 65 of one year to day 64 of the
# next year -> From the MODIS documentation
# (https://lpdaac.usgs.gov/documents/112/MOD44B_User_Guide_V6.pdf): The product
# date (of the vegetation data) refers to the start date of the annual period so
# a product with ID “2005065” was produced with data from 2005065 – 2006064. The
# start date of all MODIS VCF products is yyyy065 (where yyyy refers to the 4
# digit year). This originally derived from the first full 16-day composite
# period in the MODIS data record which begins with day of year 2000065. Hence,
# the layername of each MODIS vegetation map refers to the period that starts on
# day 65 of that year!

# Clear R's brain
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(terra)       # To handle raster data
library(lubridate)   # To handle dates
library(tidyverse)   # To wrangle data
library(pbmcapply)   # For multicore abilities with progress bar

# Load the layers we want to merge
flood <- rast("03_Data/02_CleanData/01_LandCover_WaterCoverDynamic.grd")
veget <- "03_Data/02_CleanData/00_Vegmaps" %>%
  dir(pattern = ".grd$", full.names = T) %>%
  rast()

# Separate shrub and tree cover
trees <- veget[[grepl(names(veget), pattern = "TreeCover")]]
shrub <- veget[[grepl(names(veget), pattern = "ShrubCover")]]

# Identify dates from layernames
flood_dates <- names(flood) %>%
  ymd()
tree_dates <- names(trees) %>%
  substr(start = 11, stop = 21) %>%
  ymd()
shrub_dates <- names(shrub) %>%
  substr(start = 12, stop = 22) %>%
  ymd()

# Should be the same for tree and shrubs
stopifnot(tree_dates == shrub_dates)
veget_dates <- tree_dates
rm(tree_dates, shrub_dates, veget)

# Put all into a nice tibble
veget_dat <- tibble(
    DateFrom   = veget_dates
  , DateTo     = as.Date(64, origin = paste0(year(DateFrom) + 1, "-01-01")) - days(1)
  , TreeCover  = lapply(1:nlyr(trees), function(x){trees[[x]]})
  , ShrubCover = lapply(1:nlyr(shrub), function(x){shrub[[x]]})
)
print(veget_dat)

# We extend the last date to today because we don't have any updated vegetation
# data for this
veget_dat$DateTo[veget_dat$DateTo == "2021-03-05"] <- today()

# Prepare a tibble for the flood cover as well
flood_dat <- tibble(
    Date       = flood_dates
  , FloodCover = lapply(1:nlyr(flood), function(x){flood[[x]]})
)
print(flood_dat)

# Combine floodmaps with the vegetation maps
cat("Combining floodmaps with vegetation maps...\n")
veget <- lapply(1:nrow(flood_dat), function(x) {

  # Identify into which period of the vegetation maps the current floodmap falls
  date <- flood_dat$Date[x]
  index <- which(date >= veget_dat$DateFrom & date <= veget_dat$DateTo)

  # Merge the floodmap with the respective vegetation maps
  flood_map <- flood_dat$FloodCover[[x]]
  trees_map <- veget_dat$TreeCover[[index]]
  shrub_map <- veget_dat$ShrubCover[[index]]
  trees_map <- mask(trees_map, flood_map, maskvalue = 1, updatevalue = 0)
  shrub_map <- mask(shrub_map, flood_map, maskvalue = 1, updatevalue = 0)

  # Put everything together
  veget_map <- c(trees_map, shrub_map)
  names(veget_map) <- paste(c("Trees_", "Shrubs_"), date)
  return(veget_map)
})
veget <- rast(veget)

# Separate tree cover from shrub cover
trees_maps <- veget[[grepl(names(veget), pattern = "Trees")]]
shrub_maps <- veget[[grepl(names(veget), pattern = "Shrubs")]]

# Store the layers to file
cat("Storing data...\n")
writeRaster(trees_maps, "03_Data/02_CleanData/01_LandCover_TreeCoverDynamic.grd", overwrite = T)
writeRaster(shrub_maps, "03_Data/02_CleanData/01_LandCover_ShrubCoverDynamic.grd", overwrite = T)

# We also want to create a corresponding set of maps for the "static" watermap
cat("Creating static vegetation maps...\n")
flood_static <- rast("03_Data/02_CleanData/01_LandCover_WaterCoverStatic.tif")
trees_static <- mask(mean(trees), flood_static, maskvalue = 1, updatevalue = 0)
shrub_static <- mask(mean(shrub), flood_static, maskvalue = 1, updatevalue = 0)

# Store the layers to file
cat("Storing data...\n")
writeRaster(trees_static, "03_Data/02_CleanData/01_LandCover_TreeCoverStatic.tif", overwrite = T)
writeRaster(shrub_static, "03_Data/02_CleanData/01_LandCover_ShrubCoverStatic.tif", overwrite = T)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/08_DynamicWaterAndVegetation.rds")
cat("Done :)\n")
