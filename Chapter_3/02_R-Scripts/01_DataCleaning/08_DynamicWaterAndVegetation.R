################################################################################
#### Dynamic Vegetation Maps
################################################################################
# Description: Use floodmaps to mask out anything that is covered by water in
# the vegetation maps. Specifically, we will, for each floodmap, create a
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
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(terra)       # To handle raster data
library(lubridate)   # To handle dates
library(tidyverse)   # To wrangle data
library(pbapply)     # For progress bar lapply

# Load the layers we want to merge
flood <- rast("03_Data/02_CleanData/WaterDynamic.tif")
veget <- "03_Data/02_CleanData/00_Vegmaps" %>%
  dir(pattern = ".tif$", full.names = T) %>%
  rast()

# Separate shrub and tree cover
trees <- veget[[grepl(names(veget), pattern = "TreeCover")]]
shrub <- veget[[grepl(names(veget), pattern = "ShrubCover")]]

# Identify dates from layernames
flood_dates <- names(flood) %>%
  ymd()
tree_dates <- names(trees) %>%
  substr(start = nchar(.) - 9, stop = nchar(.)) %>%
  ymd()
shrub_dates <- names(shrub) %>%
  substr(start = nchar(.) - 9, stop = nchar(.)) %>%
  ymd()

# Should be the same for tree and shrubs
stopifnot(tree_dates == shrub_dates)
veget_dates <- tree_dates
rm(tree_dates, shrub_dates, veget)

# Put all into a nice tibble
veget_dat <- tibble(
    DateFrom   = veget_dates
  , DateTo     = as.Date(64, origin = paste0(year(DateFrom) + 1, "-01-01")) - days(1)
  , Trees      = lapply(1:nlyr(trees), function(x){trees[[x]]})
  , Shrubs     = lapply(1:nlyr(shrub), function(x){shrub[[x]]})
)

# We extend the last date to today because we don't have any updated vegetation
# data for this
veget_dat$DateTo[veget_dat$DateTo == max(veget_dat$DateTo)] <- today()

# Prepare a tibble for the flood cover as well
flood_dat <- tibble(
    Date       = flood_dates
  , FloodCover = lapply(1:nlyr(flood), function(x){flood[[x]]})
)
print(flood_dat)

# Combine floodmaps with the vegetation maps
cat("Combining floodmaps with vegetation maps...\n")
pb <- txtProgressBar(min = 0, max = nrow(flood_dat), style = 3)
veget <- pblapply(1:nrow(flood_dat), function(x) {

  # Identify into which period of the vegetation maps the current floodmap falls
  date  <- flood_dat$Date[x]
  index <- which(date >= veget_dat$DateFrom & date <= veget_dat$DateTo)

  # Merge the floodmap with the respective vegetation maps
  flood_map <- flood_dat$FloodCover[[x]]
  trees_map <- veget_dat$Trees[[index]]
  shrub_map <- veget_dat$Shrubs[[index]]
  trees_map <- mask(trees_map, flood_map, maskvalue = 1, updatevalue = 0)
  shrub_map <- mask(shrub_map, flood_map, maskvalue = 1, updatevalue = 0)

  # Put everything together
  veget_map <- c(trees_map, shrub_map)
  names(veget_map) <- paste(c("Trees_", "Shrubs_"), date)
  setTxtProgressBar(pb, value = x)
  return(veget_map)
})
veget <- rast(veget)

# Separate tree cover from shrub cover
trees_maps <- veget[[grepl(names(veget), pattern = "Trees")]]
shrub_maps <- veget[[grepl(names(veget), pattern = "Shrubs")]]

# Make the layernames a bit nicer again
names(trees_maps) <- names(trees_maps) %>%
  substr(start = nchar(.) - 9, stop = nchar(.)) %>%
  ymd()
names(shrub_maps) <- names(trees_maps) %>%
  substr(start = nchar(.) - 9, stop = nchar(.)) %>%
  ymd()

# Store the layers to file
cat("Storing data...\n")
writeRaster(trees_maps
  , filename  = "03_Data/02_CleanData/TreesDynamic.tif"
  , overwrite = T
  , gdal      = "COMPRESS=NONE"
)
writeRaster(shrub_maps
  , filename  = "03_Data/02_CleanData/ShrubsDynamic.tif"
  , overwrite = T
  , gdal      = "COMPRESS=NONE"
)

# We also want to create a corresponding set of maps for the "static" watermap
cat("Creating static vegetation maps...\n")
flood_static <- rast("03_Data/02_CleanData/WaterStatic.tif")
trees_static <- mask(mean(trees), flood_static, maskvalue = 1, updatevalue = 0)
shrub_static <- mask(mean(shrub), flood_static, maskvalue = 1, updatevalue = 0)

# Store the layers to file
cat("Storing data...\n")
writeRaster(
    trees_static
  , filename  = "03_Data/02_CleanData/TreesStatic.tif"
  , overwrite = T
  , gdal      = "COMPRESS=NONE"
)
writeRaster(
    shrub_static
  , filename  = "03_Data/02_CleanData/ShrubsStatic.tif"
  , overwrite = T
  , gdal      = "COMPRESS=NONE"
)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/08_DynamicWaterAndVegetation.rds")
cat("Done :)\n")
