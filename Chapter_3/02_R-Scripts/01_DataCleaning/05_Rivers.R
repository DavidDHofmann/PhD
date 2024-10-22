################################################################################
#### Preparation of River Data
################################################################################
# Description: Prepare the Merit river dataset (stitch, reproject, etc.)

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load packages
library(terra) # For handling spatial data

################################################################################
#### MERIT Rivers
################################################################################
# Identify merit river files that need to be stitched
files <- dir(
    path        = "03_Data/01_RawData/MERIT"
  , pattern     = ".tif$"
  , full.names  = T
)

# Load them as terra rasters
dat <- lapply(files, rast)

# Mosaic them.
cat("Mosaicing MERIT river data...\n")
riv <- terra::mosaic(dat[[1]], dat[[2]], dat[[3]], dat[[4]], dat[[5]]
  , dat[[6]], dat[[7]], dat[[8]], dat[[9]], fun = "max")

# Reclassify pixel values and keep only those rivers with a desired width
rcl <- data.frame(from = c(-Inf, 10), to = c(10, Inf), new = c(0, 1))
riv <- classify(riv, rcl)

# Load reference raster
r <- rast("03_Data/02_CleanData/Raster.tif")

# Aggregate the layer to match the resolution of the reference raster
fact <- res(r)[1] / res(riv)[1]
riv <- aggregate(riv, fact = round(fact), fun = max)

# Resample the river layer to match the origin and extent of the reference
# raster
riv <- resample(riv, r, "near")

# Store the final result to file
writeRaster(riv, "03_Data/02_CleanData/Rivers.tif", overwrite = TRUE)

################################################################################
#### Some Rivers for Plotting
################################################################################
# Load shapefiles and extract rivers of interest
library(tidyterra)
riv1 <- "03_Data/01_RawData/DAVID/Rivers.shp" %>%
  vect() %>%
  tidyterra::select(Name = NAME) %>%
  filter(Name == "Okavango")
riv2 <- "03_Data/01_RawData/DAVID/Rivers2.shp" %>%
  vect() %>%
  tidyterra::select(Name) %>%
  filter(Name %in% c("Tamalakane", "Mababe-Dombo", "Boteti", "Selinda", "Selinda2", "Savuti"))
riv <- rbind(riv1, riv2)
writeVector(riv, "03_Data/02_CleanData/MajorRivers.gpkg", overwrite = T)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/05_Rivers.rds")
cat("Done :)\n")
