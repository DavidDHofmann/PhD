############################################################
#### Preparation of Globeland Data
############################################################
# Description: The script stitches together the seperated tiles downloaded from
# Globeland. After stitching the resulting raster is cropped, simplified
# (reclassified) and resampled to 250m.

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load packages
library(tidyverse)
library(raster)
library(rgdal)
library(gdalUtils)
library(parallel)
library(davidoff)

# Start cluster
beginCluster()

############################################################
#### Stitching the Tiles
############################################################
# Identify all the tifs in our folder
files <- dir(
    path        = "03_Data/01_RawData/GLOBELAND"
  , pattern     = "030.tif$"
  , full.names  = T
)

# Before we can merge the rasters they need to be reprojected to a uniform
# projection
new_files <- c()
for (i in 1:length(files)){
  new_files[i] <- tempfile(fileext = ".tif")
  gdalwarp(
      srcfile = files[i]
    , dstfile = new_files[i]
    , t_srs   = CRS("+init=epsg:4326")
  )
}

# Put the reprojected files into a virtual raster
new_vrt <- tempfile(fileext = ".vrt")
gdalbuildvrt(
    gdalfile    = new_files
  , srcnodata   = 0
  , output.vrt  = new_vrt
)

# Store the virtual raster to file
gdal_translate(
    src_dataset   = new_vrt
  , dst_dataset   = "03_Data/01_RawData/GLOBELAND/Globeland.tif"
  , output_Raster = TRUE
  , options       = c("BIGTIFFS=YES")
)

############################################################
#### Cropping the Stitched Raster
############################################################
# Load the merged file
merged <- raster("03_Data/01_RawData/GLOBELAND/Globeland.tif")

# Load the reference shapefile
s <- shapefile("03_Data/02_CleanData/00_General_Shapefile")

# Crop the merged globeland tiles to our extent
merged <- crop(merged, s)

# Store the result to file
writeRaster(
    merged
  , "03_Data/01_RawData/GLOBELAND/Globeland.tif"
  , overwrite = TRUE
)

############################################################
#### Reduce Resolution to 250m
############################################################
# Load the layer again
glo <- raster("03_Data/01_RawData/GLOBELAND/Globeland.tif")

# Aggregate the layer to 250m
coarse <- aggregate(glo, fact = round(250 / 30), fun = modal)

############################################################
#### Simplify Land Cover to Water
############################################################
# We only want to use the water layer from the globeland dataset. Let's
# reclassify the values accordingly
classes <- data.frame(
    CodesOld = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  , DescriptionOld = c(
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
  )
  , CodesNew = c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0)
  , DescriptionNew = c(
      "Dryland"
    , "Dryland"
    , "Dryland"
    , "Dryland"
    , "Water"
    , "Water"
    , "Dryland"
    , "Dryland"
    , "Dryland"
    , "Dryland"
  )
)

# Prepare a reclassification matrix
rcl <- data.frame(old = classes[, 1], new = classes[, 3])

# Run the reclassification
new <- reclassify(coarse, rcl)

############################################################
#### Resampling
############################################################
# Load the reference raster
r250 <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Resample the globeland layer to the reference raster
new <- resample(new, r250, "ngb")

# Store the result
writeRaster(
    new
  , "03_Data/02_CleanData/01_LandCoverClasses30_Globeland.tif"
  , overwrite = TRUE
)

# End cluster
endCluster()
