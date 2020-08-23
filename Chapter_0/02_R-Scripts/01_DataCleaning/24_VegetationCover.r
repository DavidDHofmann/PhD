############################################################
#### Preparation of the Modis Vegetation Data
############################################################
# Description: Stitching, resampling and preparation of the modis continuous
# vegetation maps. Note that the MODIS treecover maps are dated the 6 March
# 2017.

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_0"
setwd(wd)

# load packages
library(tidyverse)
library(raster)
library(gdalUtils)
library(parallel)
library(terra)

# Make use of multicore abilities
beginCluster()

############################################################
#### Extracting TIFF files
############################################################
# Identify files as they were downloaded
files <- dir(
    path        = "03_Data/01_RawData/MODIS/MOD44B"
  , pattern     = "hdf"
  , full.names  = T
)

# Run a loop that takes each hdf dataset and extracts the layers of interest
# (Percent Tree Cover, Percent Non-Tree Cover, Percent-Non-Vegetated)
for (i in 1:length(files)){
  subs  <- get_subdatasets(files[i])
  for (j in 1:3){
    part1 <- substr(subs, start = 67, stop = 72)
    part2 <- substr(subs, start = 113, stop = nchar(subs))
    part3 <- ".tif"
    names <- paste0(part1, part2, part3)
    names <- paste0(wd, "/03_Data/01_RawData/MODIS/MOD44B/", names)
    gdal_translate(subs[j], dst_dataset = names[j])
  }
}

############################################################
#### Stitching TIFFs together
############################################################
# Load the reference raster
r250 <- rast("03_Data/02_CleanData/00_General_Raster250.tif")

# Specify layers that need to be stitched
patterns <- c(
    "Tree_Cover"
  , "NonVegetated"
  , "NonTree_Vegetation"
)

# We also prepare the names of the merged layers
newnames <- c(
    "01_LandCover_TreeCover_MODIS.tif"
  , "01_LandCover_NonVegetated_MODIS.tif"
  , "01_LandCover_NonTreeVegetation_MODIS.tif"
)
newnames <- paste0("03_Data/01_RawData/MODIS/MOD44B/Stitched/", newnames)

# Create directory
dir.create("03_Data/01_RawData/MODIS/MOD44B/Stitched", showWarnings = F)

# Then loop through the three layers
for (i in 1:length(patterns)){

  # Load the tiles that correspond to the same layer type
  files <- dir(
      path        = "03_Data/01_RawData/MODIS/MOD44B"
    , pattern     = patterns[i]
    , full.names  = T
  )

  # Create a virtual raster using the tiles
  gdalbuildvrt(gdalfile = files, output.vrt = "00_Merge.vrt")

  # Coerce the virtual raster to a true raster
  gdal_translate(
    src_dataset   = "00_Merge.vrt",
    dst_dataset   = newnames[i],
    output_Raster = TRUE,
    options       = c("BIGTIFFS=YES")
  )

  # Remove the virtual raster
  file.remove("00_Merge.vrt")

  # Reproject, resample and store the merged files
  merged <- rast(newnames[i]) %>%
    terra::project(., r250) %>%
    terra::resample(., r250, "near")
  writeRaster(merged, newnames[i], overwrite = TRUE)

  # Print status
  cat(i, "of", length(patterns), "done...\n")
}

# Remove files that are not needed anymore
files <- dir(
    path        = "03_Data/01_RawData/MODIS/MOD44B"
  , pattern     = ".tif$"
  , full.names  = T
) %>% file.remove()

############################################################
#### Remove Water Areas
############################################################
# Load the vegetation layers again
files <- dir(
    path        = "03_Data/01_RawData/MODIS/MOD44B/Stitched"
  , pattern     = ".*MODIS.tif$"
  , full.names  = T
)
names <- substr(files, start = 55, stop = nchar(files) - 10)
modis <- rast(files)
names(modis) <- names

# Extract the separate layers
modis_shrub <- modis[[1]]
modis_noveg <- modis[[2]]
modis_trees <- modis[[3]]

# Values above 100 are water. Let's reclassify those to 0% Vegetation, i.e. 100%
# NonVegetated
values(modis_shrub)[values(modis_shrub) > 100] <- 0
values(modis_noveg)[values(modis_noveg) > 100] <- 100
values(modis_trees)[values(modis_trees) > 100] <- 0

# Visualize again
plot(c(modis_shrub, modis_noveg, modis_trees))

# Add up all of the rasters
summed <- sum(modis_shrub, modis_noveg, modis_trees)

# Check if summed values add up to 100
hist(summed)

# Load the Globeland Layer
glo <- rast("03_Data/02_CleanData/01_LandCoverClasses30_Globeland.tif")

# Also load the dynamic water layers
water <- rast("03_Data/02_CleanData/01_LandCover_Water(Merged).tif")

# Identify all dates for which we have a closest watermap
dates_globe <- names(water) %>%
  substr(start = 2, stop = nchar(.)) %>%
  as.Date(format = "%Y.%m.%d")

# For the same dates we now want to create dynamic vegetation layers. However,
# in contrast to the merged globeland layer we now don't want to rasterize
# rivers again.
files <- dir(
    path        = "03_Data/02_CleanData/00_Floodmaps/02_Resampled/"
  , full.names  = T
)

# Get their dates
dates_ori <- files %>%
  substr(start = 49, stop = 58) %>%
  as.Date("%Y.%m.%d")

# Subset to relevant dates
files <- files[dates_ori %in% dates_globe]

# Load those dynamice water layers
ori <- stack(files)

# Convert globeland to raster
glo <- raster(glo)

# Extend floodmaps so that their extent matches the one of the globeland layer
ori <- suppressMessages(
  mclapply(1:nlayers(ori), mc.cores = detectCores() - 1, function(x){
    extended <- extend(ori[[x]], glo, value = NA)
    extended <- writeRaster(extended, tempfile())
    return(extended)
  }) %>% stack()
)

# Cover missing values in the ori layers with globeland data
ori <- suppressMessages(
  mclapply(1:nlayers(ori), mc.cores = detectCores() - 1, function(x){
    covered <- cover(ori[[x]], glo)
    covered <- writeRaster(covered, tempfile())
    return(covered)
  }) %>% stack()
)

# Mask water in shrub layer
modis_shrub <- suppressMessages(
  mclapply(1:nlayers(ori), mc.cores = detectCores() - 1, function(x){
    masked <- mask(modis_shrub
      , mask        = ori[[x]]
      , maskvalue   = 1
      , updatevalue = 0
    )
    masked <- writeRaster(masked, tempfile())
    return(masked)
  }) %>% stack()
)

# Mask water in nonvegetated layer
modis_noveg <- suppressMessages(
  mclapply(1:nlayers(ori), mc.cores = detectCores() - 1, function(x){
    masked <- mask(modis_noveg
      , mask        = ori[[x]]
      , maskvalue   = 1
      , updatevalue = 100
    )
    masked <- writeRaster(masked, tempfile())
    return(masked)
  }) %>% stack()
)

# Mask water in trees layer
modis_trees <- suppressMessages(
  mclapply(1:nlayers(ori), mc.cores = detectCores() - 1, function(x){
    masked <- mask(modis_trees
      , mask        = ori[[x]]
      , maskvalue   = 1
      , updatevalue = 0
    )
    masked <- writeRaster(masked, tempfile())
    return(masked)
  }) %>% stack()
)

# Convert to terra again
modis_shrub <- rast(modis_shrub)
modis_noveg <- rast(modis_noveg)
modis_trees <- rast(modis_trees)

# Put dates as layer names
names(modis_shrub) <- dates_ori
names(modis_noveg) <- dates_ori
names(modis_trees) <- dates_ori

# Visualize layers again
plot(c(modis_shrub[[1]], modis_noveg[[1]], modis_trees[[1]]))

############################################################
#### Store the Final Maps
############################################################
# Make sure that the layers range from 0 to 1 rather than from 0 to 100
modis_shrub <- modis_shrub / 100
modis_noveg <- modis_noveg / 100
modis_trees <- modis_trees / 100

# Visualize layers again
plot(c(modis_shrub[[1]], modis_noveg[[1]], modis_trees[[1]]))

# Put the stacks into a list
modis <- list(modis_shrub, modis_noveg, modis_trees)

# Prepare filenames
names <- c(
    "03_Data/02_CleanData/01_LandCover_NonTreeVegetation_MODIS.tif"
  , "03_Data/02_CleanData/01_LandCover_NonVegetated_MODIS.tif"
  , "03_Data/02_CleanData/01_LandCover_TreeCover_MODIS.tif"
)

# Store the rasterstacks
for (i in 1:length(names)){
  terra::writeRaster(modis[[i]], names[i], overwrite = TRUE)
}

# Terminate cluster
endCluster()
