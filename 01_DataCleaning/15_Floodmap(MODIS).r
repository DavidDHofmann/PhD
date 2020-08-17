############################################################
#### Automated Floodmaps using MODIS MCD43A4 Imagery
############################################################
# Description: This script allows to download, process and classify MODIS imagery
# automatically.

# Clear R's brain
rm(list = ls())

# Load required packages. For the installation of the "getSpatialData" package,
# check here (https://github.com/16EAGLE/getSpatialData)
library(tidyverse)
library(raster)
library(gdalUtils)
library(lubridate)
library(velox)
library(getSpatialData)
library(rgeos)
library(splitstackshape)

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load custom functions
source("Functions.r")

# Make use of multiple cores
beginCluster()

############################################################
#### Set Parameters
############################################################
# Specify your login credentials for USGS
login_USGS(username = "DoDx9", password = "EarthExplorer99")

# Specify your working directory, which needs to hold the following folders:
# (1) "Masks":  Holds projected shapefiles (.shp) of the watermask and dryland mask
# (2) "Input":  Should be empty. This is where the downloaded MODIS images will go
# (3) "Output": This is where the classified images go
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"

# Name of the shapefile that indicates the extent for which you want a floodmap
extent <- "OkavangoExtent"

# Name of the dryland mask in your extent
drymask <- "MaskDryland"

# Do you want to use dynamic watermasks? Note that this will drastically
# increase the time required to obtain floodmaps
dynamic <- TRUE

# If you set the above value to FALSE, please indicate the name of the wetland
# mask for the desired extent
watmask <- "MaskWater"

# Please indicate the name of the dryland mask for the desired extent
nowatmask <- "MaskNoWater"

# Specify dates for which you want inundation maps
dates <- c(
    "2012-01-09"
  , "2012-01-17"
  , "2012-01-25"
  , "2012-02-02"
  , "2012-02-10"
  , "2015-11-17"
  , "2015-11-25"
  , "2015-12-03"
  , "2015-12-11"
  , "2015-12-19"
  , "2015-12-27"
  , "2016-01-01"
  , "2016-01-09"
  , "2016-01-17"
  , "2016-11-08"
  , "2016-11-16"
  , "2016-11-24"
  , "2016-12-02"
  , "2017-03-14"
  , "2017-03-22"
  , "2017-11-01"
  , "2017-11-09"
  , "2017-11-17"
  , "2017-11-25"
  , "2017-12-03"
  , "2017-12-11"
  , "2017-12-19"
  , "2017-12-27"
  , "2018-01-01"
  , "2018-01-09"
  , "2018-01-17"
  , "2018-03-06"
  , "2018-03-14"
  , "2018-03-22"
  , "2018-03-30"
  , "2018-04-07"
  , "2018-04-15"
  , "2018-04-23"
  , "2018-05-01"
  , "2018-05-09"
  , "2018-05-17"
  , "2018-05-25"
  , "2018-08-29"
  , "2018-09-06"
  , "2018-09-14"
  , "2018-09-22"
  , "2018-09-30"
  , "2018-10-08"
  , "2018-10-16"
  , "2018-10-24"
  , "2018-11-01"
  , "2018-11-09"
  , "2018-11-17"
  , "2018-11-25"
  , "2018-12-03"
  , "2018-12-11"
  , "2018-12-19"
  , "2019-01-17"
  , "2019-01-25"
  , "2019-02-18"
  , "2019-02-26"
  , "2019-03-06"
  , "2019-03-14"
  , "2019-03-22"
  , "2019-03-30"
  , "2019-04-07"
  , "2019-04-15"
  , "2019-04-23"
  , "2019-05-01"
  , "2019-05-17"
  , "2019-05-25"
  , "2019-06-02"
  , "2019-06-10"
  , "2019-06-18"
  , "2019-06-26"
  , "2019-07-04"
  , "2019-07-12"
  , "2019-07-20"
  , "2019-07-28"
  , "2019-08-05"
  , "2019-08-13"
  , "2019-08-21"
  , "2019-08-29"
  , "2019-09-06"
  , "2019-09-14"
  , "2019-09-22"
)

############################################################
#### For Validation Purposes Only
############################################################
# # This part contains some code that allows to compare the results of our own
# # classification algorithm against the results of the original ORI algorithm
#
# set.seed(1234)
# setwd("/home/david/ownCloud/University/15. PhD/00_WildDogs/03_Data/01_RawData/ORI/01_Maps")
# dates <- dir(pattern = ".*tif$") %>%
#   substr(., start = 1, stop = 10) %>%
#   as.Date()
#
# # Extract the months
# dates <- data.frame(Date = dates, Month = as.factor(month(dates)))
#
# # Randomly draw 48 dates (4 x 12) and make sure that each month has equal
# # representation
# dates <- stratified(dates, "Month", 4)[[1]] %>% as.character()
# dates

############################################################
#### Automated Download and Preparation of Data
############################################################
# Define directory in which your masks are located
masks <- paste0(wd, "/03_Data/01_RawData/MODIS/MCD43A4/02_Masks")

# Define directory in which downloaded MODIS files shoud go
input <- paste0(wd, "/03_Data/01_RawData/MODIS/MCD43A4/01_Maps")

# Define directory in which classified images should go
output <- paste0(wd, "/03_Data/02_CleanData/00_Floodmaps/01_Original")
# output <- paste0(wd, "/03_Data/02_CleanData/00_Floodmaps/03_Validation") # FOR VALIDATION ONLY. DELETE FOR ORI

# Identify available MODIS products
products <- getMODIS_names()

# We are interested in the dataset MCD43A4
product <- products[grep(pattern = "MCD43A4", products)]

# Load the mask for the area of investigation (AOI)
setwd(masks)
reference <- shapefile(paste0(extent, ".shp"))

# Tell the MODIS download package the area of interest
set_aoi(reference)

# Loop through the dates, download and prepare the data
setwd(input)
for (h in 1:length(dates)){

  # Prepare the MODIS query
  query <- getMODIS_query(
      time_range = c(dates[h], dates[h])
    , name = product
  )

  # Subset to desired acquisition date. Note that I subset by the day (3 digit
  # number) in the displayId because only then I am able to get the same dates as
  # are shown on ORI's website
  sub <- subset(query, substr(query$displayId, start = 14, stop = 16) ==
    sprintf("%03d", yday(dates[h])))

  # Download the hdf tiles
  data <- getMODIS_data(records = sub, dir_out = getwd())

  # Go through the tiles and extract the bands
  for (i in 1:length(data)){

    # Identify the bands
    files <- get_subdatasets(data[i])

    # Select the bands we want to keep
    files <- files[grep("Nadir.*Band[1-7]", files)]

    # Assign a tiff name to each band
    names <- paste0(
      substr(files, start = nchar(files)-65, stop = nchar(files)), ".tif"
    )

    # Convert the selected bands to tifs using the new names
    for (i in 1:length(files)){
      gdal_translate(files[i], dst_dataset = names[i])
    }

    # Load the bands and stack them into one file
    r <- stack(names)

    # Store the stacked file
    writeRaster(
        r
      , filename  = paste0("Merged", names[1])
      , format    = "GTiff"
      , overwrite = TRUE
      , options   = c("INTERLEAVE = BAND", "COMPRESS = LZW")
    )

    # Remove the seperate bands
    file.remove(names)
  }

  # Now create a virtual raster with the tiles so we can merge them
  stitch <- dir(pattern = "Merged")
  gdalbuildvrt(gdalfile = stitch, output.vrt = "00_Merge.vrt")

  # Coerce the virtual raster to a true raster.
  gdal_translate(
    src_dataset   = "00_Merge.vrt",
    dst_dataset   = paste0(dates[h], ".tif"),
    output_Raster = TRUE,
    options       = c("BIGTIFFS=YES")
  )

  # Load the raster, reproject it to WGS84, crop it and save it
  r <- stack(paste0(dates[h], ".tif"))
  r <- projectRaster(r, crs = CRS("+init=epsg:4326"), method = "ngb")
  reference <- spTransform(reference, crs(r))
  r <- crop(r, reference)
  writeRaster(
      r
    , filename = paste0(dates[h], ".tif")
    , format = "GTiff"
    , overwrite = TRUE
    , options = c("INTERLEAVE = BAND", "COMPRESS = LZW")
  )

  # Remove the virtual raster
  file.remove("00_Merge.vrt")

  # Remove the merged layers
  file.remove(dir(pattern = "Merged"))

  # Remove the MCD folders to save space
  unlink(dir(pattern = "MCD"), recursive = TRUE)
}

############################################################
#### Reload the Downloaded Images
############################################################
# Identify the downloaded and preprocessed images on which we want to run the
# classification (those that correspond to the dates that we specified in the
# beginning)
setwd(input)
files <- paste0(dates, ".tif")

# Load band 7 of all images that we want to classify into a stack
mod <- stack(files, bands = 7)

# Also load the watermask, drymask and nowater mask
setwd(masks)
water   <- shapefile(paste0(watmask, ".shp")) %>% aggregate(., dissolve = TRUE)
dryland <- shapefile(paste0(drymask, ".shp"))
nowater <- shapefile(paste0(nowatmask, ".shp"))

############################################################
#### Function to Create Dynamic Watermask
############################################################
# Maybe we want a dynamic wetmask that adapts to changes of the inundation
# extent. To do so, we use all classified images from the past 5 years and
# identify areas that have been constantly flooded. For new regions for which no
# such data is available, one could initialize the algorithm with a shapefile of
# known waterareas, which gradually looses power as more and more images are
# classified.
watMask <- function(x){

  # Identify the last date which we would include to calculate the mask
  end_date <- names(x) %>%
    substr(., start = 2, stop = 11) %>%
    as.Date(., format = "%Y.%m.%d")

  # Subtract 5 years, to get the first date we would include to calculate the mask
  start_date <- end_date - years(5)

  # Identify all possible dates between start and end dates for which we would
  # include maps to calculate the mask
  period <- seq(start_date, end_date, "days")

  # Set the working directory to the path in which all classified images are
  # stored
  setwd(output)

  # Identify the names of all classified images
  files <- dir(pattern = "*.tif$")

  # Identify the dates of all classified images
  filedate <- as.Date(files, format = "%Y.%m.%d")

  # Keep only those filenames which are within the period of interest
  files <- files[filedate %in% period]

  # Load the files into a stack
  formask <- stack(files, bands = 1)

  # Reclassify the stack so that water becomes 1, dryland and clouds 0
  rcl <- data.frame(old = c(0, 127, 255), new = c(1, 0, 0))
  formask <- mcreclassify(formask, rcl)

  # Sum the layers
  sum <- sum(formask)

  # Identify areas where there was water 99% of the time
  wetmask <- rasterToPolygons(sum
    , fun = function(x){x > 0.99 * nlayers(formask)}
    , dissolve = TRUE
  )

  # Apply a small negative buffer to avoid errors due to pixel size
  wetmask <- gBuffer(wetmask, width = -1/111*0.25)

  # Return the final watermask
  return(wetmask)
}

# Run the function on a selected layer
# watMask(mod[[1]])

############################################################
#### Compare Spectral Signatures
############################################################
# Prepare a function that allows to plot histograms of the spectral signatures
# of dryland and wetland
showSpecs <- function(x, water, dryland){

    # Make sure that the water and dryland masks have the same crs as the modis
    # layer
    water   <- spTransform(water, crs(x))
    dryland <- spTransform(dryland, crs(x))

    # Extract the spectral values of band 7 below the dryland and water polygons
    # To speed up the extraction we coerce the modis band to a velox raster
    band7 <- velox(x)

    # Extract the values and coerce the output to a dataframe
    wat <- band7$extract(water) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    dry <- band7$extract(dryland) %>%
      do.call(rbind, .) %>%
      as.data.frame()

    # Prepare a column that indicates the land cover class
    wat$Class <- "Water"
    dry$Class <- "Dryland"

    # Bind the extracted values together
    specs <- rbind(wat, dry)

    # Plot the two densities for the spectral signatures of each value
    ggplot(specs, aes(V1, fill = Class)) + geom_density(alpha = 0.2)
}

# Run the function on a selected layer
# showSpecs(mod[[3]], water = water, dryland = dryland)

############################################################
#### Bimodality check
############################################################
# Write a function that retrieves the percentiles of the spectral reflectances
# of dryland and water
percentiles <- function(x, water, dryland){

  # Make sure that the water and dryland masks have the same crs as the modis
  # layer
  water   <- spTransform(water, crs(x))
  dryland <- spTransform(dryland, crs(x))

  # Extract the spectral values of band 7 below the dryland and water polygons.
  # To speed up the extraction we coerce the modis band to a velox raster
  band7 <- velox(x)

  # Extract the values and coerce the output to a dataframe
  wat <- band7$extract(water) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  dry <- band7$extract(dryland) %>%
    do.call(rbind, .) %>%
    as.data.frame()

  # Calculate the percentiles
  wat_perc <- quantile(wat[, 1], probs = c(0.01, 0.99), na.rm = TRUE)
  dry_perc <- quantile(dry[, 1], probs = c(0.01, 0.99), na.rm = TRUE)

  # Prepare a dataframe
  percs <- data.frame(WaterPercentiles = wat_perc, DrylandPercentiles = dry_perc)

  # Define output
  return(percs)
}

# Run the function on a selected layer
# percentiles(mod[[3]], water = water, dryland = dryland)

############################################################
#### Classify MODIS images
############################################################
# Write function that classifies an image
classify <- function(image, water, dryland){

  # Check if the image is bimodal
  perc <- percentiles(image, water = water, dryland = dryland)
  bimodal <- perc$WaterPercentiles[2] - 10 / 255 < perc$DrylandPercentiles[1]

  # If the image is not bimodal, we can return a message and skip the rest
  if (!bimodal){
    return("Image not bimodal. Could not classify floodmap.")
  } else {

    # Make sure that the water and dryland masks have the same crs as the modis
    # layer
    water   <- spTransform(water, crs(image))
    dryland <- spTransform(dryland, crs(image))
    nowater <- spTransform(nowater, crs(image))

    # Extract the spectral values of band 7 below the dryland and water
    # polygons. To speed up the extraction we coerce the MODIS band to a velox
    # raster
    band7 <- velox(image)
    wat <- band7$extract(water) %>%
      do.call(rbind, .) %>%
      as.vector() %>%
      median(na.rm = TRUE)
    dry <- band7$extract(dryland) %>%
      do.call(rbind, .) %>%
      as.vector() %>%
      median(na.rm = TRUE)

    # Calculate the classification threshold
    mu <- wat + 0.3 * (dry - wat)

    # Predict/Classify the MODIS image into water (1) and dryland (0) using the
    # above calculated threshold (Clouds will become NA anyways)
    pred <- reclassify(image, c(0,mu,1, mu,1,0))

    # Get rid of areas that are most likely misclassified or only small ponds
    # that are irrelevant for our study. Set the raster values below this
    # polygon to 0
    pred <- raster::mask(pred, nowater, updatevalue = 0, inverse = TRUE)

    # Write the raster to a temporary file
    pred <- writeRaster(pred, tempfile())

    # Return the classified image
    return(pred)
  }
}

# Run the function for the desired dates
for (i in 1:length(dates)){

  # If dynamic, create new watermask
  if (dynamic){
    water <- watMask(mod[[i]])
  }

  # Classify raster
  classified <- classify(mod[[i]], water = water, dryland = dryland)

  # If the classification succeeded reclassify and store the map
  if (class(classified) == "RasterLayer"){

    # Reclassify the values so they match ORIs classification scheme (0 = water,
    # 127 = cloud, 255 = dryland)
    rcl <- data.frame(old = c(0, 1, NA), new = c(255, 0, 127))
    flood <- reclassify(classified, rcl)

    # Write the classified image to file. The name of the file will be the date
    # of the image
    filename <- names(flood) %>%
      substr(start = 2, stop = 11) %>%
      paste0(., ".tif")
    setwd(output)
    writeRaster(flood, filename, overwrite = TRUE)
  }

  # Print status of the loop
  cat(i, "of", length(dates), "done...\n")
}

# Terminate the cluster
endCluster()
