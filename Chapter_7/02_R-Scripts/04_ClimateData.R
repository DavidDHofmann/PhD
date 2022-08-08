################################################################################
#### Download and Curate Climate Data
################################################################################
# Description: Use Google Earth Engine to download hourly precipitation and
# temperature data for our study area

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
# wd <- "C:/Users/david/switchdrive/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(reticulate)   # Interface to python
library(rgee)         # Interface to google earth engine
library(terra)        # To handle spatial data
library(tidyverse)    # To wrangle data
library(lubridate)    # To handle dates
library(sf)           # To handle spatial data
library(sen2r)        # To download sentinel 2 data

# # Specify python to use (only necessary on some systems)
# Sys.setenv(RETICULATE_PYTHON = "/home/david/miniconda3/envs/rgee/bin/python")

# Specify correct python environment
# ee_install_set_pyenv(
#     py_path = "/home/david/miniconda3/envs/rgee/bin/python"
#   , py_env  = "rgee"
# )

# Make sure we have all installed for rgee
ee_check()
ee_Initialize()
# ee_clean_pyenv()
# ee_clean_credentials()

# Load the activity data
act <- read_csv("03_Data/02_CleanData/ActivityDataMoonphase.csv")

# Create an extent from the minimum and maximum x and y coordinates
xmin <- floor(range(act$x)[1])
xmax <- ceiling(range(act$x)[2])
ymin <- floor(range(act$y)[1])
ymax <- ceiling(range(act$y)[2])
ext <- ext(xmin, xmax, ymin, ymax)

# Use the extent to specify an area of interest
aoi <- ee$Geometry$Polygon(
  list(
      c(xmin, ymin)
    , c(xmin, ymax)
    , c(xmax, ymax)
    , c(xmax, ymin)
  )
)

# Instead of trying to download data for the entire range of dates at once, I
# want to download the maps by month. This will allow me to observe the download
# process better and to save files in between. Let's first identify the dates
# for which we want to download data.
dates <- act %>%
  pull(Timestamp) %>%
  range() %>%
  as.Date() %>%
  "+"(c(days(-1), days(+1)))

# Get rid of the activity data again
rm(act)
gc()

################################################################################
#### Download Precipitation Data
################################################################################
# Create a tibble that allows us to keep track of downloaded files
todownload <- seq(dates[1], dates[2], by = "day") %>%
  as_tibble() %>%
  setNames("Date") %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  group_by(Year, Month) %>%
  summarize(FirstDate = range(Date)[1], LastDate = range(Date)[2], .groups = "drop") %>%
  mutate(Filename = paste0("03_Data/02_CleanData/00_Rainmaps/", Year, "_", sprintf("%02d", Month), ".grd"))

# We only need to download files that do not exist yet
todownload <- subset(todownload, !file.exists(Filename))

# Create directory into which the downloaded maps go
dir.create("03_Data/02_CleanData/00_Rainmaps", showWarnings = F)

# Loop over the dataframe and download data
if (nrow(todownload) > 0) {
  cat("Downloading precipitation data...\n")
  for (i in 1:nrow(todownload)) {

    # Get dates
    date_from <- todownload$FirstDate[i] %>% as.character()
    date_to   <- (todownload$LastDate[i] + days(1)) %>% as.character()

    # If the date is before the "2014-03-01", we have to use a different link
    if (date_to <= "2014-03-01") {
        link <- "JAXA/GPM_L3/GSMaP/v6/reanalysis"
      } else {
        link <- "JAXA/GPM_L3/GSMaP/v6/operational"
    }

    # Query rainfall data for our study area during dispersal dates
    query_precip <- ee$
      ImageCollection(link)$
      filterDate(date_from, date_to)$
      filterBounds(aoi)$
      select("hourlyPrecipRate")$
      toBands()$
      clip(aoi)

    # Check it out
    # ee_print(query_precip)

    # Download the data
    tempfile <- tempfile(fileext = ".tif")
    ee_as_raster(query_precip, dsn = tempfile, region = aoi)

    # Also store the band names
    names_precip <- ee_print(query_precip)$img_bands_names %>%
      strsplit(., split = " ") %>%
      unlist() %>%
      data.frame(BandName = ., stringsAsFactors = F) %>%
      separate(BandName, into = c("Date", "Hour", "Variable"), sep = "_") %>%
      mutate(Date = ymd(Date), Hour = substr(Hour, start = 1, stop = 2)) %>%
      mutate(Layername = paste0(Date, "_", Hour)) %>%
      pull(Layername)

    # Load the downloaded data and adjust layernames
    precip <- rast(tempfile)
    names(precip) <- names_precip

    # Crop layer to our study area
    precip <- crop(precip, ext, snap = "out")

    # Store the file again
    writeRaster(precip, todownload$Filename[i])
  }
}

################################################################################
#### Temperature Data
################################################################################
# Create a tibble that allows us to keep track of downloaded files
todownload <- seq(dates[1], dates[2], by = "day") %>%
  as_tibble() %>%
  setNames("Date") %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  group_by(Year, Month) %>%
  summarize(FirstDate = range(Date)[1], LastDate = range(Date)[2], .groups = "drop") %>%
  mutate(Filename = paste0("03_Data/02_CleanData/00_Tempmaps/", Year, "_", sprintf("%02d", Month), ".grd"))

# We only need to download files that do not exist yet
todownload <- subset(todownload, !file.exists(Filename))

# Create directory into which the downloaded maps go
dir.create("03_Data/02_CleanData/00_Tempmaps", showWarnings = F)

# Loop over the dataframe and download data
if (nrow(todownload) > 0) {
  cat("Downloading temperature data...\n")
  for (i in 1:nrow(todownload)) {

    # Get dates
    date_from <- todownload$FirstDate[i] %>% as.character()
    date_to   <- (todownload$LastDate[i] + days(1)) %>% as.character()

    # Link to the dataset
    link <- "ECMWF/ERA5_LAND/HOURLY"

    # Query temperature data for our study area during dispersal dates
    query_temp <- ee$
      ImageCollection(link)$
      filterDate(date_from, date_to)$
      filterBounds(aoi)$
      select("temperature_2m")$
      toBands()$
      clip(aoi)

    # Check it out
    # ee_print(query_temp)

    # Download the data
    tempfile <- tempfile(fileext = ".tif")
    ee_as_raster(query_temp, dsn = tempfile, region = aoi)

    # Also store the band names
    names_temp <- ee_print(query_temp)$img_bands_names %>%
      strsplit(., split = " ") %>%
      unlist() %>%
      data.frame(BandName = ., stringsAsFactors = F) %>%
      separate(BandName, into = c("Date", "Variable", "Remove"), sep = "_") %>%
      separate(Date, into = c("Date", "Hour"), sep = "T") %>%
      mutate(Date = ymd(Date)) %>%
      mutate(Layername = paste0(Date, "_", Hour)) %>%
      pull(Layername)

    # Load the downloaded data and adjust layernames
    temp <- rast(tempfile)
    names(temp) <- names_temp

    # Crop layer to our study area
    temp <- crop(temp, ext, snap = "out")

    # Convert kelvin to celsius
    temp <- temp - 273.15

    # Store the file again
    writeRaster(temp, todownload$Filename[i])
  }
}

################################################################################
#### Cloud Data
################################################################################
# Create a tibble that allows us to keep track of downloaded files
todownload <- seq(dates[1], dates[2], by = "day") %>%
  as_tibble() %>%
  setNames("Date") %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  mutate(Product1 = "Aqua", Product2 = "Terra") %>%
  pivot_longer(Product1:Product2, names_to = "Remove", values_to = "Product") %>%
  select(-Remove) %>%
  group_by(Year, Month, Product) %>%
  summarize(FirstDate = range(Date)[1], LastDate = range(Date)[2], .groups = "drop") %>%
  mutate(Filename = paste0("03_Data/02_CleanData/00_Cloudmaps/", Year, "_", sprintf("%02d", Month), "_", Product, ".grd"))

# We only need to download files that do not exist yet
todownload <- subset(todownload, !file.exists(Filename))

# Create directory into which the downloaded maps go
dir.create("03_Data/02_CleanData/00_Cloudmaps", showWarnings = F)

# Function to extract values from specific bits
getQABits <- function(image, fromBit, toBit) {
  maskSize <- ee$Number(1)$
    add(toBit)$
    subtract(fromBit)
  mask <- ee$Number(1)$leftShift(maskSize)$
    subtract(1)
  masked <- image$rightShift(fromBit)$bitwiseAnd(mask)
  return(masked)
}

# Using getQABits we construct a single-argument function 'mod13A2_clean'
getClouds <- function(img) {
  qa_bits <- img$select("state_1km")
  qa_mask <- getQABits(qa_bits, 0, 1)$eq(1)
  return(qa_mask)
}

# Loop over the dataframe and download data
if (nrow(todownload) > 0) {

  cat("Downloading cloud cover data...\n")
  for (i in 1:nrow(todownload)) {

    # Get dates
    date_from <- todownload$FirstDate[i] %>% as.character()
    date_to   <- (todownload$LastDate[i] + days(1)) %>% as.character()

    # Link to the dataset
    if (todownload$Product[i] == "Aqua") {
        link <- "MODIS/061/MYD09GA"
      } else {
        link <- "MODIS/061/MOD09GA"
    }

    # Query temperature data for our study area during dispersal dates
    query_cloud <- ee$
      ImageCollection(link)$
      filterDate(date_from, date_to)$
      filterBounds(aoi)$
      map(getClouds)$
      toBands()$
      clip(aoi)

    # Check it out
    # ee_print(query_cloud)

    # Download the data
    tempfile <- tempfile(fileext = ".tif")
    ee_as_raster(query_cloud, dsn = tempfile, region = aoi)

    # Also store the band names
    names_cloud <- ee_print(query_cloud)$img_bands_names %>%
      strsplit(., split = " ") %>%
      unlist() %>%
      data.frame(BandName = ., stringsAsFactors = F) %>%
      separate(BandName, into = c("Year", "Month", "Day", "Remove", "Remove2"), sep = "_") %>%
      mutate(Date = ymd(paste0(Year, Month, Day))) %>%
      mutate(Layername = paste0(Date, "_", todownload$Product[i])) %>%
      pull(Layername)

    # Load the downloaded data and adjust layernames
    cloud <- rast(tempfile)

    # Reproject
    cloud <- project(cloud, "epsg:4326", method = "near")

    # Crop layer to our study area
    cloud <- crop(cloud, ext, snap = "out")

    # Store the file again
    writeRaster(cloud, todownload$Filename[i])
  }
}

################################################################################
#### Combine Datasets
################################################################################
# Rainmaps
files <- dir(path = "03_Data/02_CleanData/00_Rainmaps", pattern = ".grd$", full.names = T)
files <- lapply(files, rast)
files <- do.call(c, files)
writeRaster(files, "03_Data/02_CleanData/Precipitation.tif", overwrite = T)

# Temperature maps
files <- dir(path = "03_Data/02_CleanData/00_Tempmaps", pattern = ".grd$", full.names = T)
files <- lapply(files, rast)
files <- do.call(c, files)
writeRaster(files, "03_Data/02_CleanData/Temperature.tif", overwrite = T)

# Cloud maps (we will also coarsen them and combine aqua and terra products from
# the same date)
files <- dir(path = "03_Data/02_CleanData/00_Cloudmaps", pattern = ".grd$", full.names = T)
files <- lapply(files[1:4], rast)
files <- do.call(c, files)
dates <- unique(names(files))
files <- lapply(dates, function(x) {
  sub <- files[[names(files) %in% x]]
  sub <- mean(sub)
  sub <- aggregate(sub, fun = mean, fact = 10)
  names(sub) <- substr(x, start = 1, stop = 10)
  return(sub)
}) %>% do.call(c, .)
writeRaster(files, "03_Data/02_CleanData/CloudCover.tif", overwrite = T)

# ################################################################################
# #### Cloud Cover 2
# ################################################################################
# # Identify the sentinel tiles that intersect with our study area
# extent <- as.polygons(ext)
# crs(extent) <- "+init=epsg:4326"
# tiles <- tiles_intersects(st_as_sf(extent, crs = 4326), out_format = "sf")
# tiles <- as(tiles, "SpatVector")
# plot(tiles)
#
# # Span a vector of dates
# dat <- seq(dates[1], dates[2], by = "day") %>%
#   as_tibble() %>%
#   setNames("Date") %>%
#   mutate(Year = year(Date), Month = month(Date)) %>%
#   group_by(Year, Month) %>%
#   summarize(From = range(Date)[1], To = range(Date)[2], .groups = "drop")
#
# # Go through the rows and check for availability of the fils
# cat("Checking available Sentinel data for the area and time of interest...\n")
# pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
# dat$Files <- lapply(1:nrow(dat), function(x) {
#
#   # List available files
#   suppressMessages(
#     files <- s2_list(
#         # spatial_extent = st_as_sf(dat$Window[[x]])
#         time_interval  = c(dat$From[x], dat$To[x])
#       , tile           = tiles$tile_id
#       , level          = "auto"
#     )
#   )
#
#   # Return the available files
#   setTxtProgressBar(pb, x)
#   return(files)
#
# })
#
# # Get Metadata on each file
# dat$FilesMetaData <- lapply(dat$Files, as.data.frame)
#
# # Bind it
# info <- do.call(rbind, dat$FilesMetaData)
#
# # Keep only unique entries
# info <- distinct(info)
#
# # Remove undesired columns
# info <- select(info, name, id_tile, sensing_datetime, clouds, footprint)
#
# # Store it
# write_csv(info, "03_Data/02_CleanData/CloudCover.csv")
# info <- read_csv("03_Data/02_CleanData/CloudCover.csv")
