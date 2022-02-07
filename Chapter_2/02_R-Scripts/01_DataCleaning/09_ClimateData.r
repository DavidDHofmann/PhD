################################################################################
#### Download and Curate Climate Data
################################################################################
# Description: Use Google Earth Engine to download hourly precipitation and
# temperature data for our study area

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(reticulate)   # Interface to python
library(rgee)         # Interface to google earth engine
library(terra)        # To handle spatial data
library(tidyverse)    # To wrangle data
library(lubridate)    # To handle dates

# # Specify python to use (only necessary on some systems)
# Sys.setenv(RETICULATE_PYTHON = "/home/david/miniconda3/envs/rgee/bin/python")

# Specify correct python environment
ee_install_set_pyenv(
    py_path = "/home/david/miniconda3/envs/rgee/bin/python"
  , py_env  = "rgee"
)

# Make sure we have all installed for rgee
ee_check()
ee_Initialize()
# ee_clean_pyenv()
# ee_clean_credentials()

# Specify an area of interest
aoi <- ee$Geometry$Polygon(
  list(
      c(20, -22)
    , c(20, -16)
    , c(28, -16)
    , c(28, -22)
  )
)

# Load the reference raster
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")

################################################################################
#### Download Precipitation Data
################################################################################
# Instead of trying to download data for the entire range of dates at once, I
# want to download the maps by month. This will allow me to observe the download
# process better and to save files in between. Let's first identify the dates
# for which we want to download data.
dates <- "03_Data/02_CleanData/00_General_Dispersers.csv" %>%
  read_csv() %>%
  pull(Timestamp) %>%
  range() %>%
  as.Date() %>%
  "+"(c(days(-7), days(+7)))

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
    precip <- crop(precip, r, snap = "out")

    # Store the file again
    writeRaster(precip, todownload$Filename[i])
  }
}

################################################################################
#### Temperature Data
################################################################################
# Instead of trying to download data for the entire range of dates at once, I
# want to download the maps by month. This will allow me to observe the download
# process better and to save files in between. Let's first identify the dates
# for which we want to download data.
dates <- "03_Data/02_CleanData/00_General_Dispersers.csv" %>%
  read_csv() %>%
  pull(Timestamp) %>%
  range() %>%
  as.Date() %>%
  "+"(c(days(-7), days(+7)))

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

    # Query rainfall data for our study area during dispersal dates
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
    temp <- crop(temp, r, snap = "out")

    # Convert kelvin to celsius
    temp <- temp - 273.15

    # Store the file again
    writeRaster(temp, todownload$Filename[i])
  }
}

################################################################################
#### Aggregate Precipitation Data
################################################################################
# Prepare a tibble for easier subsetting of the data
cat("Loading precipitation data and generating aggregation table...\n")
precip <- "03_Data/02_CleanData/00_Rainmaps" %>%
  dir(pattern = ".grd$", full.names = T) %>%
  rast() %>%
  as.list() %>%
  tibble(Raster = .)

# Extract date and hour from layernames
precip <- mutate(precip
  , Name = sapply(Raster, names)
  , Date = substr(Name, start = 1, stop = 10)
  , Hour = substr(Name, start = 12, stop = 13) %>% as.numeric()
  , Timestamp = make_datetime(year(Date), month(Date), day(Date), Hour, 0, 0)
)

# Reorder columns
precip <- dplyr::select(precip, Timestamp, Date, Hour, Name, Raster)

# We don't need four hourly data. Hence, we'll aggregate the temperature data to
# four hours (except for 11 o clock for which we'll use 8 hours instead),
# corresponding to the GPS fixes. Let's create a help table for this.
help <- seq(min(precip$Timestamp) + hours(8), max(precip$Timestamp) - hours(8), by = "hour")
help <- help[hour(help) %in% c(3, 7, 15, 19, 23)]
help <- tibble(Timestamp = help)

# For each timestamp in the help table we can now identify the temperature maps
# that we need to average
averaged <- mutate(help, Average = map(Timestamp, function(x) {
  if (hour(x) != 15) {
    sub <- subset(precip, Timestamp <= x & Timestamp >= (x - hours(4)))
  } else {
    sub <- subset(precip, Timestamp <= x & Timestamp >= (x - hours(8)))
  }
  return(sub$Raster)
}))

# Check it
print(averaged)

# Each entry should either contain 5 or 9 maps. Let's check this
table(sapply(averaged$Average, length))

# For some reason we are missing the 00:00 data for one entry. I tried multiple
# time to access this data but it's not possible. Anyways, it shouldn't
# influence our results at all.

# Let's average the respective rasters
cat("Aggregating precipitation data...\n")
n <- nrow(averaged)
pb <- txtProgressBar(min = 0, max = n, style = 3)
averaged$Averaged <- lapply(1:nrow(averaged), function(x) {
  average <- rast(averaged$Average[[x]])
  average <- mean(average)
  names(average) <- averaged$Timestamp[[x]]
  setTxtProgressBar(pb, x)
  return(average)
})

# Create a single stack of rasters
averaged <- rast(averaged$Averaged)
names(averaged)

# Visualize some of the maps
plot(averaged[[sample(nlyr(averaged), 4)]])

# Store the final maps to file
writeRaster(averaged, "03_Data/02_CleanData/05_Climate_Precipitation.grd")

################################################################################
#### Aggregate Temperature Data
################################################################################
# Prepare a tibble for easier subsetting of the data
cat("Loading temperature data and generating aggregation table...\n")
temp <- "03_Data/02_CleanData/00_Tempmaps" %>%
  dir(pattern = ".grd$", full.names = T) %>%
  rast() %>%
  as.list() %>%
  tibble(Raster = .)

# Extract date and hour from layernames
temp <- mutate(temp
  , Name = sapply(Raster, names)
  , Date = substr(Name, start = 1, stop = 10)
  , Hour = substr(Name, start = 12, stop = 13) %>% as.numeric()
  , Timestamp = make_datetime(year(Date), month(Date), day(Date), Hour, 0, 0)
)

# Reorder columns
temp <- dplyr::select(temp, Timestamp, Date, Hour, Name, Raster)

# We don't need four hourly data. Hence, we'll aggregate the temperature data to
# four hours (except for 11 o clock for which we'll use 8 hours instead),
# corresponding to the GPS fixes. Let's create a help table for this.
help <- seq(min(temp$Timestamp) + hours(8), max(temp$Timestamp), by = "hour")
help <- help[hour(help) %in% c(3, 7, 15, 19, 23)]
help <- tibble(Timestamp = help)

# For each timestamp in the help table we can now identify the temperature maps
# that we need to average
averaged <- mutate(help, Average = map(Timestamp, function(x) {
  if (hour(x) != 15) {
    sub <- subset(temp, Timestamp <= x & Timestamp >= (x - hours(4)))
  } else {
    sub <- subset(temp, Timestamp <= x & Timestamp >= (x - hours(8)))
  }
  return(sub$Raster)
}))

# Check it
print(averaged)

# Each entry should either contain 5 or 9 maps. Let's check this
check <- sapply(averaged$Average, length)
table(check)

# Again there appear to be a couple of days were data is missing (all in
# September 2021). We'll have to skip those dates.
averaged$Timestamp[check < 5]

# Subset accordingly
averaged <- subset(averaged, check > 0)

# Let's average the respective rasters
cat("Aggregating temperature data...\n")
n <- nrow(averaged)
pb <- txtProgressBar(min = 0, max = n, style = 3)
averaged$Averaged <- lapply(1:nrow(averaged), function(x) {
  average <- rast(averaged$Average[[x]])
  average <- mean(average)
  names(average) <- averaged$Timestamp[[x]]
  setTxtProgressBar(pb, x)
  return(average)
})

# Create a single stack of rasters
averaged <- rast(averaged$Averaged)
names(averaged)

# Visualize some of the maps
plot(averaged[[sample(nlyr(averaged), 4)]])

# Store the final maps to file
writeRaster(averaged, "03_Data/02_CleanData/05_Climate_Temperature.grd")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/09_ClimateData.rds")
cat("Done :)\n")
