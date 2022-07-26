################################################################################
#### Download and Curate Climate Data
################################################################################
# Description: Use Google Earth Engine to download hourly precipitation and
# temperature data for our study area

# Clean environment
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
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
