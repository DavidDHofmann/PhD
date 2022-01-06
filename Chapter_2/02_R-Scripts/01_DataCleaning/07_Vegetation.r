################################################################################
#### Vegetation Data
################################################################################
# Description: use Google Earth Engine to download and process MODIS continuous
# vegetation dataset as well as NDVI data

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
library(RGISTools)    # To download satellite data

# Specify username and password to Earth Data account
username <- "DoDx9"
password <- "EarthData99"

# Load reference raster
r <- rast("03_Data/02_CleanData/00_General_Raster.tif")

################################################################################
#### MODIS Continuous Vegetation Data
################################################################################
# Specify the area of interest
aoi <- extent(c(20, 28, -22, -16))

# Check out the available products
cat("Searching for vegetation data...\n")
search <- modSearch(
    product    = "MOD44B"
  , collection = 006
  , startDate  = ymd("2011-01-01")
  , endDate    = ymd("2021-12-31")
  , extent     = extent(aoi)
)

# Function to extract information from modis filenames
modInfo <- function(urlname){
  name <- basename(urlname)
  info <- strsplit(name, split = "\\.")
  info <- lapply(info, rbind)
  info <- do.call(rbind, info)
  info <- as.data.frame(info, stringsAsFactors = F)
  names(info) <- c("Product", "AcquisitionDate", "Tile", "Version", "ProductionDate", "Format")
  info$AcquisitionDate <- substr(info$AcquisitionDate, start = 2, stop = nchar(info$AcquisitionDate))
  info$AcquisitionDate <- as.Date(info$AcquisitionDate, "%Y%j")
  info$ProductionDate <- as.Date(as.POSIXct(info$ProductionDate, format = "%Y%j%H%M%S"))
  return(info)
}

# Apply the function to the urls and create a dataframe
files <- modInfo(search$hdf)
files <- cbind(URL = search$hdf, files, stringsAsFactors = F)
head(files)

# Prepare temporary filenames for downloaded files
files$Filename <- paste0(tempdir(), "/", basename(files$URL))

# Group the files by their dates
files <- files %>% nest(Group = -AcquisitionDate)

# Let's also prepare filenames for the grouped files (we store the vegetation
# maps by date ultimately)
dir.create("03_Data/02_CleanData/00_Vegmaps", showWarnings = F)
files$FilenameGroup <- paste0(
    "03_Data/02_CleanData/00_Vegmaps/"
  , files$AcquisitionDate
  , ".grd"
)

# Subset to files that have not been downloaded yet
files <- subset(files, !file.exists(FilenameGroup))

# Check it
print(files)

# Loop through the groups and download the respective data
for (i in 1:nrow(files)) {

  # Subset to the specific group
  cat(paste0("Downloading data for ", files$AcquisitionDate[i], "\n"))
  todownload <- files$Group[[i]]

  # Download the tiles of that specific group
  for (j in 1:nrow(todownload)) {
    file <- httr::GET(
        todownload$URL[j]
      , httr::authenticate(username, password)
      , httr::progress()
      , httr::write_disk(todownload$Filename[j], overwrite = T)
    )
    cat("Tile", j, "out of", nrow(todownload), "downloaded\n")
  }

  # Load the downloaded rasters and extract the layers of interest
  downloaded <- lapply(todownload$Filename, function(x) {
    maps <- rast(x)[[1:2]]
    names(maps) <- paste0(c("TreeCover_", "ShrubCover_"), files$AcquisitionDate[i])
    return(maps)
  })

  # Separate shrubs from trees
  trees <- lapply(downloaded, function(x){
    x[[1]]
  })
  shrub <- lapply(downloaded, function(x){
    x[[2]]
  })

  # Moasic them
  cat("Mosaicing vegetation layers...\n")
  trees <- mosaic(trees[[1]], trees[[2]], trees[[3]], trees[[4]])
  shrub <- mosaic(shrub[[1]], shrub[[2]], shrub[[3]], shrub[[4]])

  # Reproject to study area
  cat("Reprojecting vegetation layers...\n")
  trees <- project(trees, r, method = "near")
  shrub <- project(shrub, r, method = "near")

  # Put layers together again
  veg <- c(trees, shrub)

  # Set water to 0
  veg <- subst(veg, 200, 0)

  # Rescale to a maximum of 1
  veg <- veg / 100

  # Write to file
  writeRaster(veg, files$FilenameGroup[[i]], overwrite = T)

}

# ################################################################################
# #### NDVI Data
# ################################################################################
# # # Specify python to use (only necessary on some systems)
# # Sys.setenv(RETICULATE_PYTHON = "/home/david/miniconda3/envs/rgee/bin/python")
#
# # Specify correct python environment
# ee_install_set_pyenv(
#     py_path = "/home/david/miniconda3/envs/rgee/bin/python"
#   , py_env  = "rgee"
# )
#
# # Make sure we have all installed for rgee
# ee_check()
# ee_Initialize()
# # ee_clean_pyenv()
# # ee_clean_credentials()
#
# # Specify an area of interest
# aoi <- ee$Geometry$Polygon(
#   list(
#       c(20, -22)
#     , c(20, -16)
#     , c(28, -16)
#     , c(28, -22)
#   )
# )
#
# ################################################################################
# #### MODIS Continuous Vegetation Dataset (MOD44B)
# ################################################################################
# # Identify Years for which we need to download vegetation data
# years <- "03_Data/02_CleanData/00_General_Dispersers.csv" %>%
#   read_csv() %>%
#   pull(Timestamp) %>%
#   range() %>%
#   as.Date() %>%
#   year()
#
# # Create a tibble that allows us to keep track of downloaded files (note that we
# # can't get a date for the last year)
# todownload <- years[1]:(years[2] - 1) %>%
#   as_tibble() %>%
#   setNames("Year") %>%
#   mutate(Filename = paste0("03_Data/02_CleanData/00_Vegmaps/", Year, ".grd"))
#
# # We only need to download files that do not exist yet
# todownload <- subset(todownload, !file.exists(Filename))
#
# # Create directory into which the downloaded maps go
# dir.create("03_Data/02_CleanData/00_Vegmaps", showWarnings = F)
#
# # Loop over the dataframe and download data
# if (nrow(todownload) > 0) {
#   cat("Downloading Vegetation data...\n")
#   for (i in 1:nrow(todownload)) {
#
#     # Get dates
#     date_from <- todownload$Year[i] %>% paste0(., "01-01") %>% ymd() %>% as.character()
#     date_to <- todownload$Year[i] %>% paste0(., "12-31") %>% ymd() %>% as.character()
#
#     # Query vegetation data for our study area during the selected years
#     query_veg <- ee$
#       ImageCollection("MODIS/006/MOD44B")$
#       filterDate(date_from, date_to)$
#       filterBounds(aoi)$
#       select("Percent_Tree_Cover", "Percent_NonTree_Vegetation")$
#       toBands()$
#       clip(aoi)
#
#     # Check it out
#     # ee_print(query_veg)
#
#     # Download the data
#     tempfile <- tempfile(fileext = ".tif")
#     ee_as_raster(query_veg, dsn = tempfile, region = aoi)
#
#     # Also store the band names
#     names_veg <- ee_print(query_veg)$img_bands_names %>%
#       strsplit(., split = " ") %>%
#       unlist() %>%
#       data.frame(BandName = ., stringsAsFactors = F) %>%
#       separate(BandName, into = c("Year", "Month", "Day", "Remove1", "Class", "Remove2"), sep = "_") %>%
#       mutate(Date = make_date(Year, Month, Day)) %>%
#       dplyr::select(Date, Class) %>%
#       mutate(Layername = paste0(Date, "_", Class)) %>%
#       pull(Layername)
#
#     # Load the downloaded data and adjust layernames
#     veg <- rast(tempfile)
#     names(veg) <- names_veg
#     #
#     # # Reproject files to the reference raster
#     # veg <- project(veg, crs(r), method = "near")
#     #
#     # # Crop them to the study area
#     # veg <- crop(veg, r, snap = "out")
#
#     # Store the file again
#     writeRaster(veg, todownload$Filename[i])
#   }
# }
#
# ################################################################################
# #### Merge Maps into "TreeCover" and "ShrubCover"
# ################################################################################
# # Once all the files are downloaded, we can merge them into a single rasterstack
# if (!file.exists("03_Data/02_CleanData/00_Vegmaps/TreeCover.grd") | !file.exists("03_Data/02_CleanData/00_Vegmaps/ShrubCover.grd")) {
#
#   # Identify all files to merge
#   veg <- "03_Data/02_CleanData/00_Vegmaps" %>%
#     dir(pattern = "[0-9].grd$", full.names = T) %>%
#     tibble(Filepath = .) %>%
#     mutate(Date = as.numeric(substr(basename(Filepath), start = 1, stop = 4))) %>%
#     mutate(Raster = map(Filepath, rast)) %>%
#     mutate(Layernames = map(Raster, names))
#
#   # Look at the dataframe
#   # print(veg)
#
#   # Let's put tree cover and nontree vegetation layers together
#   cat("Storing layers for shrub cover and tree cover...\n")
#   trees <- rast(lapply(veg$Raster, function(x) {x[[1]]}))
#   shrub <- rast(lapply(veg$Raster, function(x) {x[[2]]}))
#
#   # Put dates as layernames
#   names(trees) <- veg$Date
#   names(shrub) <- veg$Date
#
#   # Store the rasters
#   writeRaster(trees, "03_Data/02_CleanData/00_Vegmaps/TreeCover.grd", overwrite = T)
#   writeRaster(shrub, "03_Data/02_CleanData/00_Vegmaps/ShrubCover.grd", overwrite = T)
# }
#
# # To save some space, remove the old maps
# files <- dir(path = "03_Data/02_CleanData/00_Vegmaps", full.names = T)
# files <- files[!grepl(".*TreeCover.*|.*ShrubCover.*", basename(files))]
# file.remove(files)
#
# ################################################################################
# #### MODIS NDVI (MOD13Q1)
# ################################################################################
# # Identify dates for which we want the data
# dates <- "03_Data/02_CleanData/00_General_Dispersers.csv" %>%
#   read_csv() %>%
#   pull(Timestamp) %>%
#   range() %>%
#   as.Date() %>%
#   "+"(c(days(-14), days(+14)))
#
# # Create a tibble that allows us to keep track of downloaded files
# todownload <- seq(dates[1], dates[2], by = "day") %>%
#   as_tibble() %>%
#   setNames("Date") %>%
#   mutate(Year = year(Date), Month = month(Date)) %>%
#   group_by(Year, Month) %>%
#   summarize(FirstDate = range(Date)[1], LastDate = range(Date)[2], .groups = "drop") %>%
#   mutate(Filename = paste0("03_Data/02_CleanData/00_NDVI/", Year, "_", sprintf("%02d", Month), ".grd"))
#
# # We only need to download files that do not exist yet
# todownload <- subset(todownload, !file.exists(Filename))
#
# # Create directory into which the downloaded maps go
# dir.create("03_Data/02_CleanData/00_NDVI", showWarnings = F)
#
# # Loop over the dataframe and download data
# if (nrow(todownload) > 0) {
#   cat("Downloading NDVI data...\n")
#   for (i in 1:nrow(todownload)) {
#
#     # Get dates
#     date_from <- todownload$FirstDate[i] %>% as.character()
#     date_to   <- (todownload$LastDate[i] + days(1)) %>% as.character()
#
#     # Query NDVI data for our study area during dispersal dates
#     query_ndvi <- ee$
#       ImageCollection("MODIS/006/MOD13Q1")$
#       filterDate(date_from, date_to)$
#       filterBounds(aoi)$
#       select("NDVI")$
#       toBands()$
#       clip(aoi)
#
#     # Check it out
#     # ee_print(query_ndvi)
#
#     # Download the data
#     tempfile <- tempfile(fileext = ".tif")
#     ee_as_raster(query_ndvi, dsn = tempfile, region = aoi)
#
#     # Also store the band names
#     names_ndvi <- ee_print(query_ndvi)$img_bands_names %>%
#       strsplit(., split = " ") %>%
#       unlist() %>%
#       data.frame(BandName = ., stringsAsFactors = F) %>%
#       separate(BandName, into = c("Year", "Month", "Day", "Remove"), sep = "_") %>%
#       mutate(Date = make_date(Year, Month, Day)) %>%
#       pull(Date)
#
#     # Load the downloaded data and adjust layernames
#     ndvi <- rast(tempfile)
#     names(ndvi) <- names_ndvi
#
#     # Reproject files to the reference raster (I'll use the nearest neighbor
#     # method to avoid NA propagation)
#     ndvi <- project(ndvi, crs(r), method = "near")
#
#     # Crop them to the study area
#     ndvi <- crop(ndvi, r, snap = "out")
#
#     # Rescale values
#     ndvi <- ndvi * 0.0001
#
#     # Store the file again
#     writeRaster(ndvi, todownload$Filename[i])
#
#   }
# }
#
# # Once all the files are downloaded, we can merge them into a single rasterstack
# # Let's first make sure that there are no duplicates
# if (!file.exists("03_Data/02_CleanData/00_NDVI/NDVI.grd")) {
#
#   # Get all files to merge
#   ndvi <- "03_Data/02_CleanData/00_NDVI" %>%
#     dir(pattern = ".grd$", full.names = T) %>%
#     tibble(Filepath = .) %>%
#     mutate(Raster = map(Filepath, rast)) %>%
#     mutate(Layernames = map(Raster, names))
#
#   # Look at the dataframe
#   print(ndvi)
#
#   # Check for duplicates
#   unlist(ndvi$Layernames) %>%
#     duplicated() %>%
#     table()
#
#   # Put all files together and store the raster
#   ndvi <- rast(ndvi$Raster)
#
#   # Store file
#   writeRaster(ndvi, "03_Data/02_CleanData/01_LandCover_NDVI.grd", overwrite = T)
# }
#
# # To save some space, remove the old maps
# files <- dir(path = "03_Data/02_CleanData/00_NDVI", full.names = T)
# files <- files[!grepl(".*NDVI.*", basename(files))]
# print(files)
# file.remove(files)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/07_Vegetation.rds")
cat("Done :)\n")
