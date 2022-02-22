################################################################################
#### Sentinel Prerequisits for Data Download
################################################################################
# Description: Automatic download of sentinel 2 data. This data will be used to
# map the location of wet pans. Note that we use a moving window to download
# data only for the region of interest at a given date. The region of interest
# will be driven by the location of GPS data of dispersing wild dogs at that
# specific point in time.

# Clear R's brain
rm(list = ls())

# Load required packages
library(sen2r)         # To automate the download of Sentinel 2 data
library(raster)        # To handle spatial data
library(sf)            # To handle spatial data
library(terra)         # To handle spatial data
library(rgeos)         # To handle spatial data
library(tidyverse)     # For data wrangling
library(lubridate)     # To handle dates
library(pbmcapply)     # For multicore abilities with progress bar

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Login to scihub
write_scihub_login("dodx9", "Scihubbuster69_")

################################################################################
#### Check Overlap with Sentinel Tiles
################################################################################
# Get a rough overview of the different sentinel tiles and how many overlap with
# our main study area. Hence, let's load a shapefile of study area
study <- shapefile("03_Data/02_CleanData/00_General_Shapefile.shp")

# Check with which tiles of the sentinel grid our study area overlaps
tiles <- as(tiles_intersects(st_as_sf(study), out_format = "sf"), "Spatial")

# Also load dispersal data
disp <- "03_Data/02_CleanData/00_General_Dispersers.csv" %>%
  read_csv() %>%
  subset(State == "Disperser")

# Make it spatial
coordinates(disp) <- c("x", "y")
crs(disp) <- "+init=epsg:4326"

# Visualize
plot(tiles, main = paste0("Number of Tiles: ", length(tiles)))
plot(study, add = T, border = "red", col = adjustcolor("red", alpha.f = 0.2))
text(tiles, "tile_id", cex = 0.5)
plot(disp, add = T, pch = 20, cex = 0.5, col = adjustcolor("blue", alpha.f = 0.1))
axis(1)
axis(2)

################################################################################
#### Generate Spatio-Temporal Moving Window
################################################################################
# Instead of downloading satellite data for the entire study area and the entire
# range of dates, I'm going to create a spatio-temporal moving window over the
# locations of our dispersing wild dogs. This should help to reduce the amount
# of data that needs to be downloaded. Thus, let's load the dispersal data again
disp <- "03_Data/02_CleanData/00_General_Dispersers.csv" %>%
  read_csv() %>%
  subset(State == "Disperser")

# Let's group the GPS data by month, so that we can download and prepare
# satellite data for each month
disp <- disp %>%
  mutate(Year = year(Timestamp), Month = month(Timestamp)) %>%
  group_by(Year, Month) %>%
  nest() %>%
  setNames(c("Year", "Month", "GPS")) %>%
  arrange(Year, Month)

# Create a bounding box around the coordinates of each week. Note that I'm going
# to buffer each of those boxes substantially, as the buffer needs to be large
# enough such that potential random steps from the SSF analysis are covered, but
# also that the "distance to" maeasurements can be done meaningfully
cat("Create spatio-temporal moving windows...\n")
disp$Window <- pbmclapply(
    X                  = disp$GPS
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {

  # Make coordinates spatial
  pts <- vect(cbind(x$x, x$y))
  crs(pts) <- "+init=epsg:4326"

  # Buffer the points
  buff <- buffer(pts, width = 100 * 1000)
  buff <- aggregate(buff)
  buff <- disagg(buff)

  # Create extents around buffered points
  ext <- lapply(buff, function(y) {
    ex <- as.polygons(ext(y))
    crs(ex) <- "+init=epsg:4326"
    return(ex)
  })

  # Return result as SpatialPolygonsDataFrame
  if (length(ext) > 1) {
      ext <- do.call(rbind, ext)
    } else {
      ext <- ext[[1]]
  }
  ext <- as(ext, "Spatial")
  ext <- as(ext, "SpatialPolygonsDataFrame")
  return(ext)
})

# Let's check again the overlap of each moving window with the tiles of Sentinel
wind <- do.call(rbind, disp$Window)
ext <- as(extent(wind), "SpatialPolygons")
crs(ext) <- "+init=epsg:4326"
tiles <- tiles_intersects(st_as_sf(ext), out_format = "sf")
tiles <- as(tiles, "Spatial")

# Visualize
plot(tiles, lwd = 0.2, main = "Sentinel Tiles")
text(tiles, "tile_id", cex = 0.6)
plot(wind, add = T, border = "blue")

# Some windows only overlap very minimally with the sentinel tiles. Hence, let's
# identify the overlap of the moving windows with each tile
cat("Computing overlap of tiles with moving windows...\n")
disp$Tiles <- pbmclapply(
    X                  = 1:nrow(disp)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {

  # Get the window(s) of the respective date
  win <- disp$Window[[x]]
  win <- vect(win)
  til <- vect(tiles)

  # Identify with which tiles the windows overlap/intersect
  index <- colSums(relate(win, til, "intersects")) > 0
  index <- which(index)
  result <- til[index, ]

  # Compute the percentage of overlap of each tile with the moving window(s)
  percent_overlap <- sapply(index, function(y) {
    total <- expanse(til[y, ])
    overl <- expanse(crop(til[y], win))
    perc <- round(overl / total, 5)
    return(perc)
  })

  # Put everything into a dataframe and return the results
  values(result) <- data.frame(TileID = tiles$tile_id[index], Overlap = percent_overlap)
  result <- as(result, "Spatial")
  return(result)
})

# Let's remove tiles that are not overlapping more than 5 percent
disp$Tiles <- lapply(disp$Tiles, function(x) {
  subset(x, Overlap > 0.05)
})

# Visualize them
i <- sample(nrow(disp), size = 1)
plot(disp$Tiles[[i]], main = paste0(disp$Year[i], " - ", disp$Month[i]))
plot(disp$Window[[i]], add = T, border = "red", lwd = 2)
points(y ~ x, data = disp$GPS[[i]], col = "red", pch = 20)
text(disp$Tiles[[i]], "Overlap", cex = 0.8)
axis(1)
axis(2)

# Store the object to file
write_rds(disp, "/media/david/Elements/Windows.rds")

################################################################################
#### Prepare Monthly Download Dates
################################################################################
# We want to download by month, so let's generate sequences for all dates that
# we want to download. First, create all possible dates for the duration of our
# study and nest the dates by year and month.
sent <- seq(ymd("2011-01-01"), ymd("2021-12-31"), by = "day") %>%
  tibble(Year = year(.), Month = month(.), Date = .) %>%
  group_by(Year, Month) %>%
  nest()

# Obviously we don't need to download data for all of those months, but only for
# those that are present in the dispersal data. Hence, let's only keep the
# combinations of "Year" and "Week" that are also in the dispersal data
# index <- paste0(sent$Year, "_", sent$Week) %in% paste0(disp$Year, "_", disp$Week)
index <- paste0(sent$Year, "_", sent$Month) %in%
  paste0(disp$Year, "_", disp$Month)

# Subset to those dates and determine start-end dates for each month
sent <- subset(sent, index) %>%
  mutate(
      From = map(data, function(x) {pull(x, Date) %>% min()}) %>% do.call(c, .)
    , To   = map(data, function(x) {pull(x, Date) %>% max()}) %>% do.call(c, .)
  )

# Take a look at the download schedule
print(sent)

# Let's merge the "disp" and "sent" tables. Note that they have the same amount
# of rows
nrow(disp) == nrow(sent)
dat <- full_join(disp, sent, by = c("Year", "Month"))

# Remove anything that is not needed anymore to avoid confusion
rm(disp, sent, ext, study, tiles, wind, index)

# Sentinel data is only available since 26.03.2015, so let's subset to data
# since 2015
dat <- subset(dat, Year >= 2015)

# Go through the rows and check for availability of the fils
cat("Checking available Sentinel data for the area and time of interest...\n")
pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
dat$Files <- lapply(1:nrow(dat), function(x) {

  # List available files
  suppressMessages(
    files <- s2_list(
        spatial_extent = st_as_sf(dat$Window[[x]])
      , time_interval  = c(dat$From[x], dat$To[x])
      , level          = "auto"
    )
  )

  # Return the available files
  setTxtProgressBar(pb, x)
  return(files)

})

# Get Metadata on each file
dat$FilesMetaData <- lapply(dat$Files, as.data.frame)

# Store the results to file
write_rds(dat, "03_Data/03_Results/99_SentinelResults.rds")
dat <- read_rds("03_Data/03_Results/99_SentinelResults.rds")

# Check how many available files there are
dat$NumberFiles <- sapply(dat$Files, length)
sum(dat$NumberFiles)

# Check how many 1A and 1C products there are
dat %>%
  select(Year, Month, FilesMetaData) %>%
  unnest(FilesMetaData) %>%
  ungroup() %>%
  count(level)

# Let's create a clean tibble of all files that we need to download
metadata <- do.call(rbind, dat$FilesMetaData)
todownload <- do.call(c, dat$Files)
todownload <- tibble(todownload, metadata)
print(todownload)

# Store all to file
write_rds(todownload, "/media/david/Elements/Todownload.rds")
