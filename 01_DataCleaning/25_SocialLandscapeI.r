############################################################
#### Preparing the Social Landscape
############################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(data.table)
library(viridis)
library(adehabitatHR)
library(raster)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(rasterVis)
library(rgdal)
library(rgeos)
library(parallel)

# Load custom functions
source("Functions.r")

# Set the temporary directory (we need to do this because the job takes so long
# that it eventually cleans the default temporary r folder)
rasterOptions(tmpdir = "/home/david/Schreibtisch/Temp")

# Make use of multicore abilities
beginCluster()

############################################################
#### Preprocess GPS Fixes
############################################################
# Read in all gps fixes and do some cleaning
dat <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv() %>%
  subset(!is.na(CurrentPack)) %>%
  subset(State != "Disperser") %>%

  # For some packs we might have more than one gps fix taken at the same time.
  # Let's get rid of those
  distinct(CurrentPack, Timestamp, .keep_all = TRUE) %>%
  mutate(Date = as.Date(Timestamp)) %>%
  group_by(CurrentPack, Date) %>%
  nest() %>%

  # Resample fixes to one fix a day. Keep the fix that is closest to 07:00:00
  mutate(Resampled = map(data, function(x){
    refTimestamp <- update(x$Timestamp, hour = 7, min = 0, sec = 0)
    keep <- subset(x
      , Timestamp == Timestamp[which.min(abs(Timestamp - refTimestamp))]
    )
    return(keep)
  })) %>%
  unnest(Resampled) %>%
  dplyr::select(-c("Date", "data"))

# I will use this opportunity to create an extent object that shows the
# boundaries within which all of our GPS data has been collected
ext <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv()
ext <- extent(c(
    min(ext$x)
  , max(ext$x)
  , min(ext$y)
  , max(ext$y)
))
ext <- as(ext, "SpatialPolygons")
crs(ext) <- CRS("+init=epsg:4326")

# Prepare list of dates for which we want to prepare social landscapes. Let's
# check the dispersal dates for this
dates <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv() %>%
  subset(., State == "Disperser") %>%
  group_by(., DogName) %>%
  summarise(StartDate = range(Timestamp)[1], EndDate = range(Timestamp)[2])

# The easiest solution would know be to simply use the earliest and latest dates
# and then create a social landscape for each month in between the two.
easy <- seq(
    from  = min(dates$StartDate)
  , to    = max(dates$EndDate)
  , by    = "month"
)

# However, we can see that this would require us to calculate >90 social
# landscapes, of which we would actually only need a few
length(easy)

# Since only Abrahms' individuals dispersed before 2016 we can define the dates
# required for them manually
stetson   <- c("2012-01-01 00:00:00", "2012-02-01 00:00:00")
scorpion  <- c("2013-08-01 00:00:00", "2013-09-01 00:00:00")
lupe      <- c("2014-10-01 00:00:00", "2014-11-01 00:00:00")

# Put them together and convert them to posix
dates1 <- c(stetson, scorpion, lupe)
dates1 <- as.POSIXct(dates1, tz = "UTC")

# For our own dispersers we create social landscapes for each month since
# beginning of 2016 until know
first <- as.POSIXct("2016-01-01 00:00:00", tz = "UTC")
last  <- max(dates$EndDate)
dates2 <- seq(first, last, by = "month")

# Combine the dates
range <- c(dates1, dates2)

# Due to the fact that we are missing sufficient data for the periods prior to
# 2014, we are only going to calculate social landscape for dates after 2013
range <- subset(range, year(range) > 2013)

# Let's check again the number of social landscapes that we will compute
length(range)

# For each date we now prepare the corresponding data that is used for the
# social landscapes. In other words, for each date we subset to the data 6
# months prior to the date itself
sets <- list()
for (i in 1:length(range)){
  sets[[i]] <- subset(dat,
    Timestamp <= range[i] &
    Timestamp >= range[i] - months(6)
  )
}
names(sets) <- range

# Collapse the list while preserving the group name as a new column
dat <- rbindlist(sets, idcol = TRUE) %>%
  as.data.frame() %>%
  rename(Group = .id) %>%
  group_by(Group, CurrentPack) %>%
  nest()

# Look at the final tibble
print(dat, n = 50)

# Note that we can only estimate uds for packs with at least 5 fixes. Let's
# check for each pack the number of fixes at a given date and keep only those
# with at least 5 fixes
dat <- dat %>%
  mutate(NoFixes = map(data, nrow) %>%
    do.call(rbind, .) %>%
    as.vector()
  ) %>%
  subset(NoFixes >= 5) %>%
  dplyr::select(-NoFixes)

# Lastly we want to make our gps fixes spatial
dat <- dat %>% mutate(data = map(data, function(x){
    SpatialPointsDataFrame(
        x[, c("x", "y")]
      , data        = x
      , proj4string = CRS("+init=epsg:4326"))
  }))

# Look at the final table
print(dat, n = 50)

############################################################
#### Preprocess Reference Raster
############################################################
# Load reference raster to calculate utility distributions (uds) on
r <- raster("03_Data/02_CleanData/00_General_Raster250.tif")

# Check how the extent and reference raster overlap
plot(r)
plot(ext, add = T)

# Let's slightly buffer the extent
ext2 <- extent(ext) + 1/111 * c(-20, 20, -20, 20)
ext2 <- as(ext2, "SpatialPolygons")
crs(ext2) <- crs(ext)

# Plot again
plot(r)
plot(ext, add = T)
plot(ext2, add = T, border = "red")

# Crop the reference raster to this extent
r <- crop(r, ext2)

# The raster needs to be converted into a SpatialPixelsDataFrame in order to
# work with the adehabitatHR package
r <- as(r, "SpatialPixelsDataFrame")

# Reproject the gps fixes to match the crs of the reference raster
dat$data <- lapply(dat$data, function(x){spTransform(x, crs(r))})

# Now calculate all uds and hrs
dat <- dat %>%

  # Calculate UD
  mutate(ud = mclapply(data, mc.cores = (detectCores() - 1), function(x){
    ud_ref <- kernelUD(
        x
      , h     = "href"
      , grid  = r
      # , grid  = 150
      # , ext   = 4
    )

    # Calculate HR
    hr <- getverticeshr(ud_ref, percent = 95)

    # Identify number of polygons of the HR using href
    npoly <- length(hr@polygons[[1]]@Polygons)

    # Identify href
    href <- ud_ref@h$h

    # Calculate upper and lower bound for optimized href
    href_new <- (href * 0.3 + href * 1.2) / 2

    # Calculate increment
    increment <- href_new - 0.3 * href

    # Initialize index for loop (will increase every iteration)
    index <- 1

    # Run loop that optimizes href
    while (increment * 0.5 ** index > 0.001){

      # Increase or decrease href depending on number of polygons
      if (npoly > 1){
        href_new <- href_new + increment * 0.5 ** index
      } else {
        href_new <- href_new - increment * 0.5 ** index
      }

      # Update the counter
      index <- index + 1

      # Calculate the new ud for the specific pack
      ud_new <- kernelUD(
          x
        , h     = href_new
        , grid  = r
        # , grid  = 150
        # , ext   = 4
      )

      # Get the new home ranges
      hr_new <- getverticeshr(ud_new, percent = 95)

      # Update the number of polygons
      npoly <- length(hr_new@polygons[[1]]@Polygons)
    }

    # The resulting h might still be just below the h that ensures a single polygon.
    # Let's thus add 0.001
    href_new <- href_new + 0.001

    # Calculate ud again
    ud_new <- kernelUD(
        x
      , h     = href_new
      , grid  = r
      # , grid  = 150
      # , ext   = 4
    )

    # Return the optimized ud
    return(ud_new)
  })) %>%

  # Calculate the hrs from the optimized uds
  mutate(hr = mclapply(ud, mc.cores = (detectCores() - 1), function(x){
    getverticeshr(x, percent = 95)
  })) %>%

  # Coerce the uds to raster files and store them to a temporary file to prevent
  # memory overflow
  mutate(ud = mclapply(ud, mc.cores = (detectCores() - 1), function(x){
    y <- raster(x)
    y <- writeRaster(y, rasterTmpFile())
    gc()
    return(y)
  }))

# Save data to file for part II
saveRDS(dat, "99_SocialLandscape(PartI).rds")

# End the cluster
endCluster()
