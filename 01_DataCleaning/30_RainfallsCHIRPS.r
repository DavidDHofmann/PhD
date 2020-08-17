############################################################
#### Download CHIRPS rainfall data
############################################################
# Description: Download of CHIRPS data for africa

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # For spatial data analysis
library(sf)             # For spatial data analysis
library(RColorBrewer)   # To get nice colors
library(viridis)        # To get more nice colors
library(lubridate)      # To handle timestamps
library(ggdark)         # For plotting a dark ggplot theme
library(tmap)           # To plot spatial data
library(pracma)         # To find peaks in timeseries data
library(astsa)          # To visualize lagged correlation

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load custom functions
source("Functions.r")

# Make use of multicore abilities
beginCluster()

############################################################
#### Download Data
############################################################
# Write a function that downloads the desired chirps data
getCHIRPS <- function(
    tres      = c("daily", "monthly")
  , sres      = c(0.05, 0.25)
  , begin     = NULL
  , end       = NULL
  , dsn       = getwd()
  , overwrite = FALSE
  ){

  # Prepare dates
  dates <- seq(as.Date(begin), as.Date(end), by = "day")

  # Determine filenames
  if (tres == "daily"){
    filename <- paste0("chirps-v2.0.", format(dates, "%Y.%m.%d"), ".tif")

    # Daily images before 2018 are zipped
    index <- which(dates <= as.Date("2016-04-30"))

    # Adjust the url
    filename[index] <- paste0(filename[index], ".gz")

  } else {
    filename <- paste0("chirps-v2.0.", format(dates, "%Y.%m"), ".tif.gz")
  }

  # Keep only unique filenames
  filename <- unique(filename)

  # Write function to determine URL
  url_prep <- function(dates, tres, sres){

    # Define the fixed part of the url
    fixed <- "ftp://ftp.chg.ucsb.edu/pub/org/chg/products/CHIRPS-2.0/"

    # Define the variable part, depending on daily or monthly data
    if (tres == "daily"){

      # For daily images
      variable <- paste0("africa_daily/tifs/p"
        , substr(sres, start = 3, stop = 4)
        , "/"
        , lubridate::year(dates)
        , "/"
        , filename
      )
      url <- paste0(fixed, "/", variable)

      return(url)
    } else {

      # For monthly images
      variable <- paste0("africa_monthly/tifs"
        , "/"
        , filename
      )
      url <- paste0(fixed, "/", variable)
      return(url)
    }
  }

  # Get all urls (keep only uniques)
  url <- url_prep(dates, tres, sres)

  # Define the destination
  destfile <- paste0(dsn, "/", filename)

  # Define the unzipped filename
  destfile_unzip <- gsub(pattern = ".gz", "", destfile)

  # Prepare an empty vector
  exists <- rep(FALSE, length(url))

  # Download files
  for (i in 1:length(url)){

    # Check if the file already exists as tif
    if (file.exists(destfile_unzip[i]) & overwrite == FALSE){

      # If it exists, we tell the user...
      print(paste0("File already exists"))
      exists[i] <- TRUE

      # ... and we skip an iteration
      next
    } else {

      # Otherwise we download the file
      download.file(url[i], destfile = destfile[i], quiet = TRUE)
    }

    # Print status
    cat(i, "of", length(url), "done...\n")
  }

  # For those files that are zipped we need to unzip them
  index <- grep(".gz", destfile)
  index <- index[index %in% which(!exists)]
  if (length(index) > 0){
    for (i in 1:length(index)){
      R.utils::gunzip(destfile[index[i]]
        , overwrite = overwrite
        , destname = destfile_unzip[i]
      )}
  }

  # Return the file paths
  return(destfile_unzip)
}

# We want data for the time during which we have modis floodmaps
dates <- "03_Data/02_CleanData/00_Floodmaps/01_Original" %>%
  dir(. , pattern = ".tif$") %>%
  as.Date(format = "%Y.%m.%d") %>%
  range(.)

# Note that we want data prior to the first MODIS image as well
dates[1] <- dates[1] - months(12)

# Run the function
files <- getCHIRPS(
    tres  = "monthly"
  , sres  = 0.25
  , begin = dates[1]
  , end   = dates[2]
  , dsn   = "03_Data/01_RawData/CHIRPS"
)

############################################################
#### Crop Rain Data to Catchment Areas
############################################################
# Load the files as a stack
rain <- dir(
    path        = "03_Data/01_RawData/CHIRPS"
  , pattern     = ".tif$"
  , full.names  = T
) %>% stack()

# Load the catchment area shapefile and do some cleaning
catch <- "03_Data/01_RawData/DAVID/Catchments.shp" %>%
  st_read() %>%
  mutate(ID = 1:nrow(.)) %>%
  dplyr::select(-c("NAME2")) %>%
  rename(Name = NAME1, Level = LEVEL)

# Plot them
tm_shape(catch) + tm_polygons(col = "ID", palette = "viridis", n = nrow(catch))
tm_shape(catch) + tm_polygons(col = "Level", palette = "viridis")

# Crop the rainmaps to the extent of the catchment shapefile
rain <- mccrop(rain, extent(catch))

# Plot some of the resulting maps
plot(rain[[1:9]], col = rev(brewer.pal(10, "Spectral")))

# Let's plot a rainmap and the catchment areas
plot(rain[[1]], col = rev(brewer.pal(10, "Spectral")))
plot(catch$geometry, border = "black", add = T)

############################################################
#### Explore Rain Data
############################################################
# Extract dates from rainmaps
dates <- rain %>%
  names() %>%
  substr(start = nchar(.) - 6, stop = nchar(.)) %>%
  paste0(., ".15") %>%
  as.Date(format = "%Y.%m.%d")

# Extract rainfalls below each catchment area
catch <- mutate(catch, Rainfall = map(ID, function(x){
  masked <- mask(rain, subset(catch, ID == x))
  totalrain <- cellStats(masked, "sum")
  totalrain <- data.frame(Date = dates, Rain = totalrain)
  rownames(totalrain) <- NULL
  return(totalrain)
}))

# Unnest the tibble to get a nice dataframe
catch <- unnest(catch, cols = c(Rainfall))

# Remove geometries (will speed up calculations a lot)
st_geometry(catch) <- NULL

# Plot the rainfalls throughout time
catch %>%
  group_by(Date) %>%
  summarize(Rain = sum(Rain)) %>%
  ggplot(aes(x = Date, y = Rain)) +
  geom_line(color = viridis(20)[20]) +
  dark_theme_grey()

# Plot the rainfalls throughout a year
catch %>%
  group_by(Date) %>%
  summarize(Rain = sum(Rain)) %>%
  mutate(Year = year(Date)) %>%
  mutate(Month = month(Date)) %>%
  ggplot(aes(x = Month, y = Rain, col = factor(Year))) +
  geom_line() +
  scale_color_viridis_d() +
  dark_theme_grey()

# Plot the rainfalls throughout some years
catch %>%
  group_by(Date) %>%
  summarize(Rain = sum(Rain)) %>%
  mutate(Year = year(Date)) %>%
  mutate(Month = month(Date)) %>%
  subset(Year %in% c(2010:2016)) %>%
  ggplot(aes(x = Date, y = Rain)) +
  geom_line(color = viridis(20)[20]) +
  scale_color_viridis_d() +
  dark_theme_grey()

# Note: There appears to be a tendency for double-peaking, with one peak around
# new year and another one around march.

# Check if there are substantial differences between catchment areas
catch %>%
  group_by(Date, ID) %>%
  summarize(Rain = sum(Rain)) %>%
  mutate(Year = year(Date)) %>%
  mutate(Month = month(Date)) %>%
  ggplot(aes(x = Month, y = Rain, col = factor(Year))) +
  geom_line() +
  scale_color_viridis_d() +
  dark_theme_grey() +
  facet_wrap("ID")

# Plot the outcome by region
catch %>%
  group_by(Date, Name) %>%
  summarize(Rain = sum(Rain)) %>%
  mutate(Year = year(Date)) %>%
  mutate(Month = month(Date)) %>%
  ggplot(aes(x = Month, y = Rain, col = factor(Year))) +
  geom_line() +
  scale_color_viridis_d() +
  dark_theme_grey() +
  facet_wrap("Name")

# Note: It doesn't really look like there are massive differences in the
# patterns of rainfall across the regions (obviously there are differences in
# the amount of rainfall though). I will therefore simplify the dataframe and
# consolidate all regions
catch <- catch %>%
  group_by(Date) %>%
  summarize(Rain = sum(Rain))

# Store the dataframe to file
write_csv(
    catch
  , "03_Data/02_CleanData/06_EnvironmentalFeatures_RainfallsCatchment_CHIRPS.csv"
)

############################################################
#### Prepare Flood Data
############################################################
# Load all floodmaps
flood <- dir(
    path        = "03_Data/02_CleanData/00_Floodmaps/01_Original"
  , pattern     = ".tif$"
  , full.names  = T) %>% stack()

# Let's check the cloud cover (127) in each image
clouds <- freq(flood, value = 127)

# We can use this to identify outliers
plot(clouds, type = "l")

# Let's check in which images the cloud cover exceeds 20 percent
clouds_prop <- clouds / ncell(flood[[1]])
sum(clouds_prop > 0.2)
index <- which(clouds_prop > 0.2)

# Let's see how badly it looks
plot(flood[[index]])

# We can remove them
index <- which(clouds_prop <= 0.2)
flood <- flood[[index]]

# Reclassify the remaining images
rcl <- data.frame(old = c(0, 127, 255), new = c(1, 0, 0))
flood <- mcreclassify(flood, rcl)

# Mask pixels outside the delta
flood <- "03_Data/01_RawData/DAVID/Catchments.shp" %>%
  st_read() %>%
  subset(NAME1 == "Delta") %>%
  mask(flood, .)

# Calculate the number of inundated pixels in all floodmaps
flood_sum <- cellStats(flood, "sum")

# Extract the dates from the file descriptions
dates <- flood %>%
  names(.) %>%
  substr(., start = 2, stop = nchar(.)) %>%
  as.Date(., format = "%Y.%m.%d")

# Create dataframe
dat_flood <- data.frame(Flood = flood_sum, Date = dates) %>%
  mutate(Year = year(Date)) %>%
  mutate(Month = month(Date)) %>%
  mutate(Day = yday(Date))

# Remove rownames
rownames(dat_flood) <- NULL

# Plot the outcome throughout all years
dat_flood %>%
  ggplot(aes(x = Date, y = Flood)) +
  geom_line(color = viridis(20)[20]) +
  dark_theme_grey()

# Plot the outcome throughout all years
dat_flood %>%
  subset(Year %in% c(2000:2005)) %>%
  ggplot(aes(x = Date, y = Flood)) +
  geom_line(color = viridis(20)[20]) +
  dark_theme_grey()

dat_flood %>%
  subset(Year %in% c(2006:2010)) %>%
  ggplot(aes(x = Date, y = Flood)) +
  geom_line(color = viridis(20)[20]) +
  dark_theme_grey()

dat_flood %>%
  subset(Year %in% c(2011:2015)) %>%
  ggplot(aes(x = Date, y = Flood)) +
  geom_line(color = viridis(20)[20]) +
  dark_theme_grey()

dat_flood %>%
  subset(Year %in% c(2016:2020)) %>%
  ggplot(aes(x = Date, y = Flood)) +
  geom_line(color = viridis(20)[20]) +
  dark_theme_grey()

# Plot the outcome throughout a year
dat_flood %>%
  ggplot(aes(x = Day, y = Flood, col = factor(Year))) +
  geom_line() +
  scale_color_viridis_d() +
  dark_theme_grey()

############################################################
#### Putting Flood and Rain together
############################################################
# For simplicity I will rename the data here
dat_rain <- catch %>%
  as.data.frame() %>%
  dplyr::select(Date, Value = Rain) %>%
  mutate(Data = "Rain") %>%
  mutate(NormalizedValue = Value / max(Value))

dat_flood <- dat_flood %>%
  as.data.frame() %>%
  dplyr::select(Date, Value = Flood) %>%
  mutate(Data = "Flood") %>%
  mutate(NormalizedValue = Value / max(Value))

# Put them together
dat <- rbind(dat_rain, dat_flood)

# Let's plot them as time series
ggplot(dat, aes(x = Date, y = NormalizedValue, color = Data)) +
  geom_line() +
  scale_color_viridis_d(begin = 0.7) +
  dark_theme_grey() +
  facet_grid("Data")

# Let's plot them as time series (without facets)
ggplot(dat, aes(x = Date, y = NormalizedValue, color = Data)) +
  geom_line() +
  scale_color_viridis_d(begin = 0.7) +
  dark_theme_grey()

# Let's look at some subsets to see more details
subset(dat, year(Date) %in% c(2000:2005)) %>%
  ggplot(aes(x = Date, y = NormalizedValue, color = Data)) +
    geom_line() +
    scale_color_viridis_d(begin = 0.7) +
    dark_theme_grey()

subset(dat, year(Date) %in% c(2006:2010)) %>%
  ggplot(aes(x = Date, y = NormalizedValue, color = Data)) +
    geom_line() +
    scale_color_viridis_d(begin = 0.7) +
    dark_theme_grey()

subset(dat, year(Date) %in% c(2011:2015)) %>%
  ggplot(aes(x = Date, y = NormalizedValue, color = Data)) +
    geom_line() +
    scale_color_viridis_d(begin = 0.7) +
    dark_theme_grey()

subset(dat, year(Date) %in% c(2016:2020)) %>%
  ggplot(aes(x = Date, y = NormalizedValue, color = Data)) +
    geom_line() +
    scale_color_viridis_d(begin = 0.7) +
    dark_theme_grey()

############################################################
#### Finding Peaks
############################################################
# Let's find peaks in the rain data
peaks_rain <- findpeaks(dat_rain$Value, minpeakdistance = 5)
peaks_rain <- dat_rain[peaks_rain[, 2], ]

# Visualize them
ggplot(dat_rain, aes(x = Date, y = Value)) +
  geom_line(col = viridis(20)[20]) +
  geom_point(data = peaks_rain, aes(x = Date, y = Value), col = viridis(20)[16]) +
  dark_theme_grey()

# We can also check during which months most of the peaks occur
table(month(peaks_rain$Date, label = T))

# To identify the peaks in the flood data we first need to interpolate our data
# to get a regular time series (otherwise choosing distances between two peaks
# does not make sense). Let's prepare an empty dataframe to fill
dat_flood_inter <- data.frame(
    Date = seq(min(dat_flood$Date), max(dat_flood$Date), by = "days")
  , Value = NA
)

# Interpolate data linearly
dat_flood_inter$Value <- approx(
    dat_flood$Date
  , dat_flood$Value
  , n = nrow(dat_flood_inter)
)$y

# Let's plot the data to make sure its correct
plot(Value ~ Date, data = dat_flood, type = "l")
lines(Value ~ Date, data = dat_flood_inter, col = "red", lty = 3)

# Compare the number of observations
nrow(dat_flood)
nrow(dat_flood_inter)

# Because we need one peak per year, we will split the data by year
extrema <- dat_flood_inter %>%
  mutate(Year = year(Date)) %>%
  group_by(Year) %>%
  nest() %>%
  mutate(Peaks = map(data, function(x){
    peaks <- which.max(x$Value)
    peaks <- x[peaks, ]
    return(peaks)
  })) %>%
  mutate(Valleys = map(data, function(x){
    valleys <- which.min(x$Value)
    valleys <- x[valleys, ]
    return(valleys)
  }))
peaks_flood <- do.call(rbind, extrema$Peaks)
valleys_flood <- do.call(rbind, extrema$Valleys)

# Let's look at them
arrange(peaks_flood, Date)
arrange(valleys_flood, Date)

# Visualize them
ggplot(dat_flood_inter, aes(x = Date, y = Value)) +
  geom_line(col = viridis(20)[20]) +
  geom_point(data = peaks_flood, aes(x = Date, y = Value), col = viridis(20)[16]) +
  geom_point(data = valleys_flood, aes(x = Date, y = Value), col = viridis(20)[10]) +
  dark_theme_grey()

# We can also check during which months most of the peaks occur
table(month(peaks_flood$Date, label = T))
table(month(valleys_flood$Date, label = T))

############################################################
#### Identify Lag Between Rain and Flood
############################################################
# Again, we need regular time series, so we will interpolate both, the rain and
# flood data to daily values
dat_flood_inter <- data.frame(
    Date = seq(min(dat_flood$Date), max(dat_flood$Date), by = "days")
  , Value = NA
)
dat_flood_inter$Value <- approx(
    dat_flood$Date
  , dat_flood$Value
  , n = nrow(dat_flood_inter)
)$y
dat_rain_inter <- data.frame(
    Date = seq(min(dat_rain$Date), max(dat_rain$Date), by = "days")
  , Value = NA
)
dat_rain_inter$Value <- approx(
    dat_rain$Date
  , dat_rain$Value
  , n = nrow(dat_rain_inter)
)$y

# Let's resample the data to monthly data
dat_flood_inter <- dat_flood_inter %>%
  mutate(Date = update(Date, day = 15)) %>%
  group_by(Date) %>%
  summarize(Value = mean(Value))
dat_rain_inter <- dat_rain_inter %>%
  mutate(Date = update(Date, day = 15)) %>%
  group_by(Date) %>%
  summarize(Value = mean(Value))

# # We only keep dates that are present in both datasets
# dat_flood_inter <- subset(dat_flood_inter
#   , Date %in% dat_flood_inter$Date
#   & Date %in% dat_rain_inter$Date
# )
# dat_rain_inter <- subset(dat_rain_inter
#   , Date %in% dat_flood_inter$Date
#   & Date %in% dat_rain_inter$Date
# )

# Make sure the everything aligns
nrow(dat_flood_inter)
range(dat_flood_inter$Date)
nrow(dat_rain_inter)
range(dat_rain_inter$Date)

# Coerce data to timeseries objects
ts_flood <- ts(dat_flood_inter$Value, start = c(2000, 4), frequency = 12)
ts_rain <- ts(dat_rain_inter$Value, start = c(2000, 4), frequency = 12)

# Plot them
plot(ts_flood)
plot(ts_rain)

# Calculate cross correlation
ccf <- ccf(ts_rain, ts_flood, lag.max = 9)
ccf

# Visualize
lag2.plot(ts_rain, ts_flood, max.lag = 9)

# Calculate lagged rain time series
ts_rain_lag0 <- stats::lag(ts_rain, 0)
ts_rain_lag1 <- stats::lag(ts_rain, -1)
ts_rain_lag2 <- stats::lag(ts_rain, -2)
ts_rain_lag3 <- stats::lag(ts_rain, -3)
ts_rain_lag4 <- stats::lag(ts_rain, -4)
ts_rain_lag5 <- stats::lag(ts_rain, -5)
ts_rain_lag6 <- stats::lag(ts_rain, -6)
ts_rain_lag7 <- stats::lag(ts_rain, -7)
ts_rain_lag8 <- stats::lag(ts_rain, -8)
ts_rain_lag9 <- stats::lag(ts_rain, -9)

# Combine data
all_dat <- ts.intersect(ts_flood
  , ts_rain_lag0
  , ts_rain_lag1
  , ts_rain_lag2
  , ts_rain_lag3
  , ts_rain_lag4
  , ts_rain_lag5
  , ts_rain_lag6
  , ts_rain_lag7
  , ts_rain_lag8
  , ts_rain_lag9
)

# Run a model
mod <- lm(ts_flood ~
  + ts_rain_lag0
  + ts_rain_lag1
  + ts_rain_lag2
  + ts_rain_lag3
  + ts_rain_lag4
  + ts_rain_lag5
  + ts_rain_lag6
  + ts_rain_lag7
  + ts_rain_lag8
  + ts_rain_lag9
  , data = all_dat
)

# Look at the results
summary(mod)

# End cluster
endCluster()
