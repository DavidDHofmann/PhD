############################################################
#### Cleaning and Preparing Tracks from GPS Fixes
############################################################
# Description: In this script I put together all gps fixes, regardless of
# whether an individual is a disperser or not. I then clean the data and add
# additional information, such as dispersal dates and pack-affiliation.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/Schreibtisch/15. PhD/Chapter_1"
setwd(wd)

# Load packages
library(tidyverse)  # For data wrangling
library(lubridate)  # To handle dates easily
library(pbmcapply)  # To make use of multiple cores
library(davidoff)   # Custom functions
library(raster)     # To handle spatial data
library(rgdal)      # To handle spatial data

############################################################
#### Data Source 1: POPECOL
############################################################
# Identify all files containing GPS data
files <- dir(
    path        = "03_Data/01_RawData/POPECOL"
  , pattern     = "GPS_Collar"
  , full.names = T
)

# There are some inconsistencies we need to take care of. Firstly, sometimes "."
# is used as decimal, sometimes ",". Moreover, there are some files containing
# "°" symbols. We need to get rid of this stuff. This only needs to be run once
# as it overwrites the default csv files. This also makes sure that every file
# is stored exactly the same way
# lapply(files, function(x){
#   dat <- read_file(x, local = locale(encoding = "latin1"))
#   dat <- gsub(dat, pattern = ",", replacement = ".")
#   dat <- gsub(dat, pattern = "°", replacement = "")
#   dat <- gsub(dat, pattern = "/", replacement = ".")
#   write(dat, x)
# })

# Load all of them and do some cleaning
dat1 <- lapply(files, function(x){

  # Load data...
  dat <- x %>%

    # ...using readr
    read_delim(local = locale(encoding = "latin1"), delim = ";") %>%

    # Retrieve DogName from filename. For this we use regex to identify a
    # pattern preceeded by five digits (?<=\\d{5}) and a "_", but is then
    # followed by a "_" and whatever
    mutate(DogName = str_extract(x, pattern = "(?<=\\d{5}_).*(?=_)")) %>%

    # Retrieve timestamp
    mutate(Timestamp = as.POSIXct(
      paste(UTC_Date, UTC_Time), tz = "UTC", format = "%d.%m.%Y %H:%M:%S")
    ) %>%

    # Remove stuff like [°]
    setNames(gsub(names(.), pattern = " \\[*", replacement = "")) %>%
    setNames(gsub(names(.), pattern = "\\]", replacement = "")) %>%

    # Keep only desired columns
    select(.
      , DogName   = DogName
      , DOP       = DOP
      , CollarID  = CollarID
      , x         = Longitude
      , y         = Latitude
      , Timestamp = Timestamp
    ) %>%

    # Add a column indicating the data source
    mutate(Source = "Popecol")
}) %>% do.call(rbind, .)

# Check validity of columns
summary(dat1$x)
summary(dat1$y)
sum(is.na(dat1$Timestamp))
sum(is.na(dat1$CollarID))
range(dat1$Timestamp)

# The GPS devices typcially record GPS locations even before and after collaring
# an animal. These GPS fixes need to be removed since they provide wrong
# location data. To do so, we load the dataframe that Dominik prepared. It shows
# for each collar and dog the first and last fix retrieved, i.e. when an
# individual was collared or uncollared. The table also contains valuable
# information about pack-affiliations
collars <- read_csv2("03_Data/01_RawData/POPECOL/CollarSettings.csv") %>%

  # Remove rows where the Dog Names are NA since these rows do not contain any
  # useful information
  subset(., !is.na(`Dog Name`)) %>%

  # Keep only the desired columns and rename them nicely. Note that we want to
  # keep the Dog Code because it provides information about the birth pack of
  # each dog.
  select(.
    , CollarID  = `Collar Nr.`
    , DogName   = `Dog Name`
    , DogCode   = `Dog Code`
    , Sex       = `Sex`
    , FirstDate = `Collaring Date`
    , LastDate1 = `Stop recording date`
    , LastDate2 = `Last fix date`
  )

# Looking at the table we can see that we are missing exact times for some of
# the timestamps. This will cause errors when we convert those to posixct. We
# Will therefore assign very conservative times (i.e. 24:00 for first dates,
# 00:01 for last dates)
for (i in 1:nrow(collars)){
  if (!is.na(collars$FirstDate[i]) & nchar(collars$FirstDate[i]) == 10){
    collars$FirstDate[i] <- paste0(collars$FirstDate[i], " 24:00")
  }
  if (!is.na(collars$LastDate1[i]) & nchar(collars$LastDate1[i]) == 10){
    collars$LastDate1[i] <- paste0(collars$LastDate1[i], " 00:01")
  }
  if (!is.na(collars$LastDate2[i]) & nchar(collars$LastDate2[i]) == 10){
    collars$LastDate2[i] <- paste0(collars$LastDate2[i], " 00:01")
  }
}

# Taryn's collar number is also entered incorrectly
collars$CollarID[collars$DogName == "Taryn" & collars$CollarID == 22028] <- 20228

# Now we can coerce the date columns to true dates. Note that we subtract two
# hours because the original format was in local time (which is UTC + 2 hours).
# Also note that there are two possible "last dates". The first one refers to
# the last fix, the other to the date when the dog was uncollared. We need to
# combine both because depending on whether the dog is still collared we will
# use on or the other. Again we subtract two hours to make sure that the times
# are in UTC
collars <- collars %>%
  mutate(.
    , FirstDate = as.POSIXct(FirstDate
      , tz = "UTC"
      , format = "%d.%m.%Y %H:%M") - hours(2)
    , LastDate1 = as.POSIXct(LastDate1
      , tz = "UTC"
      , format = "%d.%m.%Y %H:%M") - hours(2)
    , LastDate2 = as.POSIXct(LastDate2
      , tz = "UTC"
      , format = "%d.%m.%Y %H:%M") - hours(2)
  )

# Merge the columns for the "LastDates". We only want to keep the latest of the
# two.
collars$LastDate <- pmax(collars$LastDate1, collars$LastDate2, na.rm = T)
collars$LastDate <- collars$LastDate + minutes(5)
collars$LastDate1 <- NULL
collars$LastDate2 <- NULL

# Logical check that LastDate > FirstDate
table(collars$LastDate > collars$FirstDate)

# Make sure there are no duplicates
table(paste(collars$CollarID, collars$DogName))
table(table(paste(collars$CollarID, collars$DogName)))

# Left-join the dataset to the gps locations
dat1 <- left_join(dat1, collars, by = c("CollarID", "DogName"))

# Create a column that indicates if the respective fix lies within the first and
# last date
dat1$Keep <- ifelse(
    test  = dat1$Timestamp >= dat1$FirstDate & dat1$Timestamp <= dat1$LastDate
  , yes   = T
  , no    = F
)

# Let's check how many datapoints we loose (per dog and collar)
table(dat1$Keep)
table(dat1$DogName, dat1$Keep)
table(dat1$CollarID, dat1$Keep)

# Subset data and remove unnecessary columns
dat1 <- subset(dat1, Keep)
dat1 <- select(dat1, -c("FirstDate", "LastDate", "Keep", "DogCode"))

# Prepare data for visualization
vis <- dat1 %>% group_by(DogName, CollarID) %>% summarize(
    First = range(Timestamp)[1]
  , Last  = range(Timestamp)[2]
)

# Verify that collars are not overlapping
ggplot(vis, aes(color = factor(CollarID))) +
  geom_segment(
      aes(x = First, xend = Last, y = DogName, yend = DogName)
    , size = 3
    , alpha = 0.6
  )

# Maybe look at some in more detail
ggplot(subset(vis, DogName == "Taryn"), aes(color = factor(CollarID))) +
  geom_segment(
      aes(x = First, xend = Last, y = DogName, yend = DogName)
    , size = 10
    , alpha = 0.6
  )

################################################################################
#### Data Source 2: Abrahms
################################################################################
# We also want to import the GPS fixes provided by Abrahms. Note that I will not
# clean this data anymore, as this is the cleaned data that Abrahms used for her
# publication
dat2 <- read_csv("03_Data/01_RawData/ABRAHMS/DispersalPaths.csv") %>%

  # Remove rows with missing fixes
  filter(., !is.na(`Longitude..deg.`)) %>%

  # Select only the columns of interest
  select(.,
      DogName   = `id`
    , x         = `Longitude..deg.`
    , y         = `Latitude..deg.`
    , Timestamp = `Timestamp`
  ) %>%

  # Add all columns that are also in dat1
  mutate(., CollarID = NA, DOP = NA, Sex = "M", Source = "Abrahms")

################################################################################
#### Data Source 3: Botswana Camp
################################################################################
# Lastly, we retrieved fixes from the camp in Botswana that were collected prior
# to own our project in Botswana. This is data from residents only! I only
# include it for completeness. This data requires some intensive cleaning
files <- dir(
    path        = "03_Data/01_RawData/DOMINIK"
  , pattern     = ".txt$"
  , full.names  = T
)

# Load all the files
dat3 <- lapply(files, function(x){

  # Identify the DogName
  name <- substr(basename(x), start = 6, stop = nchar(basename(x)) - 7)

  # Load the data (skip the first two rows as they contain stuff we dont care
  # about)
  read_delim(x, delim = "\t", skip = 2) %>%

    # Keep only the columns of interest
    select(
        x                   = `Longitude (deg)`
      , y                   = `Latitude (deg)`
      , Timestamp           = `UTC time (yyyy-mm-dd HH:MM:SS)`
      , Height              = `Height above MSL (m)`
      , HorizontalAccuracy  = `Horizontal accuracy (m)`
      , VerticalAccuracy    = `Vertical accuracy (m)`
    ) %>%

    # Add dog name
    mutate(DogName = name) %>%

    # Prepare columns present in dat1
    mutate(., CollarID = NA, DOP = NA, Sex = "M", Source = "Botswana") %>%

    # Remove rows with missing fixes
    filter(., !is.na(x)) %>%

    # Order the data by timestamp
    arrange(., Timestamp) %>%

    # Create Column that indicates the timelag and distance between two fixes
    # This allows us to also calculate the speed
    mutate(x,
        dt = as.numeric(Timestamp - lag(Timestamp), units = "hours")
      , dl = sqrt((x - lag(x))**2 + (y - lag(y)) ** 2) * 111
      , speed = dl / dt
    )
}) %>% do.call(rbind, .)

# There are some pretty beefy outliers that we need to remove. We can use the
# Accuracy indicators, the height and speed in order to remove any unreasonable
# gps fix (we are going to be pretty rough here)
dat3 <- dat3 %>%
  filter(., Height > quantile(Height, 0.1)) %>%
  filter(., Height < quantile(Height, 0.9)) %>%
  filter(., VerticalAccuracy < quantile(VerticalAccuracy, 0.9)) %>%
  filter(., HorizontalAccuracy < quantile(HorizontalAccuracy, 0.9)) %>%
  filter(., year(Timestamp) != 1970) %>%
  filter(., speed < quantile(speed, 0.9, na.rm = TRUE)) %>%
  select(., c("DogName", "CollarID", "x", "y", "Timestamp", "DOP", "Sex", "Source"))

# Put data of all sources together
data <- rbind(dat1, dat2, dat3)

# Check sources
table(data$Source)

# Check if there are any duplicates
dups_complete <- duplicated(data)
sum(dups_complete)

# Store data for exploration
write_csv(data, "GPS_Data.csv")

################################################################################
#### Adding Cutoff Dates
################################################################################
# For each GPS fix we want to know whether it was collected during dispersal or
# residency. We can use the following table for this.
cut <- read_csv("03_Data/01_RawData/POPECOL/CutoffDates.csv") %>%

  # Remove undesired columns
  select(c("CollarID", "DogName", "StartDate", "EndDate")) %>%

  # Add Minutes to the timestamps
  mutate(
      StartDate = paste0(StartDate, ":00")
    , EndDate   = paste0(EndDate, ":00")
  ) %>%

  # Make proper dates and subtract 2 hours
  mutate(
      StartDate = as.POSIXct(StartDate, tz = "UTC"
        , format = "%d.%m.%Y %H:%M") - hours(2)
    , EndDate = as.POSIXct(EndDate, tz = "UTC"
      , format = "%d.%m.%Y %H:%M") - hours(2)
  ) %>%

  # Sort
  arrange(DogName, StartDate)

# As you can see for some individuals there are multiple dispersal phases. Let's
# create a counter for each individual
cut <- cut %>%
  group_by(DogName) %>%
  mutate(DispersalNo = row_number()) %>%
  ungroup()

# Store the merged cutoff dates
write.csv(cut, "03_Data/02_CleanData/00_General_Dispersers_Popecol_CutoffDates.csv")

# We can use the table to create a new column that indicates the state of the
# animal. By default we say that our individuals are resident.
data$State <- "Resident"

# Now loop through all dogs and use the cutoff table to specify whether a fix
# was taken during dispersal or not
names <- unique(cut$DogName)
for (i in seq_along(names)){
  cutoff <- subset(cut, DogName == names[i])
  index <- which(data$DogName == names[i])
  for (h in 1:nrow(cutoff)){
    data$State[index][data$Timestamp[index] >= cutoff$StartDate[h] &
    data$Timestamp[index] <= cutoff$EndDate[h]] <- "Disperser"
  }
}

# Look at the final table
head(data)
tail(data)

# Check the number of data per individual
table(data$DogName, data$State)

# ############################################################
# #### Adding Pack Information
# ############################################################
# # Let's load the csv containing all pack affiliations
# packs <- read_csv("03_Data/01_RawData/POPECOL/PackAffiliation.csv")
#
# # Create a column that indicates if a fixe was recorder before or after
# # dispersal
# data$AfterDispersal <- ifelse(data$Timestamp > data$EndDate, T, F)
#
# # Lets join the birth and postdispersal pack to each fix
# data <- left_join(data, packs, by = "DogName")
#
# # During dispersal the individuals are not with their pack anymore so let's put
# # NAs for the pack during dispersal
# data$BirthPack[data$State == "Disperser"] <- NA
#
# # After dispersal (i.e. settlement) the individuals will be in a new pack. Let's
# # create a column that indicates whether a GPS fix was taken before or after
# # dispersal so that we can change the packs accordingly. Note that there are
# # individuals with multiple dispersal phases. We will only use the most recent
# # EndDate for them (we assume they haven't settled in a new pack before this
# # date).
# enddates <- cut %>%
#
#   # Group the enddates by dog
#   group_by(DogName) %>%
#
#   # Keep only the latest date
#   slice(which.max(EndDate)) %>%
#
#   # Remove the start date from the resulting table
#   select(., -StartDate)
#
# # Join these dates to our dataframe
# data <- left_join(data, enddates, by = "DogName")
#
# # Finally we merge the "BirthPack" and "PostDispersalPack" columns and create
# # a column that indicates the current pack
# data$CurrentPack <- data$BirthPack
# data$CurrentPack[data$AfterDispersal] <-
#   data$PostDispersalPack[data$AfterDispersal]
#
# # Remove any undesired column
# data <- select(data, -c("BirthPack", "PostDispersalPack", "EndDate"))
#
# # Now there are some manual adjustments that I do according to Dominik's
# # instructions. For Fatalii the pack is unknown before the 18.02.2016
# data$CurrentPack[data$DogName == "Fatalii" &
#   data$Timestamp < "2016-02-18 00:00:00"] <- NA
#
# # Mirage and Kalahri belong to BN only after the 20.09.2017. Let's make anything
# # after dispersal NA. Then add the BN pack after the 20-09-2017
# data$CurrentPack[data$AfterDispersal & data$DogName == "Kalahari"] <- NA
# data$CurrentPack[data$AfterDispersal & data$DogName == "Kalahari" &
#   data$Timestamp > "2017-09-20"] <- "BN"
# data$CurrentPack[data$AfterDispersal & data$DogName == "Mirage"] <- NA
# data$CurrentPack[data$AfterDispersal & data$DogName == "Mirage" &
#   data$Timestamp > "2017-09-20"] <- "BN"
#
# # Kalahari eventually switched the pack from BN to WA. This happened the same
# # day that Belgium switched to WA. We can use our cutoff data table to find this
# # date
# date <- cut$EndDate[cut$DogName == "Belgium"]
#
# # Let's change Kalaharis pack after this timestamp to WA
# data$CurrentPack[data$Timestamp > date & data$DogName == "Kalahari"] <- "WA"
#
# # Appalachia's "after-dispersal-pack" is called APL
# data$CurrentPack[data$AfterDispersal & data$DogName == "Appalachia"] <- "APL"
#
# # Bongwe eventually dispersed socially
# data$CurrentPack[data$DogName == "Bongwe"
#   & data$Timestamp > "2012-09-22 00:00:00"] <- "HU"
#
# # Remove unnecessary columns
# data <- select(data, -c("AfterDispersal"))
#
# # Order the dataframe according to: DogName, Timestamp
# data <- arrange(data, DogName, Timestamp)

################################################################################
#### Resampmling Data
################################################################################
# Remove NA fixes
table(is.na(data$x))
table(is.na(data$y))
data <- subset(data, !is.na(x) & !is.na(y))

# Sort data
data <- arrange(data, DogName, Timestamp)

# Some of the data has a very high resolution that we don't need. We will
# therefore subsample to a resolution we can work with.
backup <- data
data <- data %>% group_by(DogName) %>% nest()
data$data <- suppressMessages(
  pbmclapply(data$data
    , ignore.interactive =  T
    , mc.cores = detectCores() - 1
    , FUN = function(x){resFix2(x, hours = 1, start = 1)}
  )
)
data <- unnest(data)

############################################################
#### Store the Output
############################################################
# Write the data to file
write_csv(data, "03_Data/02_CleanData/00_General_Dispersers_Popecol.csv")

# Create a shapefile for visualization
data <- data %>% group_by(DogName) %>% nest()
tracks <- pbmclapply(1:nrow(data), ignore.interactive = T, mc.cores = detectCores() - 1, function(x){
  coords <- data$data[[x]]
  coords$DogName <- data$DogName[x]
  x <- coords
  coordinates(x) <- c("x", "y")
  crs(x) <- CRS("+init=epsg:4326")
  lines <- spLines(x)
  lines <- createSegments(lines)
  lines <- as(lines, "SpatialLinesDataFrame")
  lines@data <- x@data[1:(nrow(x) - 1), ]
  return(lines)
}) %>% do.call(rbind, .)

# Store them
writeOGR(tracks
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_Dispersers_Popecol"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)
