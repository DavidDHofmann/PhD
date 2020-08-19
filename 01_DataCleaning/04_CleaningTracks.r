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
    )
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
  mutate(., CollarID = NA, DOP = NA, Sex = "M")

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
    mutate(., CollarID = NA, DOP = NA, Sex = "M") %>%

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
  select(., c("DogName", "CollarID", "x", "y", "Timestamp", "DOP", "Sex"))
data <- rbind(dat1, dat2, dat3)
data <- data[!duplicated(data), ]

################################################################################
#### CONTINUE HERE!!!
################################################################################
############################################################
#### Some Exploration
############################################################
# Check number of Individuals
length(unique(data$DogName))

# Number of GPS fixes
nrow(data)

# Number of GPS fixes per Dog. You will find that there are vastly different
# amounts of fixes. This is mostly due to different fix-rates
data %>%
  group_by(., DogName) %>%
  summarize(., NoFixes = n()) %>%
  arrange(., NoFixes)

# Check the range of observation for each individual
data %>%
  group_by(., DogName) %>%
  summarize(., FirstDate = min(Timestamp), LastDate = max(Timestamp))

############################################################
#### Adding Cutoff Dates
############################################################
# For each GPS fix we want to know whether it was collected during dispersal or
# residency. For this we want to create a table that tells us the cutoff dates
# for each individual. To create this table we have to combine several sources.

# Let's load the first cutoff-date table, which contains most of the cutoff
# dates. We subtract two hours from the timestamp to get UTC times
cut1 <- read_csv("03_Data/01_RawData/POPECOL/CutoffDates1.csv") %>%

  # Set some nice column names
  set_names(., c("DogName", "StartDate", "EndDate", "Comments")) %>%

  # Remove comments
  select(., -c(Comments)) %>%

  # Coerce the StartDate to a POSIXct column
  mutate(., StartDate = as.POSIXct(paste(StartDate, "00:00:00")
    , tz = "UTC"
    , format = "%Y-%m-%d %H:%M:%S") - hours(2)) %>%

  # Coerce the EndDate to a POSIXct column
  mutate(., EndDate = as.POSIXct(paste(EndDate, "24:00:00")
    , tz = "UTC"
    , format = "%Y-%m-%d %H:%M:%S") - hours(2))

# Now we load the second cutoff-date table, which contains cutoff dates for
# those indvidiuals with multiple dispersal phases. These individuals left their
# pack but eventually returned before fully dispersing. Again we subtract two
# hours to get UTC times.
cut2 <- read_csv("03_Data/01_RawData/POPECOL/CutoffDates2.csv") %>%

  # Select only the columns of interest
  select(.,
      DogName   = `DogName`
    , StartDate = `StartDate`
    , EndDate   = `EndDate`
  ) %>%

  # Coerce the StartDate to a POSIXct column
  mutate(., StartDate = as.POSIXct(StartDate
    , tz = "UTC"
    , format = "%d.%m.%Y %H:%M") - hours(2)) %>%

  # Coerce the EndDate to a POSIXct column
  mutate(., EndDate = as.POSIXct(EndDate
    , tz = "UTC"
    , format = "%d.%m.%Y %H:%M") - hours(2))


# Put all three tables together into one big cutoff date table
cut <- rbind(cut1, cut2)

# Look at the resulting table
cut

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

# Check the number of fixes during dispersal per dog
data %>%
  subset(State == "Disperser") %>%
  group_by(DogName) %>%
  summarise(noFixes = n())

# Compare the number of fixes during dispersal and residency
data %>%
  group_by(DogName, State) %>%
  summarize(NoFixes = n()) %>%
  spread(State, NoFixes)

############################################################
#### Adding Pack Information
############################################################
# Let's load the csv containing all pack affiliations
packs <- read_csv("03_Data/01_RawData/POPECOL/PackAffiliation.csv")

# Lets join the birth and postdispersal pack to each fix
data <- left_join(data, packs, by = "DogName")

# During dispersal the individuals are not with their pack anymore so let's put
# NAs for the pack during dispersal
data$BirthPack[data$State == "Disperser"] <- NA

# After dispersal (i.e. settlement) the individuals will be in a new pack. Let's
# create a column that indicates whether a GPS fix was taken before or after
# dispersal so that we can change the packs accordingly. Note that there are
# individuals with multiple dispersal phases. We will only use the most recent
# EndDate for them (we assume they haven't settled in a new pack before this
# date).
enddates <- cut %>%

  # Group the enddates by dog
  group_by(DogName) %>%

  # Keep only the latest date
  slice(which.max(EndDate)) %>%

  # Remove the start date from the resulting table
  select(., -StartDate)

# Join these dates to our dataframe
data <- left_join(data, enddates, by = "DogName")

# Create new column that indicates whether a fix was taken after dispersal or not
data$AfterDispersal <- FALSE
data$AfterDispersal[data$Timestamp > data$EndDate] <- TRUE

# Finally we merge the "BirthPack" and "PostDispersalPack" columns and create
# a column that indicates the current pack
data$CurrentPack <- data$BirthPack
data$CurrentPack[data$AfterDispersal] <-
  data$PostDispersalPack[data$AfterDispersal]

# Remove any undesired column
data <- select(data, -c("BirthPack", "PostDispersalPack", "EndDate"))

# Now there are some manual adjustments that I do according to Dominik's
# instructions. For Fatalii the pack is unknown before the 18.02.2016
data$CurrentPack[data$DogName == "Fatalii" &
  data$Timestamp < "2016-02-18 00:00:00"] <- NA

# Mirage and Kalahri belong to BN only after the 20.09.2017. Let's make anything
# after dispersal NA. Then add the BN pack after the 20-09-2017
data$CurrentPack[data$AfterDispersal & data$DogName == "Kalahari"] <- NA
data$CurrentPack[data$AfterDispersal & data$DogName == "Kalahari" &
  data$Timestamp > "2017-09-20"] <- "BN"
data$CurrentPack[data$AfterDispersal & data$DogName == "Mirage"] <- NA
data$CurrentPack[data$AfterDispersal & data$DogName == "Mirage" &
  data$Timestamp > "2017-09-20"] <- "BN"

# Kalahari eventually switched the pack from BN to WA. This happened the same
# day that Belgium switched to WA. We can use our cutoff data table to find this
# date
date <- cut$EndDate[cut$DogName == "Belgium"]

# Let's change Kalaharis pack after this timestamp to WA
data$CurrentPack[data$Timestamp > date & data$DogName == "Kalahari"] <- "WA"

# Appalachia's "after-dispersal-pack" is called APL
data$CurrentPack[data$AfterDispersal & data$DogName == "Appalachia"] <- "APL"

# Bongwe eventually dispersed socially
data$CurrentPack[data$DogName == "Bongwe"
  & data$Timestamp > "2012-09-22 00:00:00"] <- "HU"

# Remove unnecessary columns
data <- select(data, -c("AfterDispersal"))

# Order the dataframe according to: DogName, Timestamp
data <- arrange(data, DogName, Timestamp)

############################################################
#### Create Shapefile and Store the Output
############################################################
# As a last step before we store the data we might want to resample the
# remaining fixes to 2 hours such that there is only one fix every two hours.
# This will make sure that we don't have too many fixes
data <- data %>% group_by(DogName) %>% nest()
data$data <- suppressMessages(
  mclapply(data$data, mc.cores = detectCores() - 1, function(x){
    resFix2(x, hours = 2, start = 1)
  })
)
data <- unnest(data)

# We can use this data to create shapefiles that we can use for visual
# inspection. Let's create a track
tracks <- data %>% make_track(.
    , .x    = x
    , .y    = y
    , .t    = Timestamp
    , id    = DogName
    , state = State
    , pack  = CurrentPack
    , crs   = CRS("+init=epsg:4326")) %>%

  # Transform to utm
  transform_coords(CRS("+init=epsg:32734")) %>%

  # And nest the tracks
  nest(data = -"id") %>%

  # Turn to a step representation (note that the option "keep_cols" allows us to
  # retain the information about dispersers in the resulting object)
  mutate(., data = map(data, function(x){
    x %>% steps(., keep_cols = "start")
  })) %>%

  # Unnest the tibble
  unnest() %>%

  # Make the turning angles clockwise
  mutate(., ta_ = ta_ * (-1)) %>%

  # Add a column that indicates the absolute turning angle (heading)
  mutate(., absta_ = absAngle(.))

# Use our custom function to create spatial lines
lines <- lineTrack(tracks, CRS("+init=epsg:32734"))

# Coerce the lines to WGS84
lines <- spTransform(lines, CRS("+init=epsg:4326"))

# To store the shapefiles we need to coerce the difftame column to a numeric
# column
lines$dt_ <- as.numeric(lines$dt_)

# Store the merged dataframe and shapefile to the different locations
write_csv(data, "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv")
writeOGR(lines
  , "03_Data/02_CleanData"
  , "00_General_Dispersers_Popecol(Regular)"
  , driver = "ESRI Shapefile"
  , overwrite = TRUE
)
