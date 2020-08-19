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
library(tidyverse)

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
# as it overwrites the default csv files.
# lapply(files, function(x){
#   dat <- read_file(x, local = locale(encoding = "latin1"))
#   dat <- gsub(dat, pattern = ",", replacement = ".")
#   dat <- gsub(dat, pattern = "°", replacement = "")
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

# This data still contains periods prior and after collaring of the animal. We
# need to get rid of such data.






# Visualize data
vis <- dat1 %>% group_by(DogName, CollarID) %>% summarize(
    First = range(Timestamp)[1]
  , Last = range(Timestamp)[2]
)

ggplot(dat1, aes(x = x, y = y, col = factor(DogName)))

ggplot(vis, aes(color = factor(CollarID))) +
  geom_segment(aes(x = First, xend = Last, y = DogName, yend = DogName), size = 2)

ggplot(subset(vis, DogName == "Odzala"), aes(color = factor(CollarID))) +
  geom_segment(aes(x = First, xend = Last, y = DogName, yend = DogName), size = 2, alpha = 0.6)










# We also want to import the GPS fixes provided by Abrahms
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

  # Add a column for the (nonexistent) collar id so that the dataframe contains
  # the same columns as the dat1 dataframe
  mutate(., CollarID = NA)

# Lastly, we retrieved fixes from the camp in Botswana that were collected
# prior to own our project in Botswana. This data requires some intensive
# cleaning
files <- dir(
    path        = "03_Data/01_RawData/DOMINIK"
  , pattern     = ".txt$"
  , full.names  = T
)

# Identify the DogNames from the file descriptions
names <- substr(files, start = 33, stop = nchar(files) - 7)

# Load all the files
dat3 <- lapply(files, function(x){

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
  )}) %>%

  # The resulting object is a list with one entry for each dog. So let's add the
  # dog names that we retrieved from the file descriptions
  set_names(., names) %>%

  # Bind the dataframes together, but create a new column indicating the
  # individuals (using idcol = TRUE)
  rbindlist(., idcol = TRUE) %>%

  # Make a nicer column name
  rename(., DogName = .id) %>%

  # Prepare an empty column for the (for dat3 nonexistent) CollarID
  mutate(., CollarID = NA) %>%

  # Remove rows with missing fixes
  filter(., !is.na(x)) %>%

  # Order the data by dogname and timestamp
  arrange(., DogName, Timestamp) %>%

  # Split dataframe into list by DogName
  split(., .$DogName) %>%

  # Create Column that indicates the timelag and distance between two fixes
  # This allows us to also calculate the speed
  lapply(., function(x){
    mutate(x,
        dt = as.numeric(Timestamp - lag(Timestamp), units = "hours")
      , dl = sqrt((x - lag(x))**2 + (y - lag(y)) ** 2) * 111
      , speed = dl / dt
    )
  }) %>%

  # Collapse list to single dataframe
  do.call(rbind, .)

# There are some pretty beefy outliers that we need to remove. We can use the
# Accuracy indicators, the height and speed in order to remove any unreasonable
# gps fix
dat3 <- dat3 %>%

  # Remove fixes with unreasonable high elevation
  filter(., Height > quantile(Height, 0.1)) %>%

  # Remove fixes with unreasonable low elevation
  filter(., Height < quantile(Height, 0.9)) %>%

  # Remove fixes with too low vertical accuracy
  filter(., VerticalAccuracy < quantile(VerticalAccuracy, 0.9)) %>%

  # Remove fixes with too low horizontal accuracy
  filter(., HorizontalAccuracy < quantile(HorizontalAccuracy, 0.9)) %>%

  # Remove fixes with wrong dates
  filter(., year(Timestamp) != 1970) %>%

  # Remove fixes with unreasonable speed
  filter(., speed < quantile(speed, 0.9, na.rm = TRUE)) %>%

  # Finally, we remove all columns that we don't need furthermore
  select(., c("DogName", "CollarID", "x", "y", "Timestamp"))

# Combine the three datasets
data <- rbind(dat1, dat2, dat3)

# Remove any remaining duplicates
data <- data[!duplicated(data), ]

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
#### Cut Periods of GPS Handling
############################################################
# The GPS devices typcially record GPS locations even before and after collaring
# an animal. These GPS fixes need to be removed since they provide wrong
# location data. To do so, we load the dataframe that Dominik prepared. It shows
# for each collar and dog the first and last fix retrieved, i.e. when an
# individual was collared or uncollared. The table also contains valuable
# information about pack-affiliations
collars <- read_csv("03_Data/01_RawData/POPECOL/CollarSettings.csv") %>%

  # Remove rows where the Dog Names are NA since these rows do not contain any
  # useful information
  subset(., !is.na(`Dog Name`)) %>%

  # Keep only the desired columns and rename them nicely. Note that we want to
  # keep the Dog Code because it provides information about the birth pack of
  # each dog.
  select(.
    , CollarNr  = `Collar Nr.`
    , DogName   = `Dog Name`
    , DogCode   = `Dog Code`
    , Sex       = `Sex`
    , FirstDate = `Collaring Date`
    , LastDate1 = `Stop recording date`
    , LastDate2 = `Last fix date`
  ) %>%

  # Coerce date columns to true dates. Note that we subtract two hours because the
  # original format was not in UTC time but rather local time (which is UTC + 2
  # hours). Also note that there are two possible "last dates". The first one
  # refers to the last fix, the other to the date when the dog was uncollared.
  # We need to combine both because depending on whether the dog is still
  # collared one or the other will be used as true "last date". Again we subtract
  # two hours to make sure that the times are in UTC
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
collars$LastDate <- collars$LastDate1
for (i in 1:nrow(collars)){
  collars$LastDate[i] <- max(collars$LastDate1[i], collars$LastDate2[i], na.rm = TRUE)
}

# Remove any undesired column
collars <- select(collars, -c("LastDate1", "LastDate2"))

# There was also an issue with the first dates for Belgium and Abel. Let's add
# them back manually. Note that the system will automatically subtract 2 hours
# from the times that we enter. I.e. the times will be 22:00 afterwards. (UTC
# time)
collars$FirstDate[collars$DogName == "Belgium" &
  is.na(collars$FirstDate)] <- "2017-10-28 24:00:00"
collars$FirstDate[collars$DogName == "Abel" &
  is.na(collars$FirstDate)] <- "2017-10-26 24:00:00"

# We can now get rid of GPS locations that are not valid because the locations
# were recorded before the dog was actually collared. Thus, for each dog we
# remove the GPS fixes that were collected prior to the collaring date. Let's
# get a table that shows the "first" date and "end" date for each dog and collar
dates <- select(collars, c(DogName, CollarID = CollarNr, FirstDate, LastDate))

# Left-join the dataset to the dataset containing the gps locations
data <- left_join(data, dates, by = c("CollarID", "DogName"))

# For the dogs that we did not collar ourselves (Abrahms dogs and the dogs
# collared before 2016) there is no such information. This might cause some
# problems later. We can easily identify these Dogs as their CollarID is NA. We
# anyways don't want to cut any of their locations so let's just use each dog's
# first and last recording date as substitutes.
temp <- unique(data$DogName[is.na(data$CollarID)])
for (i in 1:length(temp)){
  sub <- subset(data, DogName == temp[i] & is.na(CollarID))
  FirstDate <- range(sub$Timestamp)[1]
  LastDate  <- range(sub$Timestamp)[2]
  data$FirstDate[data$DogName == temp[i] & is.na(data$CollarID)] <- FirstDate
  data$LastDate[data$DogName == temp[i] & is.na(data$CollarID)] <- LastDate
}

# Make sure there are no NA's
sum(is.na(data$FirstDate))
sum(is.na(data$LastDate))

# Remove rows for which the recording date does not lie between the First and
# Last Date. In other words, we get rid of the recordings that are before or
# after the first and last dates. Let's get the indices of these rows
indices <- which(data$Timestamp < data$FirstDate | data$Timestamp > data$LastDate)
data <- data[-indices, ]

# Now remove the FirstDate and LastDate columns since they are not needed
# anymore
data <- select(data, -c("FirstDate", "LastDate"))
data %>%
  group_by(DogName) %>%
  summarize(FirstDate = min(Timestamp), LastDate = max(Timestamp))

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
