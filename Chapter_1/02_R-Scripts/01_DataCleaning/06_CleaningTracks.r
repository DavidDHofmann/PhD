################################################################################
#### Cleaning and Preparing Tracks from GPS Fixes
################################################################################
# Description: In this script I put together all gps fixes, regardless of
# whether an individual is a disperser or not. I then clean the data and add
# additional information, such as dispersal dates and pack-affiliation.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load packages
library(tidyverse)  # For data wrangling
library(lubridate)  # To handle dates easily
library(pbmcapply)  # To make use of multiple cores
library(davidoff)   # Custom functions
library(raster)     # To handle spatial data
library(rgdal)      # To handle spatial data
library(plotKML)    # To store the cleaned tracks to a kml
library(spacetime)  # To store the cleaned tracks to a kml (required for slider)
library(viridis)    # For nice colors

################################################################################
#### Data Source 1: POPECOL
################################################################################
# Identify all files containing GPS data
files <- dir(
    path        = "03_Data/01_RawData/POPECOL/01_GPS"
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

    # Remove special characters like [°]
    setNames(gsub(names(.), pattern = " \\[*", replacement = "")) %>%
    setNames(gsub(names(.), pattern = "\\]", replacement = "")) %>%

    # Keep only desired columns
    dplyr::select(.
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

  # Remove rows where the Dog Names are NA
  subset(., !is.na(`Dog Name`)) %>%

  # Keep only the desired columns and rename them nicely. Note that we want to
  # keep the Dog Code because it provides information about the birth pack of
  # each dog.
  dplyr::select(.
    , CollarID  = `Collar Nr.`
    , DogName   = `Dog Name`
    , DogCode   = `Dog Code`
    , Sex       = `Sex`
    , FirstDate = `Collaring Date`
    , LastDate1 = `Stop recording date`
    , LastDate2 = `Last fix date`
  )

# We are missing exact times for some of the timestamps. This will cause errors
# when we convert those to posixct. We Will therefore assign very conservative
# times (i.e. 24:00 for first dates, 00:01 for last dates)
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

# Taryn's collar number is wrong
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

# Logical check that LastDate > FirstDate and make sure there are no duplicates
table(collars$LastDate > collars$FirstDate)
table(paste(collars$CollarID, collars$DogName))
table(table(paste(collars$CollarID, collars$DogName)))

# Create a column that indicates if the respective fix lies within the first and
# last date
dat1 <- left_join(dat1, collars, by = c("CollarID", "DogName"))
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
dat1 <- dplyr::select(dat1, -c("FirstDate", "LastDate", "Keep", "DogCode"))

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
sub <- c("Taryn", "Abel")
ggplot(subset(vis, DogName %in% sub), aes(color = factor(CollarID))) +
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
# publication. Her timestamps are in UTC already.
dat2 <- read_csv("03_Data/01_RawData/ABRAHMS/DispersalPaths.csv") %>%

  # Remove rows with missing fixes
  filter(., !is.na(`Longitude..deg.`)) %>%

  # Select only the columns of interest
  dplyr::select(.,
      DogName   = `id`
    , x         = `Longitude..deg.`
    , y         = `Latitude..deg.`
    , Timestamp = `Timestamp`
  ) %>%

  # Add columns that are also in dat1 (so we can bind the data afterwards). Note
  # that all individuals from Abrahms are males.
  mutate(., CollarID = NA, DOP = NA, Sex = "M", Source = "Abrahms")

################################################################################
#### Data Source 3: Botswana Camp
################################################################################
# Lastly, we retrieved fixes from the camp in Botswana that were collected prior
# to own our project in Botswana. This is data from residents only! I only
# include it for completeness. This data requires some intensive cleaning. Also
# note that the data is in UTC already.
files <- dir(
    path        = "03_Data/01_RawData/DOMINIK/GPS"
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
    dplyr::select(
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
  dplyr::select(.
    , c("DogName", "CollarID", "x", "y", "Timestamp", "DOP", "Sex", "Source")
  )

# Put data of all sources together
data <- rbind(dat1, dat2, dat3)

# Check sources
table(data$Source)

# Check if there are any duplicates
dups_complete <- duplicated(data)
sum(dups_complete)
dups_incomplete <- duplicated(data[, c("DogName", "Timestamp")])
sum(dups_incomplete)

# Remove them
data <- subset(data, !dups_incomplete)

################################################################################
#### Adding Cutoff Dates
################################################################################
# For each GPS fix we want to know whether it was collected during dispersal or
# residency. We can use the following table for this.
cut <- read_csv("03_Data/01_RawData/POPECOL/CutoffDates.csv") %>%

  # Remove undesired columns
  dplyr::select(c("CollarID", "DogName", "StartDate", "EndDate")) %>%

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

# We can now check hof often each individual "dispersed"
cut %>%
  group_by(DogName) %>%
  summarize(max(DispersalNo))

# Store the merged cutoff dates
write.csv(cut, "03_Data/02_CleanData/00_General_CutoffDates_POPECOL.csv")

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

# Order the data nicely
data <- data %>% arrange(DogName, Timestamp)

# Note: For Abrahms individuals there is an overlap between data collected by
# Abrahms and data collected by the Staff in Botswana (its actually the same
# data). Because of a tiny mismatch in the timestamps they are not recognized as
# duplicates. However, the coordinates align perfectly and the temporal mismatch
# is minor. Thus, the issue will be resolved once we resample the data to a
# coarser resolution.

################################################################################
#### Visualize Dispersal Phases
################################################################################
# We now want to create a plot illustrating the dispersal phases. For this, we
# need to create a table indicating the start and end of each "phase"
phases <- data %>%
  group_by(DogName) %>%
  nest() %>%
  mutate(data = map(data, function(x){
    phase <- x$State != lag(x$State)
    phase <- replace_na(phase, F)
    x$Phase <- cumsum(phase)
    return(x)
  })) %>%
  unnest(cols = data) %>%
  group_by(DogName, State, Phase) %>%
  summarize(FirstDate = min(Timestamp), LastDate = max(Timestamp)) %>%
  mutate(DogName = as.factor(DogName))

# Reorder the factors
# phases$DogName <- fct_reorder(phases$DogName, phases$DogName, .desc = T)

# Visualize them
ggplot(phases, aes(x = rbind(FirstDate, LastDate), y = DogName)) +

  # Add segments for the resident phase
  geom_segment(data = subset(phases, State == "Resident"), aes(
      x     = FirstDate
    , xend  = LastDate
    , y     = DogName
    , yend  = DogName
  ), size = 2.5, color = "cornflowerblue") +

  # Add segments for the resident phase
  geom_segment(data = subset(phases, State == "Disperser"), aes(
      x     = FirstDate
    , xend  = LastDate
    , y     = DogName
    , yend  = DogName
  ), size = 1.0, color = "black") +

  # Reverse the y scale
  scale_y_discrete(limits = rev(levels(phases$DogName))) +

  # Put a useful title and axis labels
  ggtitle("GPS Observations") +
  xlab("Date") +
  ylab("Name")

# ################################################################################
# #### Adding Pack Information
# ################################################################################
# # Load pack affiliation table
# packs <- read_csv("03_Data/01_RawData/POPECOL/PackAffiliations.csv")
# packs <- subset(packs, DogName %in% toupper(unique(data$DogName)))
# write_csv(packs, "test.csv")
#
# # Coerce the dates to posixct (subtract 2 hours)
# packs$FirstDate <- update(as.POSIXct(packs$FirstDate), hours = 0, tz = "UTC") -
#   hours(2)
# packs$LastDate <- update(as.POSIXct(packs$LastDate), hours = 24, tz = "UTC") -
#   hours(2)
#
# # Let's check if all individuals are represented in the table
# dogs <- unique(data$DogName)
# present <- toupper(dogs) %in% packs$DogName
#
# # Apparently we're missing one individual
# table(present)
# dogs[!present]
#
# # Let's now join the pack information to each GPS fix
# data$Pack <- pbmclapply(
#     X                   = 1:nrow(data)
#   , mc.cores            = detectCores() - 1
#   , ignore.interactive  = T
#   , FUN                 = function(x){
#     pack <- subset(packs, DogName == toupper(data$DogName[x]))
#     index <- data$Timestamp[x] >= pack$FirstDate & data$Timestamp[x] <= pack$LastDate
#     index <- which(index)
#     pack <- pack$Pack[index]
#     pack <- ifelse(is.null(pack), NA, pack)
#     return(pack)
# }) %>% do.call(c, .)
#
# # Identify the different phases
# phases <- data %>%
#   group_by(DogName, State, Pack) %>%
#   summarize(FirstDate = min(Timestamp), LastDate = max(Timestamp)) %>%
#   mutate(DogName = as.factor(DogName))
#
# # Visualize them
# p <- ggplot(phases, aes(x = rbind(FirstDate, LastDate), y = DogName, col = Pack)) +
#
#   # Add segments for the resident phase
#   geom_segment(data = subset(phases, State == "Resident"), aes(
#       x     = FirstDate
#     , xend  = LastDate
#     , y     = DogName
#     , yend  = DogName
#   ), size = 2.5) +
#
#   # Add segments for the resident phase
#   geom_segment(data = subset(phases, State == "Disperser"), aes(
#       x     = FirstDate
#     , xend  = LastDate
#     , y     = DogName
#     , yend  = DogName
#   ), size = 1.0, color = "black") +
#
#   # Reverse the y scale
#   scale_y_discrete(limits = rev(levels(phases$DogName))) +
#
#   # Put a useful title and axis labels
#   ggtitle("GPS Observations") +
#   xlab("Date") +
#   ylab("Name")

# ################################################################################
# #### Adding Pack Information
# ################################################################################
# # Let's load the csv containing all pack affiliations
# packs <- read_csv("03_Data/01_RawData/POPECOL/PackAffiliation.csv")
#
# # Create a column that indicates if a fix was recorded before or after dispersal
# data <- data %>%
#   subset(State == "Disperser") %>%
#   group_by(DogName) %>%
#   summarize(StartDate = min(Timestamp), EndDate = max(Timestamp)) %>%
#   left_join(data, .)
# data$AfterDispersal <- ifelse(data$Timestamp > data$EndDate, T, F)
# data$AfterDispersal[is.na(data$AfterDispersal)] <- F
#
# # Lets join the birth and postdispersal pack to each fix
# data <- left_join(data, packs, by = "DogName")
#
# # During dispersal the individuals are not with their pack anymore so let's put
# # NAs for the pack during dispersal
# data$BirthPack[data$State == "Disperser"] <- NA
#
# # Finally we merge the "BirthPack" and "PostDispersalPack" columns and create
# # a column that indicates the current pack
# data$CurrentPack <- data$BirthPack
# data$CurrentPack[data$AfterDispersal] <-
#   data$PostDispersalPack[data$AfterDispersal]
#
# # Remove any undesired column
# data <- dplyr::select(
#     data
#   , -c(StartDate, EndDate, BirthPack, PostDispersalPack)
# )
#
# # Now there are some manual adjustments that I do according to Dominik's
# # instructions. For Fatalii the pack is unknown before the 18.02.2016
# data$CurrentPack[data$DogName == "Fatalii" &
#   data$Timestamp < "2016-02-17 22:00:00"] <- NA
#
# # Mirage and Kalahri belong to BN only after the 20.09.2017. Let's make anything
# # after dispersal NA. Then add the BN pack after the 20-09-2017
# data$CurrentPack[data$AfterDispersal & data$DogName == "Kalahari"] <- NA
# data$CurrentPack[data$AfterDispersal & data$DogName == "Kalahari" &
#   data$Timestamp > "2017-09-19 22:00:00"] <- "BN"
# data$CurrentPack[data$AfterDispersal & data$DogName == "Mirage"] <- NA
# data$CurrentPack[data$AfterDispersal & data$DogName == "Mirage" &
#   data$Timestamp > "2017-09-19 22:00:00"] <- "BN"
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
# data <- dplyr::select(data, -c("AfterDispersal"))
#
# # Order the dataframe according to: DogName, Timestamp
# data <- arrange(data, DogName, Timestamp)
#
# ################################################################################
# #### Plotting GPS Observation Phases
# ################################################################################
# # Identify all individuals that eventually dispersed
# dispersers <- data %>%
#   subset(State == "Disperser") %>%
#   .[["DogName"]] %>%
#   unique()
#
# # Replace NAs in the current packs after dispersal with character NA (i.e.
# # "NA"). This will avoid some errors.
# data$CurrentPack[is.na(data$CurrentPack) & data$State == "Resident"] <- "NA"
#
# # We want to plot a graph that depicts the period for which each individual was
# # collared and the time during which the individual was dispersing. The colours
# # should indicate to which pack each individual belongs. As a first step we get
# # the earliest and latest gps recording for each dog and collar, and current
# # pack
# collars <- data %>%
#
#   # Subset to dispersers only
#   subset(DogName %in% dispersers) %>%
#
#   # Group by CollarID, DogName and Currentpack
#   group_by(., CollarID, DogName, CurrentPack) %>%
#
#   # Now we can calculate the first and last observation for each group.
#   summarize(.
#     , FirstDate = min(Timestamp)
#     , LastDate  = max(Timestamp)
#   ) %>%
#
#   # We get rid of the "non-character" NAs since they refer to the either
#   # Abrahms' individuals, or the dispersal part of our individuals. We don't
#   # want these in our plot anyways so let's remove them
#   subset(., !is.na(CurrentPack)) %>%
#
#   # Order the data. We put the focus on the pack identities since we assume that
#   # members of the same pack are at the same location. The individual is
#   # therefore not really interesting for us.
#   arrange(., CurrentPack, FirstDate)
#
# # Make the dog names factorial. Note that the factors are now ordered as we specified them in the ordering above. This is helpful for plotting
# collars$DogName <- as.factor(as.character(collars$DogName))
#
# # We could now make the character "NAs" true NAs again
# collars$CurrentPack[collars$CurrentPack == "NA"] <- NA
#
# # Since we also want to depict the dispersal phase, we need to get the dispersal
# # dates. Since these dates are stored in the cutoff dates table we can load it
# # again
# cut <- read_csv("03_Data/02_CleanData/00_General_CutoffDates_POPECOL.csv")
#
# # Prepare the plot
# ggplot(collars, aes(x = rbind(FirstDate, LastDate), y = DogName)) +
#
#   # Add segments for the gps observation phase
#   geom_segment(data = collars, aes(
#       x       = FirstDate
#     , xend    = LastDate
#     , y       = DogName
#     , yend    = DogName
#     , colour  = CurrentPack
#   ), size = 3.5) +
#
#   scale_fill_viridis() +
#
#   # Add segments for the dispersal phase
#   geom_segment(data = cut, aes(
#       x       = StartDate
#     , xend    = EndDate
#     , y       = DogName
#     , yend    = DogName
#   ), size = 2, colour = "black") +
#
#   # Revert the y scale (top to bottom)
#   scale_y_discrete(limits = rev(levels(collars$DogName))) +
#
#   # Put a useful title and axis labels
#   ggtitle("GPS Observations") +
#   xlab("Date") +
#   ylab("Name") +
#
#   # Use the viridis colour scheme for the entire plot
#   scale_color_viridis(
#       discrete  = TRUE
#     , begin     = 0.3
#     , na.value  = "darkgrey"
#   )

################################################################################
#### Resampling Data
################################################################################
# Remove NA fixes
table(is.na(data$x))
table(is.na(data$y))
data <- subset(data, !is.na(x) & !is.na(y))

# Sort data
data <- arrange(data, DogName, Timestamp)

# Some of the data has a very high resolution that we don't need. We will
# therefore subsample to a resolution we can work with.
data <- data %>% group_by(DogName) %>% nest()
data$data <- suppressMessages(
  pbmclapply(data$data
    , ignore.interactive =  T
    , mc.cores = detectCores() - 1
    , FUN = function(x){resFix2(x, hours = 1, start = 1)}
  )
)
data <- unnest(data)

################################################################################
#### Store the Output (as csv and shapefile)
################################################################################
# Write the data to file
write_csv(data, "03_Data/02_CleanData/00_General_Dispersers_POPECOL.csv")

# Create SpatialLinesDataFrame for visualization
data <- data %>% group_by(DogName) %>% nest()
data$Tracks <- suppressMessages(
  pbmclapply(1:nrow(data)
    , ignore.interactive = T
    , mc.cores = detectCores() - 1
    , function(x){
      coords <- data$data[[x]]
      coords$DogName <- data$DogName[x]
      x <- coords
      coordinates(x) <- c("x", "y")
      crs(x) <- CRS("+init=epsg:4326")
      lines <- spLines(x)
      lines <- createSegments(lines)
      lines <- as(lines, "SpatialLinesDataFrame")
      lines@data <- x@data[1:(nrow(x) - 1), ]
      crs(lines) <- CRS("+init=epsg:4326")
      return(lines)
    }
  )
)

# Create SpatialPointsDataFrame for visualization
data$Points <- lapply(data$data, function(x){
  coordinates(x) <- c("x", "y")
  crs(x) <- CRS("+init=epsg:4326")
  return(x)
})

# Store shapefile
tracks <- do.call(rbind, data$Tracks)
writeOGR(tracks
  , dsn       = "03_Data/02_CleanData"
  , layer     = "00_General_Dispersers_POPECOL"
  , driver    = "ESRI Shapefile"
  , overwrite = TRUE
)

################################################################################
#### Create KML Files
################################################################################
# We may also want to create kml files for visualization in google earth. Let's
# first point to the shapes visualized in google earth.
shape1 <- "http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png"
shape2 <- "http://maps.google.com/mapfiles/kml/shapes/square.png"

# Prepare a folder into which we will store the kmls
dir.create("03_Data/02_CleanData/00_KML")

# Loop through the individuals and prepare kml files for them
lapply(1:nrow(data), function(x){

  # Access required data
  name    <- data$DogName[[x]]
  points  <- data$Points[[x]]
  track   <- data$Tracks[[x]]

  # Identify start and endpoints
  first <- as(points[1, ], "SpatialPoints")
  last  <- as(points[nrow(points), ], "SpatialPoints")

  # Make row-names valid
  row.names(points) <- as.character(1:nrow(points))
  row.names(track) <- as.character(1:nrow(track))

  # Create spacetime object from points
  points_st <-STIDF(
      sp   = as(points, "SpatialPoints")
    , time = points$Timestamp
    , data = points@data
  )

  # Create spacetime object from track
  track_st <- STIDF(
      sp   = as(track, "SpatialLines")
    , time = track$Timestamp
    , data = track@data
  )

  # Generate a name for the kml file
  filename <- paste0("03_Data/02_CleanData/00_KML/", name, ".kml")

  # Generate kml file (note that for some individuals we can't produce an
  # info-table because there are too many GPS fixes)
  kml_open(filename)
  kml_layer(first
    , colour     = "green"
    , size       = 1
    , shape      = shape2
    , LabelScale = 0
  )
  kml_layer(last
    , colour     = "red"
    , size       = 1
    , shape      = shape2
    , LabelScale = 0
  )
  kml_layer(track_st
    , colour       = State
    , colour_scale = rev(viridis(2, begin = 0.6))
    , width        = 2.5
  )
  tryCatch(kml_layer(points_st
    , colour       = State
    , colour_scale = rev(viridis(2, begin = 0.6))
    , size         = 0.75
    , balloon      = T
    , shape        = shape
    , LabelScale   = 0
  ), error = function(e){
    kml_layer(points_st
      , colour       = State
      , colour_scale = rev(viridis(2, begin = 0.6))
      , size         = 0.75
      , balloon      = F
      , shape        = shape
      , LabelScale   = 0
    )
  })
  kml_close(filename)
})