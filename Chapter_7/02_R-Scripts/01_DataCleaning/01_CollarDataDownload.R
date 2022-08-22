################################################################################
#### Download Collar Data
################################################################################
# Download GPS and Activity data stored in the "KML-Files"
# folder on dropbox

# Clear R's brain
rm(list = ls())

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_7")

# Load required packages
library(rdrop2)      # To get access to the dropbox through R
library(tools)       # To handle filenames and extensions
library(tidyverse)   # For data wrangling
library(readxl)      # To read excel files
library(lubridate)   # To handle dates

# May need to allow rdrop2 to access the dropbox
drop_auth(new_user = T)

################################################################################
#### Dropbox File Preparation
################################################################################
# Check if the KML-Files folder is in the top directory on dropbox. If not, get
# the correct path to where it is stored.
cat("Retrieving information on GPS and Activity files from dropbox...\n")
if (drop_exists("KML-Files")) {
    path <- "KML-Files"
  } else {
    path <- drop_dir(recursive = T)
    path <- path$path_display[path$name == "KML-Files"]
}

# Identify all files in the "KML-Files" folder
files <- path %>%
  drop_dir(recursive = T) %>%
  filter(.tag == "file") %>%
  mutate(filetype = file_ext(name))

# Check them
head(files, n = 20)

# Subset to GPS data
gps <- files %>%
  subset(filetype == "csv" & grepl(name, pattern = "GPS")) %>%
  subset(grepl(path_lower, pattern = "inactive")) %>%
  mutate(Type = "GPS")

# Subset to activity data
act <- files %>%
  subset(filetype == "csv" & (
        grepl(name, pattern = "ACT")
      | grepl(name, pattern = "ACTIVITY")
    )
  ) %>%
  subset(grepl(path_lower, pattern = "inactive")) %>%
  subset(!grepl(path_display, pattern = "Gabs")) %>%
  mutate(Type = "ACT")

# Put GPS and Activity file into one dataframe
dat <- rbind(gps, act)

# Keep only the columns of interest
dat <- select(dat, c(name, path_display, Type))

# Extract DogID, CollarID, and the Timestamp at which the file was generated
dat <- strsplit(dat$path_display, split = "/") %>%
  lapply(function(x) {c(x[[4]], x[[5]])}) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  setNames(c("String1", "String2")) %>%
  mutate(
      DogID     = str_match(String1, "[0-9]\\_\\s*(.*?)\\s*_collar")[, 2]
    , CollarID  = str_match(String1, "collar#\\s*(.*)")[, 2]
    , Timestamp = str_match(String2, "[0-9]\\_\\s*(.*?)\\s*.csv")[, 2]
  ) %>%
  select(-c(String1, String2)) %>%
  cbind(dat, .) %>%
  arrange(Type, DogID, CollarID)

# In some cases there are duplicate entries (indicated by a "_b" in the
# timestamp). Remove them.
dat <- subset(dat, !grepl(Timestamp, pattern = "_b"))

# One of the collars also appears to be named incorrectly. Correct it.
dat <- mutate(dat, CollarID = substr(CollarID, start = 1, stop = 5))

# Make nice filenames for download
dat <- mutate(dat, Filename = paste0(
  Type, "_", DogID, "_", CollarID, "_", Timestamp, ".csv"
))

# Depending on the type of data, they'll go into a different folder
dat <- mutate(dat, Filepath = file.path(
  getwd(), "03_Data/01_RawData/POPECOL", ifelse(Type == "GPS", "GPS", "ACT"), Filename
))

# Create the necessary folders
dir.create(file.path(getwd(), "03_Data/01_RawData/POPECOL"), showWarnings = F)
dir.create(file.path(getwd(), "03_Data/01_RawData/POPECOL/GPS"), showWarnings = F)
dir.create(file.path(getwd(), "03_Data/01_RawData/POPECOL/ACT"), showWarnings = F)

# Take a look
head(dat, n = 50)

# How many GPS and Activity files are there?
count(dat, Type)

# Ensure that there are no duplicates
table(duplicated(dat$Filepath))

################################################################################
#### Download
################################################################################
# Let's go through the files and download them
cat("Downloading GPS and Activity data from Dropbox...\n")
pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
sapply(1:nrow(dat), function(i) {
  drop_download(
      path       = dat$path_display[i]
    , local_path = dat$Filepath[i]
    , overwrite  = F
  )
  setTxtProgressBar(pb, i)
  return(dat$Filepath[i])
})

################################################################################
#### Collar Handling and Dispersal Dates
################################################################################
# We also need to download two helper files. (1) The file that contains the
# collar handling dates, and then (2) also the file that contains dates during
# which individuals were dispersing
helpers <- subset(files
  , grepl(name, pattern = "overview dispersal dates.xlsx") |
    grepl(name, pattern = "Collar Settings.xlsx")
)

# Remove undesired columns
helpers <- select(helpers, c(name, path_display))

# Make sure there are exactly two files
if (nrow(helpers) != 2) {
  stop("Something is wrong with the 'collar handling' or 'dispersal dates'
  file...\n")
}

# Give them a nice name and specify the filepath to which we will store them
helpers <- mutate(helpers
  , FilenameTemp  = c("DispersalDates.xlsx", "CollarHandling.xlsx")
  , FilepathTemp  = file.path(tempdir(), FilenameTemp)
  , FilenameFinal = c("DispersalDates.csv", "CollarHandling.csv")
  , FilepathFinal = file.path(getwd(), "03_Data/01_RawData/POPECOL", FilenameFinal)
)

# Let's download them to the temporary directory (I want to directly convert
# them into a .csv after download)
cat("Downloading helper files from dropbox...\n")
pb <- txtProgressBar(min = 0, max = nrow(dat), style = 3)
sapply(1:nrow(helpers), function(i) {
  drop_download(
      path       = helpers$path_display[i]
    , local_path = helpers$FilepathTemp[i]
    , overwrite  = T
  )
  setTxtProgressBar(pb, i)
  return(dat$FilepathTemp[i])
})

# Load files and remove undesired columns, then store as .csv
helpers$FilepathTemp[1] %>%
  read_excel() %>%
  select(DogID = DogName, CollarID, StartDate_UTC, EndDate_UTC, DispersalNo) %>%
  write_csv(., file = helpers$FilepathFinal[1])
helpers$FilepathTemp[2] %>%
  read_excel(skip = 1) %>%
  subset(!is.na(`Collar Nr.`) & !is.na(`Dog Name`)) %>%
  dplyr::select(
      DogID     = `Dog Name`
    , CollarID  = `Collar Nr.`
    , DogCode   = `Dog Code`
    , Sex       = `Sex`
    , FirstDate = `Collaring.Date.Text`
    , LastDate  = `Stop.Recording.Date.Text`
  ) %>%
  mutate(
      StartDate_UTC = ymd_hms(FirstDate) - hours(2) + minutes(5) # Subtract 2 hours to get utc time, add 5 mins for tolerance
    , EndDate_UTC   = ymd_hms(LastDate) - hours(2) + minutes(5) # Subtract 2 hours to get utc time, add 5 mins for tolerance
  ) %>%
  select(-c(FirstDate, LastDate)) %>%
  arrange(DogID, CollarID) %>%
  write_csv(., file = helpers$FilepathFinal[2])
