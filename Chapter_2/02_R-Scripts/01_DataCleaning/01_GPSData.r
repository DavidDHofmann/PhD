################################################################################
#### Preparation of GPS Data
################################################################################
# Description: Use our custom "wilddogr" package to download and clean all GPS
# data of wilddogs that is available on the dropbox

# Clean environment
rm(list = ls())

# Load required packages
library(wilddogr)   # To download GPS data from dropbox
library(raster)     # For handling spatial data
library(rgdal)      # For loading and storing spatial data
library(sf)         # For plotting spatial data

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

################################################################################
#### Download Data From Dropbox and Combine with Data From Briana
################################################################################
# Identify files on dropbox (ignore warning about expected 3 pieces)
files <- dog_files(rvc = F)

# Let's take a look at the downloadable data
head(files)

# Download the data
# downloaded <- dog_download(
#     x         = files
#   , clean     = T
#   , overwrite = T
#   , outdir    = "03_Data/01_RawData/POPECOL"
# )

# Load the downloaded data
dat1 <- read_csv("03_Data/01_RawData/POPECOL/Cleaned_GPSData.csv")

# Find all dogs that eventually dispersed
dispersers <- dat1 %>%
  subset(State == "Disperser") %>%
  pull(DogName) %>%
  unique()

# Subset to data of those individuals only
dat1 <- subset(dat1, DogName %in% dispersers) %>%
  mutate(Source = "Popecol")

# Let's also add the individuals that were observed by Briana
dat2 <- read_csv("03_Data/01_RawData/ABRAHMS/DispersalPaths.csv") %>%
  filter(., !is.na(`Longitude..deg.`)) %>%
  dplyr::select(.,
      DogName   = `id`
    , x         = `Longitude..deg.`
    , y         = `Latitude..deg.`
    , Timestamp = `Timestamp`
  ) %>%
  mutate(., CollarID = NA, DOP = NA, State = "Disperser", Sex = "M", Source = "Abrahms")

# Put data together and clean up
dat <- rbind(dat1, dat2)
rm(dat1, dat2, dispersers)

# Remove data with NAs in coordinates
dat <- subset(dat, !is.na(x) & !is.na(y))

# Some visualizations
ggplot(dat, aes(x = x, y = y)) +
  geom_point(size = 0.1) +
  theme_minimal() +
  coord_sf()
ggplot(subset(dat, State == "Disperser"), aes(x = x, y = y, col = as.factor(DogName))) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_path(size = 0.2, alpha = 0.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_sf()

# Store the data to file
write_csv(dat, "03_Data/02_CleanData/00_General_Dispersers.csv")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/01_GPSData_SessionInfo.rds")
