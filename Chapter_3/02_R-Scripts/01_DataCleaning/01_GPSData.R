################################################################################
#### Preparation of GPS Data
################################################################################
# Description: Use our custom "wilddogr" package to download and clean all GPS
# data of wilddogs that is available on the dropbox

# Clean environment
rm(list = ls())

# Load required packages
library(tidyverse)  # For data wrangling
library(wilddogr)   # To download GPS data from dropbox
library(raster)     # For handling spatial data
library(rgdal)      # For loading and storing spatial data
library(sf)         # For plotting spatial data
library(ggpubr)     # To arrange multiple plots

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

################################################################################
#### Download Data From Dropbox and Combine with Data From Briana
################################################################################
# Identify files on dropbox (ignore warning about expected 3 pieces)
cat("Checking for GPS data on dropbox...\n")
rdrop2::drop_auth(cache = F)
files <- dog_files(rvc = F)

# Let's take a look at the downloadable data
head(files)

# Download the data
if (!file.exists("03_Data/01_RawData/POPECOL/Cleaned_GPSData.csv")) {
    downloaded <- dog_download(
        x         = files
      , clean     = T
      , overwrite = T
      , outdir    = "03_Data/01_RawData/POPECOL"
      , legacy    = T
    )
  } else {
    cat("GPS data has already been downloaded and will not be downloaded again...\n")
}

# Load the downloaded data
cat("Merging POPECOL's GPS data with Abrahm's GPS data...\n")
dat1 <- read_csv("03_Data/01_RawData/POPECOL/Cleaned_GPSData.csv", show_col_types = F)

# Find all dogs that eventually dispersed
dispersers <- dat1 %>%
  subset(State == "Disperser") %>%
  pull(DogName) %>%
  unique()

# Subset to data of those individuals only
dat1 <- subset(dat1, DogName %in% dispersers) %>%
  mutate(Source = "Popecol") %>%
  rename(ID = DogName)

# Let's also add the individuals that were observed by Briana
dat2 <- read_csv("03_Data/01_RawData/ABRAHMS/DispersalPaths.csv", show_col_types = F) %>%
  filter(., !is.na(`Longitude..deg.`)) %>%
  dplyr::select(.,
      ID        = `id`
    , x         = `Longitude..deg.`
    , y         = `Latitude..deg.`
    , Timestamp = `Timestamp`
  ) %>%
  mutate(., CollarID = NA, DOP = NA, State = "Disperser", Sex = "M", Source = "Abrahms")

# Put data together and clean up
dat <- rbind(dat1, dat2)
rm(dat1, dat2, dispersers)

# Remove data with NAs in coordinates
# dat <- subset(dat, !is.na(x) & !is.na(y))

# Some visualizations
p1 <- ggplot(subset(dat, !is.na(x)), aes(x = x, y = y)) +
  geom_point(size = 0.1) +
  theme_minimal() +
  coord_sf()
p2 <- ggplot(subset(dat, State == "Disperser" & !is.na(x)), aes(x = x, y = y, col = as.factor(ID))) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_path(linewidth = 0.2, alpha = 0.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_sf()

# Arrange plots
ggarrange(p1, p2, ncol = 1)

################################################################################
#### Merge Individuals that Dispersed Together
################################################################################
# Francisco and Saturday moved together for their first exploratory movements
ggplot(subset(dat, State == "Disperser" & !is.na(x) & ID %in% c("Francisco", "Saturday")), aes(x = x, y = y, col = as.factor(ID))) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_path(linewidth = 0.2, alpha = 0.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_sf()

# Let's combine their data until they split, which is on the 14.06.2021 (or
# slightly after, but the first dispersal phase ends on that day so we can use
# it to merge the data)
indices <- dat$ID == "Francisco" & dat$State == "Disperser" & dat$Timestamp <= "2021-06-14"
# dat[indices, ]

# Override with Saturdays information. This will allow us to utilize Francisco's
# fixes to augment the track by saturday for this period. Note that this will
# introduce "fake duplicate" fixes for Saturday, in case both Saturday and
# Francisco fixed at the same time. However, we will later rarify the data.
dat$ID[indices]       <- "Saturday"
dat$CollarID[indices] <- "45228"

# Arrange the data again
dat <- arrange(dat, ID, desc(Timestamp))

# Check again
ggplot(subset(dat, State == "Disperser" & !is.na(x) & ID %in% c("Francisco", "Saturday")), aes(x = x, y = y, col = as.factor(ID))) +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_path(linewidth = 0.2, alpha = 0.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  coord_sf()

# Store the data to file
write_csv(dat, "03_Data/02_CleanData/Dispersers.csv")

# Clear graphs
graphics.off()

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/01_GPSData.rds")
cat("Done :)\n")
