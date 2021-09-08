################################################################################
#### Preparation of GPS Data
################################################################################
# Description: I'll use the custom package "wilddogr" to automatically download
# and clean all GPS data that is available to us

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

# # Identify files on dropbox
# files <- dog_files(rvc = F)
#
# # Let's take a look at the downloadable data
# head(files)
#
# # Download the data
# downloaded <- dog_download(
#     x         = files
#   , clean     = T
#   , overwrite = T
#   , outdir    = "03_Data/01_RawData/POPECOL"
# )
#
# # Let's take a look at the data
# dat <- read_csv(downloaded)
# nrow(dat)
# dat2 <- read_csv("03_Data/01_RawData/POPECOL/Dispersers.csv")
# dispersers2 <- unique(subset(dat2, State == "Disperser")$DogName)
#
# # Keep only individuals that eventually dispersed
# dispersers <- unique(subset(dat, State == "Disperser")$DogName)
# length(dispersers)
# length(dispersers2)
# dispersers2[!dispersers2 %in% dispersers]
#
# # Check how many fixes during dispersal there are
# table(dat$State)
# table(dat$DogName, dat$State)
#
# # Count individuals
# length(unique(dat$DogName))

################################################################################
#### NEED TO ADD BRIANAS INDIVIDUALS!!!
################################################################################

# Load cleaned data
dat <- read_csv("03_Data/01_RawData/POPECOL/Dispersers.csv")

# Keep only individuals that eventually dispersed
dispersers <- unique(dat$DogName[dat$State == "Disperser"])
dat <- subset(dat, DogName %in% dispersers)

# How many individuals are there?
length(unique(dat$DogName))
table(dat$DogName, dat$State)

################################################################################
#### ASSIGN SEASON TO EACH FIX!!!
################################################################################

# Plot the data
ggplot(dat, aes(x = x, y = y, col = factor(DogName))) +
  geom_point() +
  geom_path(size = 0.5) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_sf()
