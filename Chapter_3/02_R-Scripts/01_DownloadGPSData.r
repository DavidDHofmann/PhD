################################################################################
#### Cleaning and Preparing Tracks from GPS Fixes
################################################################################
# Description: In this script I put together all gps fixes. I am using a custom
# package called "wilddogr" for this. It allows automatic download and cleaning
# of the files stored in the dropbox folder.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(wilddogr) # Custom package to download GPS data from dropbox

# Download all GPS data into a folder
dir.create("03_Data/01_RawData/POPECOL")
todownload <- dog_files(rvc = T)
downloaded <- dog_download(todownload
  , clean     = T
  , outdir    = "03_Data/01_RawData/POPECOL"
  , overwrite = T
)
