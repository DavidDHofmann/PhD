################################################################################
#### Setup
################################################################################
# Description: Set up folder structure for this project

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_7")

# Let's prepare the main folders that we need
dir.create("01_General", showWarnings = F)
dir.create("02_R-Scripts", showWarnings = F)
dir.create("03_Data", showWarnings = F)
dir.create("03_Data/01_RawData", showWarnings = F)
dir.create("03_Data/02_CleanData", showWarnings = F)
dir.create("03_Data/03_Results", showWarnings = F)
dir.create("04_Manuscript", showWarnings = F)
