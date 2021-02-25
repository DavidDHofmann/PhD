################################################################################
#### Preparation of Data for Dryad
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load packages
library(tidyverse)

# Load GPS data of dispersers
data <- read_csv("03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv")
data <- subset(data, State == "Disperser")

# Generate random numbers using a seed
set.seed(12345)
data$identifier <- sample(1:nrow(data), replace = F)

# Sort data by identifier
data <- arrange(data, identifier)
names(data)[names(data) == c("x", "y")] <- c("lon", "lat")

# Prepare the same dataframe with fewer columns
data_dryad <- select(data, identifier, lon, lat)

# Round data
data_dryad$lon <- round(data_dryad$lon, 2)
data_dryad$lat <- round(data_dryad$lat, 2)

# Store dataframes
write_csv(data_dryad, "03_Data/03_Results/GPS_Dryad.csv")
write_csv(data, "03_Data/03_Results/GPS_Full.csv")
