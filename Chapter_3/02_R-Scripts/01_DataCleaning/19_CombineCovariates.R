################################################################################
#### Combine Covariates
################################################################################
# Description: Generating a single file that links all covariate, keeps track of
# their dates, etc.

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)    # For data wrangling
library(terra)        # To handle spatial rasters
library(lubridate)    # To handle dates
library(hms)          # To handle times

# Change the working directory
setwd("/home/david/Schreibtisch")

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

# List of all covariates
covs <- expand_grid(
      Type      = c("Static", "Dynamic", "DynamicAggregated")
    , Covariate = c("Humans", "Forest", "Trees", "Shrubs", "Water", "NDVI", "DistanceToWater", "DistanceToPans", "Temperature", "Precipitation"
  )
) %>% subset(!(Type == "DynamicAggregated" & Covariate == "DistanceToPans"))

# Only run this if the respective file has not been created yet
if (!file.exists("03_Data/02_CleanData/Covariates.rds")) {

  # Go through the files and assign a link to their raster layers if they exist
  covs$Filename <- sapply(seq_len(nrow(covs)), function(i) {
    file1 <- file.path("03_Data/02_CleanData", paste0(covs$Covariate[i], covs$Type[i], ".tif"))
    file2 <- file.path("03_Data/02_CleanData", paste0(covs$Covariate[i], ".tif"))
    if (file.exists(file1)) {
      return(file1)
    } else if (file.exists(file2)) {
      return(file2)
    }
  })

  # Extract layer dates for each covariate
  covs$Dates <- lapply(seq_len(nrow(covs)), function(x) {

    # Get the respective covariate, load it, and extract dates
    cov <- paste0(covs$Covariate[x], covs$Type[x])
    if (cov %in% c("TreesDynamic", "ShrubsDynamic", "NDVIDynamic", "WaterDynamic", "DistanceToWaterDynamic", "TreesDynamicAggregated", "ShrubsDynamicAggregated", "NDVIDynamicAggregated", "WaterDynamicAggregated", "DistanceToWaterDynamicAggregated")) {
        covariate <- rast(paste0("03_Data/02_CleanData/", cov, ".tif"))
        covariate_dates <- covariate %>%
          names() %>%
          substr(start = nchar(.) - 9, stop = nchar(.)) %>%
          ymd(tz = "UTC") %>%
          tibble(Layerdate = .) %>%
          mutate(Layerindex = 1:n())
      } else if (cov %in% c("PrecipitationDynamic", "TemperatureDynamic", "PrecipitationDynamicAggregated", "TemperatureDynamicAggregated")) {
        covariate <- rast(paste0("03_Data/02_CleanData/", cov, ".tif"))
        covariate_dates <- covariate %>%
          names() %>%
          substr(start = nchar(.) - 18, stop = nchar(.)) %>%
          ymd_hms(tz = "UTC") %>%
          tibble(Layerdate = .) %>%
          mutate(Layerindex = 1:n())
      } else if (cov == "DistanceToPansDynamic") {
        covariate <- rast(paste0("03_Data/02_CleanData/", cov, ".tif"))
        covariate_dates <- covariate %>%
          names() %>%
          ym(tz = "UTC") %>%
          tibble(Layerdate = .) %>%
          mutate(Layerindex = 1:n())
      } else if (cov %in% c("PrecipitationStatic", "TemperatureStatic")) {
        covariate <- rast(paste0("03_Data/02_CleanData/", cov, ".tif"))
        covariate_dates <- covariate %>%
          names() %>%
          parse_hms() %>%
          as.character() %>%
          tibble(Layerdate = .) %>%
          mutate(Layerindex = 1:n()) %>%
          expand_grid(., Date = seq(ymd("2000-01-01"), ymd("2022-03-01"), by = "day")) %>%
          mutate(Layerdate = ymd_hms(paste(Date, Layerdate))) %>%
          dplyr::select(-Date)
      } else {
        covariate_dates <- now(tzone = "UTC") %>%
          tibble(Layerdate = .) %>%
          mutate(Layerindex = 1:n())
    }

    # Return the derived dates
    return(covariate_dates)
  })

  # Store to file
  write_rds(covs, "03_Data/02_CleanData/Covariates.rds")

}

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/01_DataCleaning/19_CombineCovariates.rds")
cat("Done :)\n")
