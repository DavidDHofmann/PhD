############################################################
#### Preparing the Social Landscape
############################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(data.table)
library(viridis)
library(adehabitatHR)
library(raster)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(rasterVis)
library(rgdal)
library(rgeos)
library(parallel)

# Load custom functions
source("Functions.r")

# Make use of multicore abilities
beginCluster()

# Reload data from Part I
dat <- readRDS("99_SocialLandscape(PartI).rds")

############################################################
#### Process UDs
############################################################
# Stack the uds that belong to the same period
ud_final <- dat %>%
  group_by(Group) %>%
  nest()

# Extract the Groups for later
names <- ud_final$Group

# Stack the nested raster together and do some calculations
ud_final <- ud_final %>%

  # Create new column
  mutate(ud = mclapply(data, mc.cores = (detectCores() - 1), function(x){

    # Stack them
    x <- stack(x$ud)

    # Calculate probability of not encountering pack i
    values(x) <- 1 - values(x) / sum(values(x))

    # Calculate probability of encounter at all
    x <- 1 - calc(x, prod)

    # Store the result to a temporary file to save memory
    x <- writeRaster(x, tempfile())

    # Collect garbage
    gc()

    # Return the stack
    return(x)
  }))

# Stack all rasters into one stack
ud_final <- stack(ud_final$ud)

# Assign some nice layer names
names(ud_final) <- names

# Write the stack to file
writeRaster(
    ud_final
  , filename  = "03_Data/02_CleanData/05_SocialFeatures_SocialLandscape(UtilisationDistributions)"
  , format    = "raster"
  , overwrite = TRUE
  , options   = c("INTERLEAVE = BAND", "COMPRESS = LZW")
)

############################################################
#### Process HRs
############################################################
# Put all home ranges into one SpatialPolygonsDataFrame
hr_final <- dat$hr %>% do.call(rbind, .)

# Put the group and date into the dataframe of the polygons
hr_final@data <- data.frame(
    Date = dat$Group %>% as.POSIXct(., tz = "UTC")
  , Pack = dat$CurrentPack
)

# Store HRs to file
writeOGR(
    hr_final
  , "03_Data/02_CleanData"
  , "05_SocialFeatures_SocialLandscape(HomeRanges)"
  , driver = "ESRI Shapefile"
  , overwrite = TRUE
)

# Save the data to file for part III
saveRDS(hr_final, "99_SocialLandscape(PartII).rds")

# Remove the rds object from part I
file.remove("99_SocialLandscape(PartI).rds")

# Terminate the cluster
endCluster()
