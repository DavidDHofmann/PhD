################################################################################
####
################################################################################
# Description: Use floodmaps to mask out anything that is covered by water in
# the vegetation maps. Specifically, we will for each floodmap create a
# corresponding vegetation map where water has been masked out. This requires us
# to match the date of the floodmaps with the date of the vegetation maps. Note:
# vegetation maps are based on data from day 65 of one year to day 64 of the
# next year -> From the MODIS documentation
# (https://lpdaac.usgs.gov/documents/112/MOD44B_User_Guide_V6.pdf): The product
# date (of the vegetation data) refers to the start date of the annual period so
# a product with ID “2005065” was produced with data from 2005065 – 2006064. The
# start date of all MODIS VCF products is yyyy065 (where yyyy refers to the 4
# digit year). This originally derived from the first full 16-day composite
# period in the MODIS data record which begins with day of year 2000065. Hence,
# the layername of each MODIS vegetation map refers to the period that starts on
# day 65 of that year!

# Clear R's brain
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(terra)       # To handle raster data
library(lubridate)   # To handle dates
library(tidyverse)   # To wrangle data

# Let's prepare tibbles for the floodmaps and the vegetation maps so that we can
# more easily subset the data. First for the flood data
flood_dates <- dir(
      path       = "03_Data/02_CleanData/00_Floodmaps/02_Resampled"
    , pattern    = ".tif$"
  ) %>% ymd()

# Now for the shrubs/grassland cover
grass_dates <- "03_Data/02_CleanData/00_Vegmaps/ShrubCover.grd" %>%
  rast() %>%
  names() %>%
  tibble(
      Year     = .
    , FromDate = as.Date(64, origin = paste0(., "-01-01"))
    , ToDate = as.Date(63, origin = paste0(., "-01-01")) + years(1)
  )

# Same for tree cover
tree_dates <- "03_Data/02_CleanData/00_Vegmaps/TreeCover.grd" %>%
  rast() %>%
  names() %>%
  tibble(
      Year     = .
    , FromDate = as.Date(64, origin = paste0(., "-01-01"))
    , ToDate = as.Date(63, origin = paste0(., "-01-01")) + years(1)
  )

# Get dispersal dates too
disp_dates <- "03_Data/02_CleanData/00_General_Dispersers.csv" %>%
  read_csv() %>%
  pull(Timestamp) %>%
  as.Date() %>%
  unique()

# Find closest floodmap for each dispersal date
closest1 <- as.Date(NA)
closest2 <- as.Date(NA)
for (i in 1:length(disp_dates)){
  closest1[i] <- flood_dates[which(abs(disp_dates[i] - flood_dates) ==
    min(abs(disp_dates[i] - flood_dates)))][1]
  closest2[i] <- flood_dates[which(abs(disp_dates[i] - flood_dates) ==
    min(abs(disp_dates[i] - flood_dates)))][2]
}

# Put the dates together
dates <- data.frame(
    Dispersal   = disp_dates
  , Closest1    = closest1
  , Closest2    = closest2
  , Difference  = abs(disp_dates - closest1)
)
arrange(dates, -Difference)

# Use floodmaps to mask water in the vegetation maps of that year
