################################################################################
#### Dynamic Vegetation Maps
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

# Load the layers we want to merge
flood <- rast("03_Data/02_CleanData/01_LandCover_WaterCoverDynamic.grd")
shrub <- rast("03_Data/02_CleanData/00_Vegmaps/ShrubCover.grd")
trees <- rast("03_Data/02_CleanData/00_Vegmaps/TreeCover.grd")

# Extract dates from layernames (again, not that vegetation layers correspond to
# a period)
flood_dates <- ymd(names(flood))
veg_dates <- names(shrub) %>%
  tibble(
      Year     = .
    , FromDate = as.Date(64, origin = paste0(., "-01-01"))
    , ToDate   = as.Date(63, origin = paste0(., "-01-01")) + years(1)
  ) %>% mutate(Period = map2(FromDate, ToDate, function(x, y) {
    seq(x, y, by = "day")
  }))

# Let's assess for each floodmap into which period of the vegation they fall
merge <- sapply(flood_dates, function(x) {

  # Check if the floodmap is in one of those periods
  period <- sapply(veg_dates$Period, function(y) {
    x %in% y
  }) %>% which() %>% first()

  # Might fall into a date that is not covered by those periods, so let's
  # subtract days until the date falls into one of the periods
  i <- 1
  while (is.na(period)) {
    period <- sapply(veg_dates$Period, function(y) {
      (x - days(i)) %in% y
    }) %>% which() %>% first()
    i <- i + 1
  }

  # Return the year of the respective period
  return(veg_dates$Year[period])
}) %>% tibble(FloodDate = flood_dates, VegDate = .)

# Loop through the dates and merge the corresponding layers
lapply(1:nrow(merge), function(x) {

  # Get the floodmap and the vegetation layer of that date
  flood_map <- flood[[which(flood_dates == merge$FloodDate[x])]]
  shrub_map <- shrub[[which(veg_dates$Year == merge$VegDate[[x]])]]
  trees_map <- trees[[which(veg_dates$Year == merge$VegDate[[x]])]]

  # Change anything that is covered by water to 0% vegetation
  trees_map <- mask(trees_map, flood_map, maskvalues = 1, updatevalue = 0)
  shrub_map <- mask(shrub_map, flood_map, maskvalues = 1, updatevalue = 0)
  plot(trees_map)
  plot(shrub_map)

})
