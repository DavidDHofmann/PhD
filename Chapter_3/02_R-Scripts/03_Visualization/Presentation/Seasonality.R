################################################################################
#### Seasonality
################################################################################
# Description: Visualization of seasonality patterns in our study area.

# Clear R's brain
rm(list = ls())

# Change the working directory.
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(terra)          # To handle raster files
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(ggpubr)         # To arrange multiple plots

# Load temperature data
temp <- "03_Data/02_CleanData/00_Tempmaps" %>%
  dir(pattern = ".grd$", full.names = T) %>%
  rast()

# Load precipitation data
prec <- "03_Data/02_CleanData/00_Rainmaps" %>%
  dir(pattern = ".grd$", full.names = T) %>%
  rast()

# Load floodmaps
flood <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled" %>%
  dir(pattern = ".tif$", full.names = T) %>%
  rast()

# Load NDVI
ndvi <- "03_Data/02_CleanData/01_LandCover_NDVI.grd" %>%
  rast()

# Identify dates for the temperature data
temp_dates <- names(temp) %>%
  tibble(
      Date     = substr(., start = 1, stop  = 10)
    , Hour     = substr(., start = 12, stop = 13)
    , Datetime = make_datetime(year(Date), month(Date), day(Date), Hour, 0, 0)
  ) %>%
  pull(Datetime)

# Identify dates for the precipitation data
prec_dates <- names(prec) %>%
  tibble(
      Date     = substr(., start = 1, stop  = 10)
    , Hour     = substr(., start = 12, stop = 13)
    , Datetime = make_datetime(year(Date), month(Date), day(Date), Hour, 0, 0)
  ) %>%
  pull(Datetime)

# Identify dates for the floodmaps
flood_dates <- "03_Data/02_CleanData/00_Floodmaps/02_Resampled" %>%
  dir(pattern = ".tif$") %>%
  basename() %>%
  ymd()

# Identify dates for the NDVI maps
ndvi_dates <- names(ndvi) %>% ymd()

# Assign dates as layernames
names(temp) <- temp_dates
names(prec) <- prec_dates
names(flood) <- flood_dates
names(ndvi) <- ndvi_dates

################################################################################
#### Flood
################################################################################
# Identify cover of dryland (0), water (1) and clouds (2) in each floodmap
flood_summary <- freq(flood, bylayer = T) %>%
  as.data.frame() %>%
  mutate(count = count / ncell(flood[[1]])) %>%
  pivot_wider(
    , id_cols     = layer
    , names_from  = value
    , values_from = count
    , values_fill = 0
  ) %>%
  setNames(c("Layer", "Dryland", "Flood", "Cloud")) %>%
  mutate(Date = flood_dates)

# Keep only floodmaps with cloud cover below 5%
flood_summary <- subset(flood_summary, Cloud < 0.05)

# Put data into a nice tibble for plotting
flood_summary <- tibble(
    Date  = flood_summary$Date
  , Flood = flood_summary$Flood
  , Year  = year(Date)
  , Month = month(Date)
  , Week  = week(Date)
  , Day   = yday(Date)
)

# Visualize flood extent by day of the year
p1 <- ggplot(flood_summary, aes(x = Day, y = Flood, col = Year, group = Year)) +
  geom_line() +
  theme_minimal() +
  scale_color_viridis_c() +
  ggtitle("Flood Extent")

################################################################################
#### Temperature
################################################################################
# Compute average temperature in temperature map
temp_summary <- global(temp, mean)

# Put data into a nice dataframe
temp_summary <- tibble(
    Date  = date(temp_dates)
  , Temp  = temp_summary$mean
  , Year  = year(Date)
  , Month = month(Date)
  , Week  = week(Date)
  , Day   = yday(Date)
  , Hour  = hour(temp_dates)
)

# Visualize temperature by day of the year
p2 <- temp_summary %>%
  group_by(Year, Day) %>%
  summarize(Temp = mean(Temp)) %>%
  ggplot(aes(x = Day, y = Temp, col = Year, group = Year)) +
  geom_line() +
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Temperature")

# Visualize temperature by time of the day
p3 <- ggplot(temp_summary, aes(x = Hour, y = Temp, col = Month)) +
  geom_jitter(size = 0.1) +
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Temperature")

################################################################################
#### Precipitation
################################################################################
# Compute average precipitation in precipitation maps
prec_summary <- global(prec, mean)

# Put data into a nice dataframe
prec_summary <- tibble(
    Date  = date(prec_dates)
  , Prec  = prec_summary$mean
  , Year  = year(Date)
  , Month = month(Date)
  , Week  = week(Date)
  , Day   = yday(Date)
  , Hour  = hour(prec_dates)
)

# Visualize temperature by day of the year
p4 <- prec_summary %>%
  group_by(Year, Day) %>%
  summarize(Prec = mean(Prec)) %>%
  ggplot(aes(x = Day, y = Prec, col = Year, group = Year)) +
  geom_line() +
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Precipitation")

# Visualize temperature by time of the day
p5 <- ggplot(prec_summary, aes(x = Hour, y = Prec, col = Month)) +
  geom_jitter(size = 0.1) +
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Precipitation")

################################################################################
#### NDVI
################################################################################
# Compute average NDVI in ndvi maps
ndvi_summary <- global(ndvi, mean)

# Put data into a nice dataframe
ndvi_summary <- tibble(
    Date  = date(ndvi_dates)
  , NDVI  = ndvi_summary$mean
  , Year  = year(Date)
  , Month = month(Date)
  , Week  = week(Date)
  , Day   = yday(Date)
)

# Visualize temperature by day of the year
p6 <- ndvi_summary %>%
  group_by(Year, Day) %>%
  summarize(NDVI = mean(NDVI)) %>%
  ggplot(aes(x = Day, y = NDVI, col = Year, group = Year)) +
  geom_line() +
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("NDVI")

################################################################################
#### Combine Plots
################################################################################
# Put plots together
ggarrange(p1, p2, p4, p6, ncol = 1, labels = "auto")
