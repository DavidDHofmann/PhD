################################################################################
#### Daytime vs. Nighttime Activity
################################################################################
# Description: Under what lighting conditions is nighttime activity comparable
# to daytime activity?

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(hms)         # To handle times

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load data
dat <- read_csv("03_Data/02_CleanData/ActivityDataCovariates.csv")

# There are two columns we can work with. We'll take ToD2 as it's binary
table(dat$ToD)
table(dat$ToD2)

# Nest data by dog
dat <- nest(dat, Data = -DogID)

# Create tuples of night and day
dat <- mutate(dat, Data = map(Data, function(x) {
  x <- dat$Data[[1]]
  x$Switch <- x$ToD2 != lag(x$ToD2)
  x$Switch[1] <- F
  x$Switch <- x$Switch & x$ToD2 == "Day"
  x$Cycle  <- cumsum(x$Switch) + 1
  x$Switch <- NULL
  return(x)
}))

# Nest by ID and cycle
dat <- unnest(dat, Data)
print(names(dat))

# Compute statistics by ID and cycle
dat %>%
  group_by(DogID, Cycle, ToD2) %>%
  summarize(Act = mean(Act), maxMoonlight = max(maxMoonlightIntensity), .groups = "drop") %>%
  mutate(maxMoonlight = cut(maxMoonlight, breaks = 5)) %>%
  ggplot(aes(x = maxMoonlight, y = Act)) +
    geom_boxplot() +
    facet_wrap(~ToD2)
