################################################################################
#### Temporal Resolution of Covariates
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(tidyverse)    # For plotting and data wrangling
library(terra)        # To handle spatial data
library(lubridate)    # To handle dates
library(hms)          # To handle times
library(ggpubr)       # To arrange multiple plots
library(sf)           # For spatial data
library(stars)        # For spatial data
library(ggnewscale)   # For multiple colorscales
library(colorspace)   # For colors
library(ggdark)       # For dark themes

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load all spatial covariates
covs1 <- "03_Data/02_CleanData/Covariates.rds" %>%
  read_rds() %>%
  subset(Type == "Dynamic") %>%
  mutate(Dates = map(Dates, function(x) {
    x$Layerdate
  })) %>%
  dplyr::select(-c(Filename, Type))

# Also load the nightly statistics
covs2 <- "03_Data/02_CleanData/Moonlight.rds" %>%
  read_rds() %>%
  pull(Timestamp) %>%
  expand_grid(Date = ., Covariate = c("MoonIllumination", "TimeOfDay")) %>%
  nest(Dates = -Covariate) %>%
  mutate(Dates = map(Dates, function(x) {
    x$Date
  }))

# Put them together
covs <- rbind(covs1, covs2)
print(covs)

################################################################################
#### Covariate Date-Ranges
################################################################################
# Function that provides a "stylized" set of timestamps given our covariate
# timestamps. This is just to get a nice figure and doesn't need to be 100%
# perfect
styledDates <- function(dates, spacing = 0.1, max_days = 365, dt_min = 1) {

  # If there is only a single date, span the entire year
  if (length(dates) == 1) {
    newdays = tibble(Start = 1, End = 366)
    return(newdays)
  }

  # Drop dates coming from leap years
  remove <- yday(dates) > 365
  dates  <- dates[!remove]

  # Convert dates to days
  dates_days <- yday(dates) + hour(dates) / 24

  # Find median time difference between dates
  dt <- lead(dates_days) - dates_days
  dt <- median(dt, na.rm = T)

  # Prepare dataframe of stylized dates
  newdays <- tibble(
      Start = seq(min(dates_days), max(dates_days), by = dt)
    , End   = Start + (1 - spacing) * dt
  )

  # Return stylized dates
  return(newdays)
}

# Define spacings
covs$Spacing = c(0.1, 0.1, 0.02, 0.02, 0.02, 0.01, 0.02, 0.005, rep(0.3, 4))
covs$Styleddates <- lapply(1:nrow(covs), function(i) {
  styledDates(covs$Dates[[i]], spacing = covs$Spacing[i], dt_min = 1)
})

# Prepare figure
p1 <- covs %>%
  dplyr::select(Covariate, Styleddates) %>%
  unnest(Styleddates) %>%
  mutate(Covariate = factor(Covariate
    , levels = c("Humans", "Forest", "Trees", "Shrubs", "Water", "DistanceToWater", "DistanceToPans", "NDVI", "Temperature", "Precipitation", "MoonIllumination", "TimeOfDay")
  )) %>%
  ggplot(aes(x = Start, xend = End, y = Covariate, ystart = Covariate, yend = Covariate)) +
    geom_segment(size = 3, col = "gray75") +
    # geom_point(aes(x = (Day + ToDay) / 2), pch = "|") +
    dark_theme_minimal() +
    xlab("Day of the Year") +
    ylab("Covariate") +
    scale_x_continuous(breaks = seq(0, 350, by = 5), limits = c(0, 366), position = "bottom") +
    scale_y_discrete(limits = rev) +
    coord_cartesian(xlim = c(10, 75)) +
    theme(
        panel.grid       = element_blank()
      , plot.background  = element_blank()
      , legend.position  = "none"
      , strip.background = element_rect(fill = "gray95", color = "white")
      , axis.line.x      = element_line()
      , axis.ticks.x     = element_line()
      , axis.text.y      = element_text(size = 10)
    )

# Store the plot
ggsave("05_Presentation/TemporalResolution.png"
  , plot   = p1
  , width  = 10
  , height = 2.5
  , scale  = 1
  , bg     = "transparent"
  , device = png
)
