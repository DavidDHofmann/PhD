################################################################################
#### Analysis
################################################################################
# Description: Analysis of the Activity data

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
wd <- "C:/Users/david/switchdrive/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(hms)         # To handle times
library(runner)      # To apply moving windows

# Load cleaned activity data (note that the timestamps are all in UTC)
dat <- read_csv("03_Data/02_CleanData/ActivityDataWithCovariates.csv")

# Compute some useful time metrics
dat$Hour <- as_hms(dat$Timestamp)

# Subset for now
# dat <- subset(dat, DogID == "Abel")
dat <- subset(dat, State == "Resident")

# Nest the data by dog, collar, and date
dat_nested <- dat %>% nest(Data = -c(DogID, CollarID, Date))

# Check the number of datapoints per day
dat_nested$N <- sapply(dat_nested$Data, function(x) {nrow(x)})

# Let's only keep days where we have at least 280 datapoints
dat_nested <- subset(dat_nested, N > 280)

# Can do some plots
dat_nested %>%
  slice(58) %>%
  unnest(Data) %>%
  ggplot(aes(x = Timestamp, y = ActX, col = ToD)) +
    geom_point() +
    geom_vline(aes(xintercept = max(Sunset))) +
    geom_vline(aes(xintercept = max(Sunset) - hours(2)), lty = 2) +
    geom_vline(aes(xintercept = max(Sunset) + hours(2)), lty = 2)

# Compute Average activity within 2 hours before and after sunset
dat_nested <- mutate(dat_nested, EveningActivity = map(Data, function(x) {
  tsunset <- max(x$Sunset)
  evening <- subset(x
    , Timestamp >= tsunset - hours(2) &
      Timestamp <= tsunset + hours(2)
  )
  evening <- summarize(evening
    , meanActX               = mean(ActX)
    , SDActX                 = sd(ActX)
    , Rain                   = ifelse(any(Rain), 1, 0)
    , Season                 = unique(Season)
    , minMoonlightIntensity  = min(minMoonlightIntensity)
    , meanMoonlightIntensity = mean(meanMoonlightIntensity)
    , maxMoonlightIntensity  = max(maxMoonlightIntensity)
    , minTemperature         = min(Temperature)
    , meanTemperature        = mean(Temperature)
    , maxTemperature         = max(Temperature)
    , minPrecipitation       = min(Precipitation)
    , meanPrecipitation      = mean(Precipitation)
    , maxPrecipitation       = max(Precipitation)
    , meanCloudCover         = mean(CloudCover)
  )
  return(evening)
}))

# Extract evening data
evening <- dat_nested %>%
  select(-c(N, Data)) %>%
  unnest(EveningActivity)

# Normalize values to a range between 0 and 1 (I won't standardize them because
# they are not normally distributed)
normalize <- function(x) {(x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))}
evening <- mutate(evening
  , minMoonlightIntensity  = normalize(minMoonlightIntensity)
  , meanMoonlightIntensity = normalize(meanMoonlightIntensity)
  , maxMoonlightIntensity  = normalize(maxMoonlightIntensity)
  , minTemperature         = normalize(minTemperature)
  , meanTemperature        = normalize(meanTemperature)
  , maxTemperature         = normalize(maxTemperature)
  , minPrecipitation       = normalize(minPrecipitation)
  , meanPrecipitation      = normalize(meanPrecipitation)
  , maxPrecipitation       = normalize(maxPrecipitation)
  , meanCloudCover         = normalize(meanCloudCover)
)

################################################################################
#### Exploratory
################################################################################
# Plot relationship of evening activity with every continuous covariate
evening %>%
  pivot_longer(minMoonlightIntensity:meanCloudCover, names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = Value, y = meanActX)) +
    geom_point(size = 0.4, alpha = 0.5) +
    facet_wrap(~ Variable, scale = "free") +
    stat_smooth() +
    theme_minimal()

# Also plot catgorical relationships
ggplot(evening, aes(x = as.factor(Rain), y = meanActX)) +
  geom_boxplot(size = 0.4, alpha = 0.5) +
  theme_minimal()
ggplot(evening, aes(x = as.factor(Season), y = meanActX)) +
  geom_boxplot(size = 0.4, alpha = 0.5) +
  theme_minimal()

# Nest data by individual
evening_nested <- nest(evening, Data = -c(DogID))

# How many datapoints per individual?
evening_nested$N <- sapply(evening_nested$Data, function(x) {nrow(x)})
ggplot(evening_nested, aes(x = DogID, y = N)) + geom_col()

# Only work with dogs that have at least 50+ datapoints
evening_nested <- subset(evening_nested, N >= 50)

# Decide on a model formula
formula <- quote(meanActX ~
    + maxMoonlightIntensity
    # + maxPrecipitation
    + Rain
    + maxTemperature
    + meanCloudCover
    + maxMoonlightIntensity:meanCloudCover
    # + I(maxPrecipitation ** 2)
    # + I(meanTemperature ** 2)
)

# Run a linear model for each dog separately
evening_nested$Model <- lapply(evening_nested$Data, function(x) {
  lm(formula, data = x)
})

# Extract model results
evening_nested$ModelCoefficients <- lapply(evening_nested$Model, function(x) {
  broom::tidy(x)
})

# Calculate average movement of the past 24 hours

# Use formulas to compute...

# Visualize
evening_nested %>%
  select(DogID, ModelCoefficients) %>%
  unnest(ModelCoefficients) %>%
  ggplot(aes(x = term, y = estimate, col = DogID)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(position = position_jitter(width = 0.2)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90), legend.position = "none")
