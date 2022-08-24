################################################################################
#### Analysis of Average Activity
################################################################################
# Description: Analysis of average activity at different times of the day

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(hms)         # To handle times
library(broom)       # To clean model outputs

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load cleaned activity data (note that the timestamps are all in UTC)
dat <- read_csv("03_Data/02_CleanData/ActivityDataWithCovariatesAggregated.csv")
dat$NFixes <- NULL

# For now, ignore dispersers
dat <- subset(dat, State == "Resident")

# Normalize covariate values to a range between 0 and 1 (I won't standardize
# them because they are not normally distributed)
dat <- mutate(dat
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
  , meanActX6              = normalize(meanActX6)
  , meanActX12             = normalize(meanActX12)
  , meanActX18             = normalize(meanActX18)
  , meanActX24             = normalize(meanActX24)
  , meanActX24             = normalize(meanCloudCoverNight)
  , maxMoonDelay           = normalize(maxMoonDelay)
)

# Nest data by individual and state
dat_nested <- dat %>% nest(Data = -c(DogID, ToD))

# Only work with rows for which there are at least 50 datapoints
dat_nested$N <- sapply(dat_nested$Data, function(x) {nrow(x)})
dat_nested <- subset(dat_nested, N >= 50)

# Decide on a model formula
formula <- quote(meanActX ~
    + maxMoonlightIntensity
    + Rain
    # + maxPrecipitation
    # + Rain:maxPrecipitation
    + maxTemperature
    # + meanCloudCover
    # + meanCloudCoverNight
    # + maxMoonlightIntensity:meanCloudCoverNight
    + maxMoonlightIntensity:maxMoonDelay
    + meanActX24
    # + I(maxPrecipitation ** 2)
    # + I(meanTemperature ** 2)
)

# Run a linear model for each dog separately for the different times of the day
dat_nested$Model <- lapply(dat_nested$Data, function(x) {
  lm(formula, x)
})

# Extract model results
dat_nested$ModelCoefficients <- lapply(dat_nested$Model, function(x) {
  broom::tidy(x)
})

# Visualize
dat_nested %>%
  select(DogID, ToD, ModelCoefficients) %>%
  unnest(ModelCoefficients) %>%
  ggplot(aes(x = ToD, y = estimate, col = DogID)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_point(position = position_jitter(width = 0.2)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
    facet_wrap(~ term, scale = "free") +
    scale_color_viridis_d()

















# Visualize how the mean activity depends on moonlight intensity
ggplot(dat, aes(x = maxMoonlightIntensity, y = meanActX)) +
  geom_point(alpha = 0.2) +
  theme_minimal() +
  facet_wrap(~ State + ToD) +
  stat_smooth()

# Let's be a bit more elaborate and run the regression for each individual
# separately
mods <- dat %>%
  nest(Data = -c(ToD, DogID, State)) %>%
  mutate(Model = map(Data, function(x) {
    lm(meanActX ~ maxMoonlightIntensity, data = x)
  })) %>%
  mutate(Coefs = map(Model, function(x) {
    tidy(x)
  }))

# Compute mean and se of the model estimates
coefs_means <- mods %>%
  select(-c(Data, Model)) %>%
  unnest(Coefs) %>%
  group_by(ToD, State, term) %>%
  summarize(
      MeanEstimate = mean(estimate)
    , SE           = sd(estimate) / sqrt(n())
    , LWR          = MeanEstimate - qt(1 - 0.05 / 2, n() - 1) * SE  # TAKE A LOOK AT MURTAUGH AGAIN
    , UPR          = MeanEstimate + qt(1 - 0.05 / 2, n() - 1) * SE  # TAKE A LOOK AT MURTAUGH AGAIN
    , .groups      = "drop"
  )

# Unnest the coefficients and plot them
mods %>%
  select(-c(Data, Model)) %>%
  unnest(Coefs) %>%
  ggplot(aes(x = ToD, y = estimate, col = DogID)) +
    geom_point(
        data        = coefs_means
      , mapping     = aes(x = ToD, y = MeanEstimate)
      , inherit.aes = F
      , size        = 4
    ) +
    geom_errorbar(
        data        = coefs_means
      , mapping     = aes(x = ToD, ymin = LWR, ymax = UPR)
      , inherit.aes = F
      , width       = 0.3
      , lwd         = 1.5
    ) +
    # geom_violin(col = "black") + # Could also use a violin plot
    geom_jitter(width = 0.1) +
    facet_wrap(~ State + term, scales = "free") +
    theme_minimal() +
    scale_color_viridis_d() +
    theme(legend.position = "none")

# Is there an issue of correlation of the residuals?
mods <- mods %>%
  mutate(Residuals = map(Model, function(x) {
    resids <- residuals(x)
    return(resids)
  })) %>%
  mutate(ACF = map(Residuals, function(x) {
    acfvalues <- acf(x, plot = F)
    return(acfvalues)
  }))

# Let's put all residuals together and plot one single acf
mods %>%
  select(Residuals) %>%
  unnest(Residuals) %>%
  pull(Residuals) %>%
  acf()
