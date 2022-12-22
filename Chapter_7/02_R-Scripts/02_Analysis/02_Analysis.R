################################################################################
#### Analysis of Average Activity and Time of Starting Moving
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
dat         <- read_csv("03_Data/02_CleanData/ActivityDataCovariatesAggregated.csv")
dat$NFixes  <- NULL
print(names(dat))

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

# Nest data by individual and tod
dat_nested <- dat %>% nest(Data = -c(DogID, ToD))

# Only work with rows for which there are at least 50 datapoints
dat_nested$N <- sapply(dat_nested$Data, function(x) {nrow(x)})
dat_nested <- subset(dat_nested, N >= 50)
print(dat_nested)

################################################################################
#### Exploratory
################################################################################
# How does activity depend on moonlight intensity?
ggplot(dat, aes(x = maxMoonlightIntensity, y = meanActX)) +
  geom_point() +
  facet_wrap(~ ToD) +
  theme_minimal() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2))

# How does activity depend on precipitation?
ggplot(dat, aes(x = maxPrecipitation, y = meanActX)) +
  geom_point() +
  facet_wrap(~ ToD) +
  theme_minimal() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2))

# How does activity depend on rain (binary)?
ggplot(dat, aes(x = as.factor(Rain), y = meanActX)) +
  geom_boxplot() +
  facet_wrap(~ ToD) +
  theme_minimal()

# Is there an interaction between the two?
ggplot(dat, aes(x = maxPrecipitation, y = meanActX, col = as.factor(Rain))) +
  geom_point() +
  facet_wrap(~ ToD) +
  theme_minimal() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2))

# How does activity depend on temperature?
ggplot(dat, aes(x = maxTemperature, y = meanActX)) +
  geom_point() +
  facet_wrap(~ ToD) +
  theme_minimal() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2))

# How does activity depend on cloud cover?
ggplot(dat, aes(x = meanCloudCover, y = meanActX)) +
  geom_point() +
  facet_wrap(~ ToD) +
  theme_minimal() +
  geom_smooth()

# How does activity depend on cloud cover in the following night?
ggplot(dat, aes(x = meanCloudCoverNight, y = meanActX)) +
  geom_point() +
  facet_wrap(~ ToD) +
  theme_minimal() +
  geom_smooth()

# Is there an interaction between cloud cover and moonlight intensity?
ggplot(dat, aes(x = maxMoonlightIntensity, y = meanActX, col = Cloudy)) +
  # geom_point() +
  facet_wrap(~ ToD) +
  theme_minimal() +
  geom_smooth()

# How does activity depend on the level of activity in the previous xx hours?
dat %>%
  select(meanActX, ToD, meanActX6:meanActX24) %>%
  pivot_longer(meanActX6:meanActX24) %>%
  ggplot(aes(x = value, y = meanActX, col = name)) +
    # geom_point() +
    facet_wrap(~ ToD) +
    theme_minimal() +
    geom_smooth()

################################################################################
#### Analysis I - Average Activity by Time of the Day
################################################################################
# Decide on a model formula
formula <- quote(meanActX ~
    + maxMoonlightIntensity
    + I(maxMoonlightIntensity ** 2)
    # + Rain
    # + maxPrecipitation
    # + Rain:maxPrecipitation
    + maxTemperature
    # + I(maxTemperature ** 2)
    # + meanCloudCover
    + meanCloudCoverNight
    + maxMoonlightIntensity:meanCloudCoverNight
    # + maxMoonlightIntensity:maxMoonDelay
    + meanActX6
    # + I(maxPrecipitation ** 2)
    # + I(meanTemperature ** 2)
)

# Run a linear model for each dog separately for the different times of the day
dat_nested$Model <- lapply(dat_nested$Data, function(x) {
  lm(formula, x)
  # Could use a tobit regression as well
  # coef(tobit(formula, data = x, left = 0, right = 255)
})

# Extract model results
dat_nested$ModelCoefficients <- lapply(dat_nested$Model, function(x) {
  broom::tidy(x)
  # If tobit
  # tidy(summary(x)$coefficients)
})

# Unnest
coefs <- dat_nested %>%
  select(DogID, ToD, ModelCoefficients) %>%
  unnest(ModelCoefficients)

# Write this to file
write_rds(coefs, "03_Data/03_Results/99_MeanActivityModelCoefficients.rds")

################################################################################
#### Analysis II - Time when Becoming Moving
################################################################################
# For this, we'll only work with "evening" data
evening <- subset(dat, ToD == "Evening")

# Compute time of the day as regular seconds since midnight
evening$StartMovingNumeric <- as.numeric(as_hms(evening$StartMoving))

# Nest data by individual and tod
evening_nested <- evening %>% nest(Data = -c(DogID, ToD))

# Only work with rows for which there are at least 50 datapoints
evening_nested$N <- sapply(evening_nested$Data, function(x) {nrow(x)})
evening_nested <- subset(evening_nested, N >= 50)
print(evening_nested)

# Some exploratory plots. How does the time of activity depend on the moonlight
# intensity of the upcoming night and the delay of the moon?
ggplot(evening, aes(x = maxMoonlightIntensity, y = as.numeric(as_hms(StartMoving)), col = maxMoonDelay)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  geom_smooth(method = "lm") +
  scale_color_viridis_c(option = "magma")

# Let's separate between cloudy and clear nights
evening %>%
  mutate(Cloudy = meanCloudCoverNight > 0.5) %>%
  ggplot(aes(x = maxMoonlightIntensity, y = as.numeric(as_hms(StartMoving)), col = Cloudy)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    geom_smooth(method = "lm") +
    scale_color_manual(values = c("orange", "cornflowerblue"))

# How does the time of becoming active change throghout the year?
ggplot(evening, aes(x = yday(Date), y = as_hms(StartMoving))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  theme_minimal()

# In some cases the dogs were already moving at 12:00. What if we ignore those?
ggplot(subset(evening, as_hms(StartMoving) > as_hms("14:00:00")), aes(x = yday(Date), y = as_hms(StartMoving))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
  theme_minimal()

# Decide on a model formula
formula <- quote(StartMovingNumeric ~
    + maxMoonlightIntensity
    + I(maxMoonlightIntensity ** 2)
    # + Rain
    # + maxPrecipitation
    # + Rain:maxPrecipitation
    + maxTemperature
    + maxTemperature:maxMoonlightIntensity
    # + I(maxTemperature ** 2)
    # + meanCloudCover
    + meanCloudCoverNight
    + maxMoonlightIntensity:meanCloudCoverNight
    + maxMoonlightIntensity:maxMoonDelay
    + meanActX12
    # + I(maxPrecipitation ** 2)
    # + I(meanTemperature ** 2)
)

# Run a linear model for each dog separately for the different times of the day
evening_nested$Model <- lapply(evening_nested$Data, function(x) {
  lm(formula, x)
})

# Extract model results
evening_nested$ModelCoefficients <- lapply(evening_nested$Model, function(x) {
  broom::tidy(x)
})

# Unnest
coefs <- evening_nested %>%
  select(DogID, ToD, ModelCoefficients) %>%
  unnest(ModelCoefficients)

# Compute weight for each coefficient
coefs_means_weighted <- coefs %>%
  nest(Coefs = -c(ToD, term)) %>%
  mutate(Coefs = map(Coefs, function(x) {
    be <- x$estimate
    se <- x$std.error
    we <- (1 / (se ** 2)) / sum(1 / (se ** 2))
    be_mean <- sum(we * be)
    se_mean <- sqrt(sum(we * (be - be_mean) ** 2) / (length(be) - 1))
    result <- tibble(
        estimate  = be_mean
      , std.error = se_mean
      , LWR       = estimate - 1.96 * std.error
      , UPR       = estimate + 1.96 * std.error
    )
    return(result)
  })) %>% unnest(Coefs)

# Compute mean and se of the model estimates (unweighted)
coefs_means_unweighted <- coefs %>%
  group_by(ToD, term) %>%
  summarize(
      estimate_mean = mean(estimate)
    , estimate_sd   = sd(estimate) / sqrt(n())
    , LWR           = estimate_mean - qt(1 - 0.05 / 2, n() - 1) * estimate_sd  # UNWEIGHTED: TAKE A LOOK AT MURTAUGH AGAIN
    , UPR           = estimate_mean + qt(1 - 0.05 / 2, n() - 1) * estimate_sd  # UNWEIGHTED: TAKE A LOOK AT MURTAUGH AGAIN
    , .groups       = "drop"
  )

# Visualize
coefs %>%
  ggplot(aes(x = term, y = estimate, col = DogID)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_jitter(width = 0.1, alpha = 0.8) +
    geom_point(
        data        = coefs_means_unweighted
      , mapping     = aes(x = term, y = estimate_mean)
      , inherit.aes = F
      , size        = 2
    ) +
    geom_errorbar(
        data        = coefs_means_unweighted
      , mapping     = aes(x = term, ymin = LWR, ymax = UPR)
      , inherit.aes = F
      , width       = 0.2
      , lwd         = 1
    ) +
    theme_minimal() +
    scale_color_viridis_d() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("04_Manuscript/99_TimeBecomingActive.png", plot = last_plot(), bg = "white")
