################################################################################
#### Exploratory Analysis
################################################################################
# Description: Exploration of the data, prior to fitting the movement model

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Load required packages
library(tidyverse)    # For data wrangling

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Loading Data
################################################################################
# Load data
dat <- read_rds("03_Data/02_CleanData/SSFExtracted.rds")
print(names(dat))

# Keep only desired columns
dat <- rename(dat, id = ID, ta = relta)

# Add log of the step length, cos of the turning angle, and sqrt of the
# distance to water, and a binary indicator of night
dat <- dat %>% mutate(
    log_sl              = log(sl)
  , cos_ta              = cos(ta)
)

# Let's also move all movement metrics to the front
dat <- dat %>% select(c(
  id, step_id, step_id_within, case, sl, log_sl, ta, cos_ta, inactive, everything()
))

dat <- mutate(dat, Hour = as.factor(hour(TimestampRounded)))

################################################################################
#### Exploration
################################################################################
# Prepare a dataframe of observed steps
dat_obs <- subset(dat, case == 1)

# Are there differences between the sexes?
dat_obs %>%
  group_by(Sex) %>%
  summarize(sl = mean(sl), cos_ta = mean(cos_ta))

# Show step length in relation to time of the day
ggplot(dat_obs, aes(x = Hour, y = sl, col = Hour, fill = Hour)) +
  geom_violin(alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_awesome()

# Do they move particularly directional at a certain time?
ggplot(dat_obs, aes(x = Hour, y = cos_ta, col = Hour, fill = Hour)) +
  geom_violin(alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_awesome()

# Quantify step lengths and turning angles by phases
dat_obs %>%
  group_by(Hour) %>%
  summarize(sl = mean(sl), cos_ta = mean(cos_ta))

# Compare overall activity during day and night
ggplot(dat_obs, aes(x = Night, y = sl, col = Night, fill = Night)) +
  geom_violin(alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_awesome()

# Compare overall activity during day and night
ggplot(dat_obs, aes(x = Night, y = cos_ta, col = Night, fill = Night)) +
  geom_violin(alpha = 0.5) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_awesome()

# Let's get some numbers for this
dat_obs %>%
  group_by(Night) %>%
  summarize(sl = mean(sl), cos_ta = mean(cos_ta))

# Let's try to distinguish moonlit, from pitch-dark nights (one can pick a
# different threshold)
dat_obs %>%
  subset(Night) %>%
  mutate(Illuminated = meanMoonlight > 0.2) %>%
  group_by(Illuminated) %>%
  summarize(sl = mean(sl), cos_ta = mean(cos_ta))

# Let's try to distinguish moonlit, from pitch-dark nights, but also consider if
# there are clouds or not
dat_obs %>%
  mutate(Illuminated = meanMoonlight > 0.2) %>%
  mutate(Rainy = PrecipitationDynamic > 0) %>%
  group_by(Night, Illuminated, Rainy) %>%
  summarize(sl = mean(sl), cos_ta = mean(cos_ta), .groups = "drop") %>%
  mutate(LightType = case_when(
      !Night ~ "Day"
    , Night & Illuminated ~ "BrightNight"
    , .default = "DarkNight"
  )) %>%
  ggplot(aes(x = Rainy, y = sl)) +
    geom_col() +
    facet_wrap(~ LightType, nrow = 1) +
    theme_awesome()

# Is directionality correlated with step-length? Let's do some non-parametric
# tests.
cor.test(dat_obs$sl, dat_obs$cos_ta, method = "spearman")
cor.test(dat_obs$sl, dat_obs$cos_ta, method = "kendall")

# What if directionality is binary (directional / not directional)?
dat_obs %>%
  mutate(Directional = cos_ta > 0) %>%
  group_by(Directional) %>%
  summarize(sl = mean(sl))

# Are they more willing to venture into human-dominated areas at night?
dat_obs %>%
  group_by(Night) %>%
  summarize(Humans = mean(HumansStatic))

# Are they more willing to venture into human-dominated areas at night when
# there is moonlight?
dat_obs %>%
  mutate(Illuminated = meanMoonlight > 0.2) %>%
  group_by(Night, Illuminated) %>%
  summarize(Humans = mean(HumansStatic), .groups = "drop")

# Check how step-length and turning-angle are affected by water
ggplot(dat_obs, aes(x = WaterDynamic, y = sl)) +
  geom_point() +
  geom_smooth() +
  theme_awesome()

ggplot(dat_obs, aes(x = WaterDynamic, y = cos_ta)) +
  geom_point() +
  geom_smooth() +
  theme_awesome()

# Check how step-length and turning-angle are affected by woodland
ggplot(dat_obs, aes(x = ForestStatic, y = sl)) +
# ggplot(dat_obs, aes(x = TreesDynamic, y = sl)) +
  geom_point() +
  geom_smooth() +
  theme_awesome()

ggplot(dat_obs, aes(x = ForestStatic, y = cos_ta)) +
# ggplot(dat_obs, aes(x = TreesDynamic, y = cos_ta)) +
  geom_point() +
  geom_smooth() +
  theme_awesome()

# Check how step-length and turning-angle are affected by "DistanceToPans"
ggplot(dat_obs, aes(x = DistanceToPansDynamic, y = sl)) +
  geom_point() +
  geom_smooth() +
  theme_awesome()

ggplot(dat_obs, aes(x = DistanceToPansDynamic, y = cos_ta)) +
  geom_point() +
  geom_smooth() +
  theme_awesome()

# Check how step-length and turning-angle are affected by "DistanceToWater"
ggplot(dat_obs, aes(x = DistanceToWaterDynamic, y = sl)) +
  geom_point() +
  geom_smooth() +
  theme_awesome()

ggplot(dat_obs, aes(x = DistanceToWaterDynamic, y = cos_ta)) +
  geom_point() +
  geom_smooth() +
  theme_awesome()

# What about the NDVI?
ggplot(dat_obs, aes(x = NDVIDynamic, y = sl)) +
  geom_point() +
  geom_smooth() +
  theme_awesome()

ggplot(dat_obs, aes(x = NDVIDynamic, y = cos_ta)) +
  geom_point() +
  geom_smooth() +
  theme_awesome()

# Check how step-length and turning-angle are affected by precipitation
ggplot(dat_obs, aes(x = PrecipitationDynamic, y = sl)) +
  geom_point() +
  geom_smooth() +
  theme_awesome() +
  facet_wrap(~ Night)

ggplot(dat_obs, aes(x = PrecipitationDynamic, y = cos_ta)) +
  geom_point() +
  geom_smooth() +
  theme_awesome() +
  facet_wrap(~ Night)

# Check how step-length and turning-angle are affected by temperature
ggplot(dat_obs, aes(x = TemperatureDynamic, y = sl)) +
  geom_point() +
  geom_smooth() +
  theme_awesome()

ggplot(dat_obs, aes(x = TemperatureDynamic, y = cos_ta)) +
  geom_point() +
  geom_smooth() +
  theme_awesome() +
  facet_wrap(~ Night)

# I imagined that forest / woodland should be preferred when it is hot
ggplot(dat_obs, aes(x = TemperatureDynamic, y = ForestStatic)) +
  geom_point() +
  geom_smooth() +
  theme_awesome() +
  facet_wrap(~ Night)

ggplot(dat_obs, aes(x = TemperatureDynamic, y = TreesDynamic)) +
  geom_point() +
  geom_smooth() +
  theme_awesome() +
  facet_wrap(~ Night)

# Similarly, I'd expect dispersers to spend more time near water when it is hot
ggplot(dat_obs, aes(x = TemperatureDynamic, y = DistanceToWaterDynamic)) +
  geom_point() +
  geom_smooth() +
  theme_awesome() +
  facet_wrap(Night ~ SeasonClimate)

ggplot(dat_obs, aes(x = TemperatureDynamic, y = DistanceToPansDynamic)) +
  geom_point() +
  geom_smooth() +
  theme_awesome() +
  facet_wrap(Night ~ SeasonClimate)

# Are there any obvious differences between the seasons?
dat_obs %>%
  group_by(SeasonClimate) %>%
  summarize(sl = mean(sl))

dat_obs %>%
  group_by(SeasonClimate) %>%
  summarize(cos_ta = mean(cos_ta))

# Is there a tendency for dispersers to reduce their step-lengths over time?
dat_obs %>%
  group_by(id) %>%
  mutate(step_number = 1:n()) %>%
  ggplot(aes(x = step_number, y = sl)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~ id, scales = "free") +
    theme_awesome()

# Maybe there is a tendency to move less directed?
dat_obs %>%
  group_by(id) %>%
  mutate(step_number = 1:n()) %>%
  ggplot(aes(x = step_number, y = cos_ta)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~ id, scales = "free") +
    theme_awesome()
