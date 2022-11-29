################################################################################
#### Exploratory Analysis
################################################################################
# Exploratory Analysis of the Activity Data

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(hms)         # To handle times
library(pracma)      # To identify peaks
library(broom)       # To clean model summary
library(ggpubr)      # To arrange multiple plots

# Reload cleaned activity data (note that the timestamps are all in UTC)
dat <- read_csv("03_Data/02_CleanData/ActivityDataCovariates.csv")
print(names(dat))

# How many datapoints for residents and dispersers are there?
table(dat$State)
prop.table(table(dat$State)) * 100

# For now, only look at residents
dat <- subset(dat, State != "Disperser")
gc()

# Get a feeling of the data
dim(dat)
nrow(dat) / 1e6
str(dat)
glimpse(dat)

# How many individuals?
length(unique(dat$DogID))

# How many datapoints per individual?
dat %>%
  count(DogID) %>%
  ggplot(aes(x = DogID, y = n)) +
  geom_col(fill = NA, col = "cornflowerblue", width = 0.1) +
  geom_point(col = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

# Compute summary
dat %>%
  count(DogID) %>%
  summarize(
      MeanNumberFixes   = mean(`n`)
    , SDNumberFixes     = sd(`n`)
    , MedianNumberFixes = median(`n`)
    , IQRNumberFixes    = IQR(`n`)
    , MinNumberFixes    = min(`n`)
    , MaxNumberFixes    = max(`n`)
  )

# How active where the dogs on average?
dat %>%
  group_by(DogID) %>%
  summarize(
      MeanActX   = mean(ActX)
    , SDActX     = sd(ActX)
    , MedianActX = median(ActX)
    , IQRActX    = IQR(ActX)
  )

# Aggregate activity by day and see if there are differences between the dogs
dat %>%
  group_by(DogID, Date) %>%
  summarize(AverageDailyActivity = mean(ActX), n = n(), .groups = "drop") %>%
  subset(n > 250) %>%
  ggplot(aes(x = DogID, y = AverageDailyActivity, fill = DogID, col = DogID)) +
    geom_boxplot(alpha = 0.3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
    scale_fill_viridis_d() +
    scale_color_viridis_d()

# Aggregate activity by month
dat %>%
  subset(DogID %in% c("Abel", "Amacuro", "Aramis", "Aspen", "Atwood")) %>%
  mutate(Month = month(Timestamp), Hour = hour(Timestamp)) %>%
  group_by(Month, Hour, DogID) %>%
  summarize(MeanHourlyAct = mean(ActX), .groups = "drop") %>%
  ggplot(aes(x = Hour, y = MeanHourlyAct, col = as.factor(Month))) +
    geom_line(alpha = 0.75) +
    facet_wrap(~DogID, ncol = 1) +
    theme_minimal() +
    scale_color_viridis_d() +
    theme(legend.position = "bottom")

# Remove undesired columns
dat <- select(dat, -c(DOP, minMoonPhase, meanMoonPhase, maxMoonPhase))

# Look at the distributions of the different numerical variables (I'll not use
# ggplot as it takes way too long)
par(mfrow = c(3, 3), mar = c(5, 5, 1, 1))
with(dat, expr = {
  hist(ActX, col = "cornflowerblue", border = "white", main = "")
  hist(minMoonlightIntensity, col = "cornflowerblue", border = "white", main = "")
  hist(meanMoonlightIntensity, col = "cornflowerblue", border = "white", main = "")
  hist(maxMoonlightIntensity, col = "cornflowerblue", border = "white", main = "")
  hist(Precipitation, col = "cornflowerblue", border = "white", main = "")
  hist(Temperature, col = "cornflowerblue", border = "white", main = "")
  hist(CloudCover, col = "cornflowerblue", border = "white", main = "")
})

# We will clearly have to aggregate or transform data at some point! Otherwise
# the distributions will be very nasty to work with.

################################################################################
#### Active vs Inactive
################################################################################
# We deem an animal active whenever its activity is above 20
dat$Active <- dat$ActX > 20

# How many points during activity and inactivity?
table(dat$Active)
prop.table(table(dat$Active)) * 100

# Make sure the data is arranged properly
dat <- arrange(dat, DogID, CollarID, Timestamp)

# Nest the data by dog, collar, and date
dat <- dat %>% nest(Data = -c(DogID, CollarID, Date))

# Check the number of datapoints per day
dat$N <- sapply(dat$Data, function(x) {nrow(x)})

# Let's only keep days where we have at least 280 datapoints
dat <- subset(dat, N > 280)

# Unnest the data for further calculations
dat <- dat %>% unnest(Data)

# Compute some useful time metrics
dat$Month <- month(dat$Timestamp)
dat$Year  <- year(dat$Timestamp)
dat$Hour  <- hour(dat$TimestampRounded)
dat$Time  <- as_hms(dat$Timestamp)

# Let's see how well the data is distributed in time
dat %>%
  count(Month, Year) %>%
  ggplot(aes(x = Month, y = Year, fill = n)) +
    geom_raster() +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal() +
    scale_y_continuous(breaks = c(2016:2022))

# We can also count by month
dat %>%
  count(Month, Year) %>%
  mutate(Date = ym(paste(Year, Month))) %>%
  ggplot(aes(x = Date, y = `n`)) +
    geom_col(fill = "cornflowerblue") +
    theme_minimal()

# How does activity change during different times of the day?
dat %>%
  group_by(Hour) %>%
  summarize(MeanActX = mean(ActX), MedianActX = median(ActX)) %>%
  pivot_longer(MeanActX:MedianActX, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Hour, y = Value, fill = as.factor(Metric))) +
    geom_col(position = position_dodge()) +
    theme_minimal() +
    scale_fill_manual(values = c("cornflowerblue", "orange"), name = "Metric")

# Plot activity against moonlight intensity
ggplot(dat, aes(x = Time, y = maxMoonlightIntensity, z = ActX)) +
  stat_summary_2d(fun = "mean", bins = 100) +
  scale_fill_viridis_c(option = "magma", name = "Activity") +
  theme_minimal()

# Plot activity against date
ggplot(dat, aes(x = Time, y = Date, z = ActX)) +
  stat_summary_2d(fun = "mean", bins = 100) +
  scale_fill_viridis_c(option = "magma", name = "Activity") +
  theme_minimal() +
  scale_y_date(date_breaks = "3 months")

# Plot activity against month, regardless of the year
ggplot(dat, aes(x = Time, y = yday(Timestamp), z = ActX)) +
  stat_summary_2d(fun = "mean", bins = 100) +
  scale_fill_viridis_c(option = "magma", name = "Activity") +
  theme_minimal()

# Let's get an idea of the usual activity patterns by computing average activity
# for the different times.
means_hour <- dat %>%
  group_by(Time) %>%
  summarize(
      Mean    = mean(ActX)
    , SD      = sd(ActX)
    # , LWR     = mean(ActX) - qt(1 - 0.05 / 2, (n() - 1)) * sd(ActX) / sqrt(n())
    # , UPR     = mean(ActX) + qt(1 - 0.05 / 2, (n() - 1)) * sd(ActX) / sqrt(n())
    # , LWR     = quantile(ActX, 0.025)
    # , UPR     = quantile(ActX, 0.975)
    , .groups = "drop"
  ) %>%
  ungroup()

# Visualize it
ggplot(means_hour, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Activity")

# We can try to identify the peaks and valleys in this plot (ignoring that these
# will change as the year progresses)
peaks <- findpeaks(means_hour$Mean, minpeakdistance = 50, npeaks = 2)
peaks <- subset(means_hour, Mean %in% peaks[, 1])
valls <- findpeaks(-means_hour$Mean, minpeakdistance = 50, npeaks = 4)
valls <- subset(means_hour, -Mean %in% valls[, 1])
valls <- valls[-nrow(valls), ]

# Let's use them to distinguish morning and evening activity (i'll make the
# endtime for the morning burst a bit earlier to have a gap between morning and
# evening burst)
phase <- tibble(
    Burst     = factor(c("Morning", "Evening", "Night", "Night")
      , levels = c("Morning", "Evening", "Night"))
  , Starttime = c(valls$Time[1], valls$Time[2], valls$Time[3], as_hms("00:00:00"))
  , Endtime   = c(as_hms(valls$Time[2] - duration(2, units = "hours")), valls$Time[3], as_hms("24:00:00"), valls$Time[1])
)

# Visualize again
ggplot(means_hour, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_rect(data = phase, inherit.aes = F, aes(xmin = Starttime, xmax = Endtime, ymin = -Inf, ymax = Inf, fill = Burst), alpha = 0.15) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_vline(data = peaks, aes(xintercept = Time), col = "gray20", lty = 2) +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = "bottom") +
  scale_fill_manual(values = c("cornflowerblue", "orange", "purple"))

# How does the activity pattern change depending on the month?
means_month <- dat %>%
  group_by(Month, Time) %>%
  summarize(
      Mean    = mean(ActX)
    , SD      = sd(ActX)
    , .groups = "drop"
  ) %>%
  ungroup()

# Find peaks each month
peaks <- means_month %>%
  nest(Data = -c(Month)) %>%
  mutate(Data = map(Data, function(x) {
    peaks <- findpeaks(x$Mean, minpeakdistance = 50, npeaks = 2)
    peaks <- subset(x, Mean %in% peaks[, 1])
    return(peaks)
  })) %>%
  unnest(Data)

# Visualize
ggplot(means_month, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_rect(data = phase, inherit.aes = F, aes(xmin = Starttime, xmax = Endtime, ymin = -Inf, ymax = Inf, fill = Burst), alpha = 0.15) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_vline(data = peaks, aes(xintercept = Time), col = "gray20", lty = 2) +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  facet_wrap(~ factor(month.abb[Month], levels = month.abb), ncol = 2, dir = "v") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = "bottom") +
  scale_fill_manual(values = c("cornflowerblue", "orange", "purple"))

# How does the activity pattern change depending on the year?
means_year <- dat %>%
  group_by(Year, Time) %>%
  summarize(
      Mean    = mean(ActX)
    , SD      = sd(ActX)
    , .groups = "drop"
  ) %>%
  ungroup()

# Find peaks each month
peaks <- means_year %>%
  nest(Data = -c(Year)) %>%
  mutate(Data = map(Data, function(x) {
    peaks <- findpeaks(x$Mean, minpeakdistance = 50, npeaks = 2)
    peaks <- subset(x, Mean %in% peaks[, 1])
    return(peaks)
  })) %>%
  unnest(Data)

# Visualize
ggplot(means_year, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_rect(data = phase, inherit.aes = F, aes(xmin = Starttime, xmax = Endtime, ymin = -Inf, ymax = Inf, fill = Burst), alpha = 0.15) +
  geom_ribbon(alpha = 0.2) +
  geom_line() +
  geom_vline(data = peaks, aes(xintercept = Time), col = "gray20", lty = 2) +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  facet_wrap(~ Year, ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = "bottom") +
  scale_fill_manual(values = c("cornflowerblue", "orange", "purple"))

# Use the cutoff times to assign to each activity fix "Morning", "Evening",
# "Night", or NA
dat <- dat %>% mutate(ActivityPhase = case_when(
    Time >= phase$Starttime[1] & Time <= phase$Endtime[1] ~ "Morning"
  , Time >= phase$Starttime[2] & Time <= phase$Endtime[2] ~ "Evening"
  , Time >= phase$Starttime[3] & Time <= phase$Endtime[3] ~ "Night"
  , Time >= phase$Starttime[4] & Time <= phase$Endtime[4] ~ "Night"
  , TRUE ~ "Ignore"
))

# Make sure it worked
sort(unique(dat[dat$ActivityPhase == "Morning", ]$Hour))
sort(unique(dat[dat$ActivityPhase == "Evening", ]$Hour))
sort(unique(dat[dat$ActivityPhase == "Night", ]$Hour))
sort(unique(dat[dat$ActivityPhase == "Ignore", ]$Hour))

# Let's compute summary statistics during those phases by dog, collar and day
dat_byphase <- subset(dat, ActivityPhase != "Ignore") %>%
  group_by(DogID, CollarID, Date, ActivityPhase) %>%
  summarize(
      meanActivity           = mean(ActX)
    , medianActivity         = median(ActX)
    , totalActivity          = sum(ActX)
    , maxActivity            = max(ActX)
    , minActivity            = min(ActX)
    , minMoonlightIntensity  = min(minMoonlightIntensity)
    , meanMoonlightIntensity = mean(meanMoonlightIntensity)
    , maxMoonlightIntensity  = max(maxMoonlightIntensity)
    , Temperature            = mean(Temperature)
    , Precipitation          = mean(Precipitation)
    , CloudCover             = mean(CloudCover)
    , .groups                = "drop"
  ) %>%
  mutate(ActivityPhase = factor(ActivityPhase, levels = c("Morning", "Evening", "Night")))

# Plot the different summary statistics by daytime phase
dat_byphase %>%
  pivot_longer(meanActivity:minActivity, names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = ActivityPhase, y = Value, fill = ActivityPhase)) +
    geom_violin(alpha = 0.5, col = NA) +
    facet_wrap(~ Variable, scales = "free") +
    theme_minimal() +
    scale_fill_manual(values = c("cornflowerblue", "orange", "gray"))

# We will either work with the mean or the total activity, let's remove the rest
dat_byphase <- select(dat_byphase, -c(maxActivity, medianActivity, minActivity))

# Plot the different summary statistics by daytime phase
dat_byphase %>%
  count(ActivityPhase)
dat_byphase %>%
  ggplot(aes(x = ActivityPhase, y = meanActivity, fill = ActivityPhase)) +
    geom_violin(alpha = 0.5, col = NA) +
    theme_minimal() +
    scale_fill_manual(values = c("cornflowerblue", "orange", "gray"))

# How well is the correlation between mean and total activity?
cor(dat_byphase$meanActivity, dat_byphase$totalActivity)

# Let's bin the nights by moon illumination
dat_byphase <- mutate(dat_byphase, MoonBin = cut(maxMoonlightIntensity, breaks = 5))
dat_byphase$MoonBinNumeric <- as.numeric(dat_byphase$MoonBin)

# Visualize how the mean activity depends on moonlight intensity
p1 <- dat_byphase %>%
  group_by(MoonBinNumeric, ActivityPhase) %>%
  summarize(
      n       = n()
    , mean    = round(mean(meanActivity))
    , .groups = "drop"
  ) %>%
  pivot_longer(n:mean, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = as.factor(MoonBinNumeric), y = Metric, label = Value)) +
    geom_text(size = 2) +
    facet_wrap(~ActivityPhase) +
    theme_void() +
    theme(axis.text.y = element_text())
p2 <- ggplot(dat_byphase, aes(x = as.factor(MoonBinNumeric), y = meanActivity)) +
  geom_boxplot(alpha = 0.2) +
  theme_minimal() +
  facet_wrap(~ ActivityPhase) +
  xlab("") +
  xlab("MoonBinNumeric")
ggarrange(p2, p1, nrow = 2, heights = c(0.8, 0.2), align = "hv")

# Visualize how the mean activity depends on moonlight intensity
ggplot(dat_byphase, aes(x = maxMoonlightIntensity, y = meanActivity)) +
  geom_point(alpha = 0.2) +
  theme_minimal() +
  facet_wrap(~ ActivityPhase) +
  stat_smooth()

# Let's try to split the data into dark / illuminated nights based on the
# quantiles and see whether activity changes depending on the illumination of
# the night
dat %>%
  select(minMoonlightIntensity, meanMoonlightIntensity, maxMoonlightIntensity) %>%
  apply(2, summary) %>%
  round(6)

# Let's split the data into "dark" and "illuminated" nights
quantiles <- quantile(dat$maxMoonlightIntensity, c(0.20, 0.80))
dat <- mutate(dat, NightIllumination = case_when(
    maxMoonlightIntensity < quantiles[1] ~ "Dark"
  , maxMoonlightIntensity > quantiles[2] ~ "Illuminated"
  , TRUE ~ "Ignore"
))
prop.table(table(dat$NightIllumination))

# Compute averages again
means <- dat %>%
  subset(NightIllumination != "Ignore") %>%
  group_by(Time, NightIllumination) %>%
  summarize(
      Mean    = mean(ActX)
    , SD      = sd(ActX)
    , .groups = "drop"
  ) %>%
  ungroup()

# Visualize them above each other
ggplot(means, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD, col = NightIllumination, fill = NightIllumination)) +
  geom_ribbon(alpha = 0.1, lwd = 0.2) +
  geom_line() +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  theme(axis.text.x = element_text(angle = 45)) +
  facet_wrap(~ NightIllumination, nrow = 2) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("cornflowerblue", "orange")) +
  scale_fill_manual(values = c("cornflowerblue", "orange"))

# Visualize them on top of each other
ggplot(means, aes(x = Time, y = Mean, ymin = Mean - SD, ymax = Mean + SD, col = NightIllumination, fill = NightIllumination)) +
  geom_ribbon(alpha = 0.1, lwd = 0.4) +
  geom_line() +
  theme_minimal() +
  scale_x_time(breaks = scales::breaks_width("2 hours")) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("cornflowerblue", "orange")) +
  scale_fill_manual(values = c("cornflowerblue", "orange"))

################################################################################
#### Some Models
################################################################################
# Let's run some simple linear models to see the relationship
summary(lm(meanActivity ~ maxMoonlightIntensity, data = subset(dat_byphase, ActivityPhase == "Morning")))
summary(lm(meanActivity ~ maxMoonlightIntensity, data = subset(dat_byphase, ActivityPhase == "Evening")))
summary(lm(meanActivity ~ maxMoonlightIntensity, data = subset(dat_byphase, ActivityPhase == "Night")))

# Let's be a bit more elaborate and run the regression for each individual
# separately
mods <- dat_byphase %>%
  nest(Data = -c(ActivityPhase, DogID)) %>%
  mutate(Model = map(Data, function(x) {
    lm(meanActivity ~ maxMoonlightIntensity, data = x)
  })) %>%
  mutate(Coefs = map(Model, function(x) {
    tidy(x)
  }))

# Compute mean and se of the model estimates
coefs_means <- mods %>%
  select(-c(Data, Model)) %>%
  unnest(Coefs) %>%
  group_by(ActivityPhase, term) %>%
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
  ggplot(aes(x = ActivityPhase, y = estimate, col = DogID)) +
    geom_point(
        data        = coefs_means
      , mapping     = aes(x = ActivityPhase, y = MeanEstimate)
      , inherit.aes = F
      , size        = 4
    ) +
    geom_errorbar(
        data        = coefs_means
      , mapping     = aes(x = ActivityPhase, ymin = LWR, ymax = UPR)
      , inherit.aes = F
      , width       = 0.3
      , lwd         = 1.5
    ) +
    # geom_violin(col = "black") + # Could also use a violin plot
    geom_jitter(width = 0.1) +
    facet_wrap(~ term, scales = "free") +
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
