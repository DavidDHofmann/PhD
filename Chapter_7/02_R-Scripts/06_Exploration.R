################################################################################
#### Exploratory Analysis
################################################################################
# Testing

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
library(pracma)      # To identify peaks
library(broom)       # To clean model summary
library(jagsUI)      # To run models in bayesian framework
library(runner)      # To apply moving windows

# Reload cleaned activity data (note that the timestamps are all in UTC)
dat <- read_csv("03_Data/02_CleanData/ActivityDataCovariates.csv")
print(names(dat))

################################################################################
#### TESTING
################################################################################
# Subset to one of the individuals
unique(dat$DogID)
sub <- subset(dat, DogID == "Calvin")
table(sub$CollarID)

# Let's compute average activity for different moving windows
windows <- c("60 minutes", "120 minutes")
for (i in windows) {
  ran <- runner(sub
    , k   = duration(i)
    , idx = sub$Timestamp
    , f   = function(x) {mean(x$ActX)}
  )
  sub <- cbind(sub, ran)
  names(sub)[ncol(sub)] <- paste0("ActX_", gsub(i, pattern = " ", replacement = ""))
}

# Active or inactive
# Transition probabilities
x<-read.table("http://www.rolandlangrock.com//OF.dat")





# Pivot and visualize
sub %>%
  pivot_longer(ActX_60minutes:ActX_120minutes, names_to = "Window", values_to = "ActXWindow") %>%
  drop_na() %>%
  subset(Date == min(Date) + days(6)) %>%
  ggplot(aes(x = hms::as_hms(Timestamp), y = ActXWindow)) +
    geom_line(lwd = 1) +
    geom_line(aes(y = ActX)) +
    facet_wrap(~ Window, ncol = 1)

sub %>%
  drop_na() %>%
  subset(Date == min(Date) + days(28)) %>%
  ggplot(aes(x = hms::as_hms(Timestamp), y = Running)) +
    geom_line(lwd = 2) +
    geom_line(aes(y = ActX))

################################################################################
#### NEED TO APPLY SOME FILTERS
################################################################################
# Look at residents only?
# What to do with individuals that were together?
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
  geom_col(fill = "cornflowerblue", col = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

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

# We will clearly have to aggregate data at some point! Otherwise the
# distributions will be very nasty to work with. Anyways, let's continue. We
# deem an animal active whenever its activity is above 20
dat$Active <- dat$ActX > 20

# How many points during activity and inactivity?
table(dat$Active) / 1e6

# Make sure the data is arranged properly
dat <- arrange(dat, DogID, CollarID, Timestamp)

# Compute temporal lag between consecutive recordings and assess whenever the
# temporal lag is above 5 mins, thus initiating a new burst. Also assess the
# number of recordings in each burst, as well as the start and enddate of it
dat <- dat %>%
  nest(Data = -c(DogID, CollarID)) %>%
  mutate(Data = map(Data, function(x) {
    x$Timelag  <- difftime(x$Timestamp, lag(x$Timestamp), unit = "mins")
    x$Newburst <- as.numeric(x$Timelag) > 5
    x$Newburst <- ifelse(is.na(x$Newburst), F, x$Newburst)
    x$BurstID  <- cumsum(x$Newburst)
    x$Newburst <- NULL
    return(x)
  })) %>%
  unnest(Data) %>%
  nest(Data = -c(DogID, CollarID, BurstID)) %>%
  mutate(
      NRecords = sapply(Data, nrow)
    , Start    = lapply(Data, function(x) {min(x$Timestamp)}) %>% do.call(c, .)
    , End      = lapply(Data, function(x) {max(x$Timestamp)}) %>% do.call(c, .)
  )

# Let's take a look at the data
print(dat)

# Make sure that each burst has a decent amount of recordings/fixes
summary(dat$NRecords)
hist(dat$NRecords, col = "cornflowerblue", border = "white")

# We only want to keep bursts were we have at least one day of data!
NRowOld <- sum(dat$NRecords)
dat <- subset(dat, NRecords >= 24 * 60 / 5)
NRowNew <- sum(dat$NRecords)
NRowOld - NRowNew
rm(NRowOld, NRowNew)

# Take a look at the data again
print(dat)

# Unnest the data for further calculations
dat <- dat %>%
  select(DogID, CollarID, BurstID, Data) %>%
  unnest(Data)

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
  group_by(ToD) %>%
  summarize(MeanActX = mean(ActX), MedianActX = median(ActX)) %>%
  pivot_longer(MeanActX:MedianActX, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = ToD, y = Value, fill = as.factor(Metric))) +
    geom_col(position = position_dodge()) +
    theme_minimal() +
    scale_fill_manual(values = c("cornflowerblue", "orange"), name = "Metric")

# How does activity change during different times of the day depending on the
# season?
dat %>%
  group_by(ToD, Season) %>%
  summarize(MeanActX = mean(ActX)) %>%
  ggplot(aes(x = ToD, y = MeanActX, fill = as.factor(Season))) +
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
ggplot(dat, aes(x = Time, y = update(Date, year = 2000), z = ActX)) +
  stat_summary_2d(fun = "mean", bins = 100) +
  scale_fill_viridis_c(option = "magma", name = "Activity") +
  theme_minimal() +
  scale_y_date(date_breaks = "3 months") +
  ylab("FAKE YEAR!!!")

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

# We can try to identify the peaks and valleys in this plot
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
  scale_fill_manual(values = c("cornflowerblue", "orange", "gray"))

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
  facet_wrap(~ Month, ncol = 2, dir = "v") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = "bottom") +
  scale_fill_manual(values = c("cornflowerblue", "orange", "gray"))

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
  scale_fill_manual(values = c("cornflowerblue", "orange", "gray"))

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

# dat %>%
#   subset(ActX > 0) %>%
#   group_by(ActivityPhase) %>%
#   summarize(
#       mean = quantile(ActX, 0.25)
#     , min  = min(ActX)
#     , q1   = quantile(ActX, 0.25)
#     , q2   = quantile(ActX, 0.50)
#     , q3   = quantile(ActX, 0.75)
#     , max  = max(ActX)
#   )
# ggplot(subset(dat, ActX > 0), aes(x = ActivityPhase, y = ActX)) +
#   geom_boxplot()

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
    geom_violin(alpha = 0.5) +
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
    geom_violin(alpha = 0.5) +
    theme_minimal() +
    scale_fill_manual(values = c("cornflowerblue", "orange", "gray"))

# How well is the correlation between mean and total activity?
cor(dat_byphase$meanActivity, dat_byphase$totalActivity)

# Let's bin the nights by moon illumination
dat_byphase <- mutate(dat_byphase, MoonBin = cut(maxMoonlightIntensity, breaks = 5))
dat_byphase$MoonBinNumeric <- as.numeric(dat_byphase$MoonBin)

# Visualize how the mean activity depends on moonlight intensity
ggplot(dat_byphase, aes(x = as.factor(MoonBinNumeric), y = meanActivity)) +
  geom_boxplot(alpha = 0.2) +
  theme_minimal() +
  facet_wrap(~ ActivityPhase)

# Visualize how the mean activity depends on moonlight intensity
ggplot(dat_byphase, aes(x = maxMoonlightIntensity, y = meanActivity, fill = ActivityPhase, col = ActivityPhase)) +
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
quantiles <- quantile(dat$maxMoonlightIntensity, c(0.10, 0.90))
dat <- mutate(dat, NightIllumination = case_when(
    maxMoonlightIntensity < quantiles[1] ~ "Dark"
  , maxMoonlightIntensity > quantiles[2] ~ "Illuminated"
  , TRUE ~ "Ignore"
))

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

# For now, subset the data
dat_backup <- dat
dat <- subset(dat, DogID %in% c("Abel", "Oolong"))

# We now want to assess at what time of the day the dogs usually become active.
# Hence, we want to nest the data by dog/collar first
dat <- nest(dat, Data = -c("DogID", "CollarID", "BurstID"))
print(dat)

# Based on this, we now create a 30 min (i.e. 6 activity periods) moving window
# to find the time at which they start becoming active
periods <- 6
dat <- mutate(dat, Data = map(Data, function(x) {
  x$Status <- NA
  for (i in 1:nrow(x)) {
    if (i < periods | as.numeric(x$Lag[i]) > 5 | is.na(x$Lag[i])) {
      x$Status[i] <- NA
    } else {
      if (all(x$Active[(i - (periods - 1)):i])) {
        x$Status[i] <-"Moving"
      } else if (all(!x$Active[(i - (periods - 1)):i])) {
        x$Status[i] <- "Resting"
      } else {
        x$Status[i] <- x$Status[i-1]
      }
    }
  }
  return(x)
}))

# Visualize status
dat %>%
  unnest(Data) %>%
  subset(!is.na(Status)) %>%
  ggplot(aes(x = Time, y = Status, col = Date)) +
    geom_jitter(alpha = 0.5, pch = 20) +
    theme_minimal() +
    scale_color_viridis_c(option = "magma")

# Identify switchpoints at which individuals go from moving to resting or from
# resting to moving. Also determine the direction of the switch
dat <- mutate(dat, Data = map(Data, function(x) {
  x$Switch <- lag(x$Status) != x$Status
  x$SwitchTo <- NA
  indices <- which(x$Switch)
  for (i in indices) {
    if (x$Status[i-1] == "Resting" & x$Status[i] == "Moving") {
      x$SwitchTo[i] <- "ToMoving"
    } else if (x$Status[i-1] == "Moving" & x$Status[i] == "Resting") {
      x$SwitchTo[i] <- "ToResting"
    } else {
      x$SwitchTo[i] <- NA
    }
  }
  return(x)
}))

# Visualize the switches
dat %>%
  unnest(Data) %>%
  subset(!is.na(SwitchTo)) %>%
  ggplot(aes(x = Time, y = SwitchTo)) +
    geom_jitter(pch = 20) +
    theme_minimal()

# Let's figure out at what time of the day the individuals start moving for the
# first time in the evening (afternoon)
first <- mutate(dat, FirstMovement = map(Data, function(x) {

  # Find all days through which we need to loop
  days <- unique(x$Date)
  firstmoves <- lapply(days, function(y) {
    # first <- x[x$Date == y & x$Time > phase$Endtime[phase$Burst == "Morning"] & x$Status == "Moving", ]
    first <- x[x$Date == y & x$Time > as_hms("12:00:00") & x$Status == "Moving", ]
    first <- first[1, ]
    return(first)
  }) %>% do.call(rbind, .)

  # Return the first moves
  firstmoves <- tibble(Day = days, FirstMovement = firstmoves$Timestamp)
  return(firstmoves)
}))

# Let's visualize them
first %>%
  select(DogID, FirstMovement) %>%
  unnest(FirstMovement) %>%
  subset(!is.na(FirstMovement)) %>%
  ggplot(aes(x = as_hms(FirstMovement))) +
    geom_histogram(bins = 30) +
    facet_wrap(~ DogID) +
    theme_minimal()

# Unnest the data again
dat <- unnest(dat, Data)

############ CONTINUE HERE!!!!

# Join this data with nightly statistics
stats <- dat %>%
  subset(Hour > 14) %>%
  select(DogID, Date, maxMoonlightIntensity) %>%
  mutate(maxMoonlightIntensity = round(maxMoonlightIntensity, 5)) %>%
  distinct()
test <- first %>%
  select(DogID, FirstMovement) %>%
  unnest(FirstMovement) %>%
  subset(!is.na(FirstMovement)) %>%
  subset(hour(FirstMovement) > 14) %>%
  left_join(., stats, by = c("DogID", "Day" = "Date"))
test$FirstMovement <- as.numeric(minutes(as_hms(test$FirstMovement)))
plot(as_hms(test$FirstMovement) ~ test$maxMoonlightIntensity)
abline(mod)
mod <- lm(FirstMovement ~ maxMoonlightIntensity, data = test)
summary(mod)

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
    , LWR          = MeanEstimate - qt(1 - 0.05 / 2, n() - 1) * SE
    , UPR          = MeanEstimate + qt(1 - 0.05 / 2, n() - 1) * SE
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

# For comparison, first run a regular model
mod <- lm(meanActivity ~ maxMoonlightIntensity, data = sub)

# Bundle data
dat_jags <- list(
    meanActivity          = sub$meanActivity
  , maxMoonlightIntensity = sub$maxMoonlightIntensity
  , n                     = nrow(sub)
)

# Write JAGS model file
file <- tempfile(fileext = ".txt")
cat(file = file, "model {
  # Priors
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  sd    ~ dunif(0, 100)
  variance  <- sd * sd
  precision <- 1 / variance

  # Likelihood
  for (i in 1:n){
    meanActivity[i] ~ dnorm(mu[i], precision)
    mu[i] <- beta0 + beta1 * maxMoonlightIntensity[i]
  }
}")

# Function to sample initial values
inits <- function() {
  list(
      beta0 = rnorm(1)
    , beta1 = rnorm(1)
    , sd    = rlnorm(1)
  )
}

# Define parameters to be monitored (i.e. estimated)
params <- c("beta0", "beta1", "sd")

# MCMC settings (usually defined using trial and error)
na <- 1000        # Number of iterations in the adaptive phase
ni <- 3000        # Number of draws from the posterior (in each chain)
nb <- 1000        # Number of draws to discard as burn-in
nc <- 5           # Number of chains
nt <- 1           # Thinning rate (nt = 1 means we do not thin)

# Run the model
mod_jags <- jags(
    data               = dat_jags
  , inits              = inits
  , parameters.to.save = params
  , model.file         = file
  , n.iter             = ni
  , n.burnin           = nb
  , n.chains           = nc
  , n.thin             = nt
  , n.adapt            = na
  , parallel           = T
)

# Show traceplots
par(mfrow = c(2, 2))
jagsUI::traceplot(mod_jags)

# Summary of output, rounded to 3 digits
print(mod_jags, 3)

# Let's look at the estimates and compare them to the frequentist approach
cbind(
    Frequentist = coef(mod)
  , Bayesian = unlist(mod_jags$mean)[1:2]
)
