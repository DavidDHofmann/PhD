################################################################################
#### Prepare Activity Data by Time of the Day
################################################################################
# Description: Aggregate dctivity data by time of the day. Remember that we
# prepared nightly statistics for the previous night if a fix was taken before
# 12:00, otherwise we computed statistics for the following night

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(hms)         # To handle times
library(pbmcapply)   # To run stuff in parallel
library(gggibbous)   # For plotting moon phases
library(ggpubr)      # To arrange multiple plots

# Load cleaned activity data (note that the timestamps are all in UTC)
dat <- read_csv("03_Data/02_CleanData/ActivityDataCovariates.csv")

# Let's also generate a column indicating the time of the day
dat$Time <- as_hms(dat$Timestamp)

# Nest the data by dog, collar, and date
dat_nested <- dat %>% nest(Data = -c(DogID, CollarID, State, Date))

# Check the number of datapoints per day
dat_nested$N <- sapply(dat_nested$Data, function(x) {nrow(x)})

# Let's only keep days where we have at least 280 datapoints
dat_nested <- subset(dat_nested, N > 280)

# Can do some plots using the function below (note that a darker moon means that
# it is more illuminated)
plotActivityMoon <- function(index) {
  subdat <- dat_nested %>%
    slice(index) %>%
    unnest(Data)
  p1 <- ggplot(subdat, aes(x = Timestamp, y = MoonAngle)) +
      geom_hline(yintercept = 0, col = "gray", lty = 2) +
      geom_line() +
      geom_moon(data = subdat[10, ], aes(x = Timestamp, y = 0, ratio = 1), col = "black", fill = "white") +
      geom_moon(data = subdat[10, ], aes(x = Timestamp, y = 0, ratio = MoonPhase), fill = "black") +
      theme_minimal() +
      theme(panel.grid.minor = element_blank()) +
      ylim(c(-90, 90))
  p2 <- ggplot(subdat, aes(x = Timestamp, y = Act, col = ToD)) +
      geom_point() +
      geom_vline(aes(xintercept = max(Sunset))) +
      geom_vline(aes(xintercept = max(Sunset) - hours(2)), lty = 2) +
      geom_vline(aes(xintercept = max(Sunset) + hours(2)), lty = 2) +
      geom_vline(aes(xintercept = max(Sunrise))) +
      geom_vline(aes(xintercept = max(Sunrise) - hours(2)), lty = 2) +
      geom_vline(aes(xintercept = max(Sunrise) + hours(2)), lty = 2) +
      theme_minimal() +
      scale_color_viridis_d() +
      theme(legend.position = c(0.1, 0.8)) +
      xlab("") +
      ylim(c(0, 510))
  ggarrange(p2, p1, nrow = 2, align = "hv", heights = c(0.8, 0.2))
}

# Try it
plotActivityMoon(index = 1)

# Function to summarize data between two dates
summarizeDat <- function(data, dogname, collar, state, timestamp_from, timestamp_to) {

  # Subset to data of interest
  subdat <- subset(data,
    DogID     == dogname & CollarID == collar &
    State     == State &
    Timestamp >= timestamp_from &
    Timestamp < timestamp_to # No "equal" to avoid double counting
  )

  # Count number of fixes
  NFixes <- nrow(subdat)

  # Summarize values
  subdat <- summarize(subdat
    , meanAct                = mean(Act)
    , SDAct                  = sd(Act)
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
    , maxMoonTime            = max(maxMoonTime)
    , maxMoonDelay           = max(maxMoonDelay)
  )

  # Assign number of fixes
  subdat$NFixes <- NFixes

  # Also compute average activity in the previous x hours
  for (i in c(6, 12, 18, 24)) {
    meanAct <- subset(data
      , DogID     == dogname &
        CollarID  == collar &
        Timestamp >= timestamp_from - hours(i + 2) &
        Timestamp < timestamp_from - hours(2)
    ) %>% pull(Act) %>% mean()
    subdat <- cbind(subdat, meanAct)
    names(subdat)[ncol(subdat)] <- paste0("meanAct", i)
  }

  # Return it
  return(subdat)

}

# Summarize data by time of the day
dat_nested$Data <- pbmclapply(
    X                  = 1:nrow(dat_nested)
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {

  # Extract relevant info
  subdata      <- dat_nested$Data[[x]]
  dogname      <- dat_nested$DogID[[x]]
  collar       <- dat_nested$CollarID[[x]]
  state        <- dat_nested$State[[x]]
  date         <- unique(as_date(subdata$Timestamp))
  sunrise      <- min(subdata$Sunrise[as_date(subdata$Sunrise) == date])
  sunset       <- max(subdata$Sunset[as_date(subdata$Sunset) == date])
  sunrise_next <- min(subdata$SunriseNext[as_date(subdata$SunriseNext) == date + days(1)])

  # Loop through the three times of interest
  summarized <- lapply(c("Morning", "Evening", "Night"), function(i) {

    # Determine cutoff timestamps
    if (i == "Morning") {
      time_from <- sunrise - hours(2)
      time_to   <- sunrise + hours(2)
    } else if (i == "Evening") {
      time_from <- sunset - hours(2)
      time_to   <- sunset + hours(2)
    } else {
      time_from <- sunset + hours(2)
      time_to   <- sunrise_next - hours(2)
    }

    # Summarize data falling within those timestamps
    summarized <- summarizeDat(dat, dogname, collar, state
      , timestamp_from = time_from
      , timestamp_to   = time_to
    )

    # Assign some more information
    summarized <- cbind(
        ActivityStart = time_from
      , ActivityEnd   = time_to
      , ToD           = i
      , summarized
    )

    # Return the summarized values
    return(summarized)
  })

  # Bind
  summarized <- do.call(rbind, summarized)

  # Return everything
  return(summarized)
})

# Unnest again
dat_tod <- unnest(dat_nested, Data)
dat_tod$N <- NULL

# For morning and evening bursts we also need to compute the average cloud cover
# of he preceeding/following night (although cloud cover is anyways not that
# greatly resolved)
meanCloudCoverNight <- sapply(1:nrow(dat_tod), function(x) {
  if (dat_tod$ToD[x] == "Night") {
    meanCloudCoverNight <- dat_tod$meanCloudCover[x]
  } else if (dat_tod$ToD[x] == "Morning") {
    subdat <- subset(dat_tod,
      DogID    == dat_tod$DogID[x] &
      CollarID == dat_tod$CollarID[x] &
      State    == dat_tod$State[x] &
      Date     == dat_tod$Date[x] - days(1) &
      ToD      == "Night"
    )
    meanCloudCoverNight <- ifelse(nrow(subdat) == 0, NA, mean(subdat$meanCloudCover))
  } else {
    subdat <- subset(dat_tod,
      DogID    == dat_tod$DogID[x] &
      CollarID == dat_tod$CollarID[x] &
      State    == dat_tod$State[x] &
      Date     == dat_tod$Date[x] &
      ToD      == "Night"
    )
    meanCloudCoverNight <- ifelse(nrow(subdat) == 0, NA, mean(subdat$meanCloudCover))
  }
})
dat_tod$meanCloudCoverNight <- meanCloudCoverNight

# Store the data
write_csv(dat_tod, "03_Data/02_CleanData/ActivityDataCovariatesAggregated.csv")
