################################################################################
#### Plot of the Activity Data Against the Moonangle
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(gggibbous)   # For plotting moon phases
library(ggpubr)      # To arrange multiple plots
library(hms)         # To handle times

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
  p1 <- ggplot(subdat, aes(x = Timestamp, y = ActX, col = ToD)) +
      geom_point(alpha = 0.5) +
      geom_vline(aes(xintercept = max(Sunset))) +
      geom_vline(aes(xintercept = max(Sunset) - hours(2)), lty = 2) +
      geom_vline(aes(xintercept = max(Sunset) + hours(2)), lty = 2) +
      theme_minimal() +
      scale_color_viridis_d() +
      theme(legend.position = "bottom") +
      xlab("") +
      ylim(c(0, 255))
  p2 <- ggplot(subdat, aes(x = Timestamp, y = MoonAngle)) +
      geom_hline(yintercept = 0, col = "gray", lty = 2) +
      geom_line() +
      geom_moon(data = subdat[10, ], aes(x = Timestamp, y = 0, ratio = 1), col = "black", fill = "white") +
      geom_moon(data = subdat[10, ], aes(x = Timestamp, y = 0, ratio = MoonPhase), fill = "black") +
      theme_minimal() +
      theme(panel.grid.minor = element_blank()) +
      ylim(c(-90, 90))
  p3 <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p <- ggarrange(p1, p2, nrow = 2, align = "hv", heights = c(0.8, 0.2))
  p <- ggarrange(p3, p, nrow = 2, heights = c(0.1, 0.9))
  p
}

# Try it
plotActivityMoon(index = 1019)

# Store it
ggsave("04_Manuscript/99_ActivityByMoongangle.png"
  , plot = last_plot()
  , bg   = "white"
)
