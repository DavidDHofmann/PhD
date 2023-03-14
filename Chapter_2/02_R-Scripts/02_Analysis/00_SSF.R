################################################################################
#### Step Selection Function - Generation of Random Steps
################################################################################
# Description: In this script we coerce our gps data to steps and generate
# random steps for the step selection analysis.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load packages
library(tidyverse)    # For data wrangling
library(lubridate)    # To handle dates
library(pbmcapply)    # For multicore abilities
library(fitdistrplus) # To fit distributions
library(terra)        # To handle spatial data
library(ggpubr)       # To combine multiple plots

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Define the number of cores to use
cores <- detectCores() - 1

################################################################################
#### Computing Step Metrics
################################################################################
# Load the gps data of dispersers and create rounded timestamps. Then, nest the
# data by individual
data <- "03_Data/02_CleanData/Dispersers.csv" %>%
  read_csv(show_col_types = F) %>%
  subset(State == "Disperser") %>%
  mutate(TimestampRounded = round_date(Timestamp, "1 hour")) %>%
  rename(ID = DogName) %>%
  nest(GPS = -ID)

# For some individuals we have only few observations during dispersal
print(data, n = 30)

# For some individuals the temporal resolution of the data is much higher than
# for others. Hence, we will want to equalize the fixrate scheme for all
# individuals  Go through the data of all individuals. Thus, let's resample the
# fixes, remove fixes not collected on the anticipated schedule, compute bursts,
# and reproject coordinates
cat("Computing bursts from GPS data...\n")
data$Bursts <- pbmclapply(data$GPS, ignore.interactive = T, mc.cores = cores, function(x) {
  res <- resampleFixes(x, hours = 2, start = 1, tol = 0.5)
  if (length(res) == 1) {
    return(NA)
  }
  sub <- subset(res, hour(TimestampRounded) %in% c(3, 7, 15, 19, 23))
  if (length(sub) == 1) {
    return(NA)
  }
  bur <- computeBursts(sub)
  bur[, c("x", "y")] <- reprojectCoords(cbind(bur$x, bur$y)
    , from = "+init=epsg:4326"
    , to   = "+init=epsg:32734"
  )
  return(bur)
})

# Check for duplicated timestamps (there should be none)
data %>%
  dplyr::select(ID, Bursts) %>%
  unnest(Bursts) %>%
  group_by(ID, TimestampRounded) %>%
  filter(n() > 1)

# We can only work with bursts that contain at least 3 fixes. Let's compute step
# metrics for those steps.
steps <- data %>%
  dplyr::select(ID, Bursts) %>%
  unnest(Bursts) %>%
  nest(data = -c(ID, burst_id)) %>%
  mutate(Nrow = map_dbl(data, nrow)) %>%
  subset(Nrow > 3) %>%
  mutate(Steps = map(data, stepMetrics)) %>%
  dplyr::select(ID, BurstID = burst_id, Steps) %>%
  mutate(BurstID = 1:n()) %>%
  unnest(Steps)

# Remove steps where the relative turning angle is NA
steps <- subset(steps, !is.na(relta))

# We can't work with 0 step lengths so we'll force each step to at least 1m
steps$sl <- ifelse(steps$sl == 0, 1, steps$sl)

# Assign unique step ids
steps$step_id <- 1:nrow(steps)

# Indicate if a fix was taken during time of activity or inactivity (we define
# inactivity as anything between 07:00 and 15:00)
steps$inactive <- hour(steps$TimestampRounded) == 7

# Remove unnecessary columns
steps <- dplyr::select(steps, -c(DOP, State, Source))

################################################################################
#### Some Summary Statistics
################################################################################
# When are the fixes collected?
steps %>%
  count(hour(TimestampRounded)) %>%
  plot(type = "l")
steps %>%
  count(month(TimestampRounded)) %>%
  plot(type = "l")
steps %>%
  count(year(TimestampRounded)) %>%
  plot(type = "l")

# Visualize the steps
ggplot(steps, aes(x = x, y = y, col = as.factor(BurstID), group = BurstID)) +
  geom_path(size = 0.1) +
  geom_point(size = 0.5) +
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "none")

# Visualize step lengths and turning angles during activity and inactivity
p1 <- ggplot(steps, aes(x = sl, col = inactive, fill = inactive)) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  ggtitle("Step Lengths") +
  xlab("Step Length (m)") +
  ylab("Density") +
  scale_fill_manual(values = c("orange", "cornflowerblue")) +
  scale_color_manual(values = c("orange", "cornflowerblue"))
p2 <- ggplot(steps, aes(x = relta, col = inactive, fill = inactive)) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  ggtitle("Turning Angles") +
  xlab("Turning Angle (rad)") +
  ylab("Density") +
  scale_fill_manual(values = c("orange", "cornflowerblue")) +
  scale_color_manual(values = c("orange", "cornflowerblue"))

# Visualize step lengths and turning angles for different sexes
p3 <- ggplot(steps, aes(x = sl, col = Sex, fill = Sex)) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  ggtitle("Step Lengths") +
  xlab("Step Length (m)") +
  ylab("Density") +
  scale_fill_manual(values = c("orange", "cornflowerblue")) +
  scale_color_manual(values = c("orange", "cornflowerblue"))
p4 <- ggplot(steps, aes(x = relta, col = Sex, fill = Sex)) +
  geom_density(alpha = 0.2) +
  theme_minimal() +
  ggtitle("Turning Angles") +
  xlab("Turning Angle (rad)") +
  ylab("Density") +
  scale_fill_manual(values = c("orange", "cornflowerblue")) +
  scale_color_manual(values = c("orange", "cornflowerblue"))

# Put plots together
ggarrange(p1, p2, p3, p4)

################################################################################
#### Fitting Step-Length Distribution
################################################################################
# Fit a Gamma distribution to the step length. To avoid scaling issues, we'll
# fit the gamma distribution to step lengths in kilometers
sl <- fitdist(as.numeric(steps$sl), "gamma", method = "mle", lower = 0)
sl <- list(
    shape = sl$estimate[["shape"]]
  , scale = 1 / sl$estimate[["rate"]]
)

# Let's visualize the fit
x <- seq(0, max(steps$sl), length.out = 1000)
y <- dgamma(x, shape = sl$shape, scale = sl$scale)
hist(steps$sl, freq = F, breaks = 100, col = "cornflowerblue", border = "white")
lines(y ~ x, col = "orange", lwd = 3)

################################################################################
#### Generating Random Steps
################################################################################
# Generate random steps
set.seed(12345)
cat("Generating random steps...\n")
ssf <- randomSteps(steps
  , n_rsteps = 100
  , scale    = sl$scale
  , shape    = sl$shape
  , slmax    = max(steps$sl)
)

# I'd like to work exclusively with lon-lat data instead of utm, hence, let's
# reproject start and end coordinates to lon-lat
ssf[, c("x", "y")] <- reprojectCoords(
    xy   = cbind(ssf$x, ssf$y)
  , from = "+init=epsg:32734"
  , to   = "+init=epsg:4326"
)
ssf[, c("x_to", "y_to")] <- reprojectCoords(
    xy   = cbind(ssf$x_to, ssf$y_to)
  , from = "+init=epsg:32734"
  , to   = "+init=epsg:4326"
)

# Visualize the random steps of a single individual
ssf %>%
  subset(ID == "Abel") %>%
  ggplot(aes(x = x, y = y, xend = x_to, yend = y_to, col = as.factor(case))) +
    geom_segment(lwd = 0.1) +
    theme_minimal() +
    coord_sf() +
    scale_color_manual(values = c("orange", "cornflowerblue"))

# Write the step selection data to file
write_csv(ssf, "03_Data/02_CleanData/SSF.csv")

################################################################################
#### Session Information
################################################################################
# Create folder into which the session info goes to
dir.create("02_R-Scripts/99_SessionInformation/02_Analysis", showWarnings = F)

# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/02_Analysis/00_SSF.rds")
cat("Done :)\n")
