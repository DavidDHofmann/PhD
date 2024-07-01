################################################################################
#### Visualization of Moonlight
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)    # For data-wrangling
library(lubridate)    # To handle dates
library(suncalc)      # To identify sunrise and sunset times
library(hms)          # To work with times
library(ggpubr)       # To arrange multiple ggplots
library(gggibbous)    # To plot the moon

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

# Define coordinates we want to get astronomical data for
lon <- 23.5
lat <- -19

################################################################################
#### Moon Illumination at Night
################################################################################
# Define a cutoff (in lux) that we classify as being "dark"
cutoff     <- 0.02
timestamps <- seq(
    from = ymd_hms("2000-06-01 00:00:00")
  , to   = ymd_hms("2000-07-01 00:00:00")
  , by   = "5 mins"
)

# Compute statistics for each timestamp
moon                   <- moonlightIntensity(timestamps, lon = lon, lat = lat, e = 0.21)
moon$moonlightModelLux <- moon$moonlightModel * 0.32
moon$LightType         <- ifelse(moon$moonlightModelLux < cutoff, "Dark", "Bright")

# Reload aggregate statistics
moon_aggr           <- read_rds("03_Data/02_CleanData/Moonlight.rds")
moon_aggr$LightType <- ifelse(moon_aggr$Night > 0.5 & moon_aggr$meanMoonlightLux < cutoff, "Dark", "Bright")

# Subset to nights which would be considered "bright"
moon_aggr <- subset(moon_aggr
  , LightType == "Bright"
  & Night > 0.5 & date(Timestamp) %in% date(moon$date)
)

# Visualize moon illumination over time
p1 <- ggplot(moon, aes(x = date, y = moonAltDegrees)) +
  geom_line(aes(color = moonlightModelLux)) +
  geom_rug(data = moon_aggr, aes(x = Timestamp), inherit.aes = F) +
  geom_moon(data = moon[seq(1, nrow(moon), length.out = 20), ], aes(x = date, y = 0, ratio = 1), col = "black", fill = "black", size = 4) +
  geom_moon(data = moon[seq(1, nrow(moon), length.out = 20), ], aes(x = date, y = 0, ratio = moonPhase), fill = "white", size = 4) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = seq(-90, 90, by = 30)) +
  scale_color_distiller(palette = "Spectral", name = "Illumination (lx)") +
  xlab("") +
  ylab("Moon Altitude °")

# Visualize light-type over time
p2 <- ggplot(moon, aes(x = date, y = moonAltDegrees, group = 1)) +
  geom_line(aes(color = LightType)) +
  geom_rug(data = moon_aggr, aes(x = Timestamp), inherit.aes = F) +
  geom_moon(data = moon[seq(1, nrow(moon), length.out = 20), ], aes(x = date, y = 0, ratio = 1), col = "black", fill = "black", size = 4) +
  geom_moon(data = moon[seq(1, nrow(moon), length.out = 20), ], aes(x = date, y = 0, ratio = moonPhase), fill = "white", size = 4) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = seq(-90, 90, by = 30)) +
  scale_color_manual(values = c("orange", "cornflowerblue"), name = "Light Type") +
  xlab("") +
  ylab("Moon Altitude °")

# Store plots
p <- ggarrange(p1, p2, nrow = 2, align = "hv")

# Arrange the plots and store
ggsave("04_Manuscript/Figures/LightTypes.png"
  , plot   = p
  , device = png
  , width  = 8
  , height = 5
  , bg     = "white"
)

################################################################################
#### Check Step-Length under Different Moonlight Conditions (REMOVE LATER)
################################################################################
# Load the step-selection data and cut moonlight into groups
dat <- "03_Data/02_CleanData/SSFExtracted.rds" %>%
  read_rds() %>%
  subset(case == 1 & Night > 0.5) %>%
  mutate(MoonlightGroup = cut(meanMoonlight, breaks = 40))

# Function to obtain difference in step-lengths between bright and dark nights
deltaSL <- function(threshold) {
  dat$Bright <- dat$meanMoonlight > threshold
  table(dat$Bright)
  stat <- dat %>%
    group_by(Bright) %>%
    summarize(Count = n(), meanSL = mean(sl))
  return(stat)
}

# Apply it using different thresholds for bright vs. dark
act <- lapply(seq(0.005, 0.2, by = 0.005), function(x) {
  del <- deltaSL(x)
  del$Threshold <- x
  return(del)
}) %>% do.call(rbind, .)

# Visualize
ggplot(act, aes(x = Threshold, y = meanSL, col = Bright)) +
  geom_line() +
  theme_awesome()

# Get ratios
act %>%
  dplyr::select(-Count) %>%
  pivot_wider(names_from = Bright, values_from = meanSL) %>%
  mutate(Ratio = `TRUE` / `FALSE`) %>%
  ggplot(aes(x = Threshold, y = Ratio)) +
    geom_line() +
    theme_awesome()

################################################################################
#### Moonlight Intensity
################################################################################
# Span a vector of timestamps for which we want to identify moonlight and sun
# metrics
timestamps <- seq(
    from = ymd_hms("2023-01-01 03:00:00")
  , to   = ymd_hms("2023-03-31 22:00:00")
  , by   = "5 mins"
)

# Compute sunrise and sunset times (this is mainly for visualization purposes)
sun <- getSunlightTimes(
    date = unique(as_date(timestamps))
  , lon = lon
  , lat = lat
  ) %>%
  select(date, nightEnd, nauticalDawn, dawn, sunrise, sunset, dusk, nauticalDusk, night) %>%
  mutate(
      sunrise      = as_hms(substr(sunrise, start      = 12, stop = 20))
    , sunset       = as_hms(substr(sunset, start       = 12, stop = 20))
    , dusk         = as_hms(substr(dusk, start         = 12, stop = 20))
    , dawn         = as_hms(substr(dawn, start         = 12, stop = 20))
    , nauticalDusk = as_hms(substr(nauticalDusk, start = 12, stop = 20))
    , nauticalDawn = as_hms(substr(nauticalDawn, start = 12, stop = 20))
    , nightEnd     = as_hms(substr(nightEnd, start     = 12, stop = 20))
    , night        = as_hms(substr(night, start        = 12, stop = 20))
  ) %>%
  pivot_longer(nightEnd:night, names_to = "phase", values_to = "time")

# Compute moonlight statistics for all datapoints
moon <- moonlightIntensity(
    date = timestamps
  , lon  = 23.5
  , lat  = -19
  , e    = 0.21
  ) %>% mutate(
    time = as_hms(substr(date, start = 12, stop = 20))
  , date = as_date(date)
)

# Compute moonlight in lux
moon$moonlightModelLux   <- moon$moonlightModel * 0.32
moon$moonlightModel24Lux <- moon$moonlightModel24 * 0.32

# Visualize moon-light
p <- ggplot() +
  geom_raster(
      data    = subset(moon, year(date) == 2023)
    , mapping = aes(x = date, y = time, fill = moonlightModel24Lux)
  ) +
  geom_path(
      data    = subset(sun, phase %in% c("sunrise", "sunset") & year(date) == 2023)
    , mapping = aes(x = date, y = time, col = phase)
  ) +
  scale_color_manual(values = c("orange", "cornflowerblue"), name = "Event") +
  scale_fill_viridis_c(name = "Moonlight-Intensity (in lux)") +
  scale_x_date(date_labels = "%Y-%m-%d") +
  guides(
      fill  = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(7, "cm"))
    , color = guide_legend(title.position   = "top", title.hjust = 0.5)
  ) +
  theme_awesome() +
  facet_wrap(~ "Moonlight Intensity") +
  xlab("Date") +
  ylab("Time") +
  theme(
      axis.title.y = element_text(angle = 90)
    , axis.text.x  = element_text(angle = 40)
  )

# Arrange the plots and store
ggsave("04_Manuscript/Figures/Moonlight.png"
  , plot   = p
  , device = png
  , width  = 5
  , height = 5
  , scale  = 1.2
  , bg     = "white"
)
