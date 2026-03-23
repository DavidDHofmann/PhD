################################################################################
#### A Few General Metrics on the Dispersal Events
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- ifelse(Sys.info()["sysname"] == "Linux"
  , "/media/david/SharedSpace/SwitchDrive/02_Academia/02_PhD/Chapter_3"
  , "D:/SwitchDrive/02_Academia/02_PhD/Chapter_3"
)
setwd(wd)

# Load required packages
library(tidyverse)    # For plotting and data wrangling
library(lubridate)    # To handle dates
library(hms)          # To handle times
library(ggpubr)       # To arrange multiple plots

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load the raw and cleaned gps data (for some metrics we'll use the raw data,
# for others we'll use the cleaned one)
dat_raw   <- read_csv("03_Data/02_CleanData/Dispersers.csv")
dat_clean <- read_rds("03_Data/02_CleanData/SSFExtracted.rds") %>% subset(case == 1)

# Data so we can overlay the seasons
seasons <- tibble(
    Season = c("Wet", "Dry", "Wet")
  , Start  = yday(dmy(c("01.01.2000", "16.04.2000", "16.10.2000")))
  , End    = yday(dmy(c("15.04.2000", "15.10.2000", "31.12.2000")))
)

# Arrange dataframes
dat_raw <- dat_raw %>%
  distinct() %>%
  arrange(desc(Sex), desc(ID)) %>%
  mutate(ID = if_else(Sex == "F"
    , paste0(ID, " \u2640")
    , paste0(ID, " \u2642")
  )) %>%
  mutate(ID = factor(ID, levels = unique(ID))) %>%
  mutate(Timestamp = Timestamp + hours(2))

dat_clean <- dat_clean %>%
  distinct() %>%
  arrange(desc(Sex), desc(ID)) %>%
  mutate(ID = if_else(Sex == "F"
    , paste0(ID, " \u2640")
    , paste0(ID, " \u2642")
  )) %>%
  mutate(ID = factor(ID, levels = unique(ID))) %>%
  mutate(Timestamp = Timestamp + hours(2))

# Create overview plot
p1 <- gpsOverview(dat_raw, id = "ID", timestamp = "Timestamp")
p2 <- gpsOverview(dat_clean, id = "ID", timestamp = "Timestamp")

# Put the plots together and store to file
p  <- ggarrange(p1, p2, ncol = 1, labels = c("a", "b"))
ggsave("04_Manuscript/Figures/GPSCleanup.png"
  , plot   = p
  , width  = 5
  , height = 6
  , bg     = "white"
  , scale  = 1.75
  , device = png
)

################################################################################
#### Collar Duration by Individual and Sex
################################################################################
# Plot data periods
p1 <- ggplot(dat_clean, aes(x = yday(Timestamp), y = ID, color = ID)) +
  geom_rect(data = subset(seasons, Season == "Wet"), aes(xmin = Start, xmax = End, ymin = 0.4, ymax = Inf, group = Season), inherit.aes = F, alpha = 0.1) +
  geom_point(shape = 15, size = 4) +
  scale_color_viridis_d() +
  xlab("Day of the Year") +
  theme_awesome() +
  scale_x_continuous(breaks = seq(50, 350, by = 50)) +
  theme(
      legend.position  = "none"
    , panel.grid.minor = element_blank()
  ) +
  facet_wrap(~ "Dispersal Periods")

################################################################################
#### Number of Fixes by Season
################################################################################
# Calculate the number of fixes per dog
number_fixes <- dat_clean %>%
  count(ID, SeasonClimate, name = "NumberOfFixes")
mean_number_fixes <- mean(number_fixes$NumberOfFixes)

# Prepare a plot
p2 <- ggplot(number_fixes, aes(x = NumberOfFixes, y = ID, fill = ID, color = ID, alpha = SeasonClimate)) +
  geom_col(width = 0.5) +
  geom_vline(xintercept = mean_number_fixes, linetype = "22", col = "gray30") +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  scale_x_continuous(labels = function(x) {format(x, big.mark = ",")}) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Number of Fixes") +
  facet_wrap(~ "Number of Fixes") +
  theme_awesome() +
  theme(
      legend.position  = "none"
    , panel.grid.minor = element_blank()
    , axis.text.y      = element_blank()
    , axis.title.y     = element_blank()
  )

################################################################################
#### Average Daily Distance Dispersed
################################################################################
# Compute the average daily distance dispersed
average_distance <- dat_clean %>%
  mutate(Date = as.Date(Timestamp)) %>%
  group_by(ID, Date) %>%
  summarize(Distance = mean(sl) * 5, .groups = "drop") %>%
  group_by(ID) %>%
  summarize(MeanDistance = mean(Distance), SDDistance = sd(Distance))
mean_average_distance <- mean(average_distance$MeanDistance)

# Visualize it
p3 <- ggplot(average_distance, aes(x = MeanDistance, y = ID, fill = ID, color = ID)) +
  geom_col(alpha = 0.75, width = 0.5) +
  geom_errorbarh(aes(xmin = MeanDistance - SDDistance, xmax = MeanDistance + SDDistance), height = 0) +
  geom_vline(xintercept = mean_average_distance, linetype = "22", col = "gray30") +
  scale_x_continuous(labels = function(x) {format(x / 1000, big.mark = ",")}, breaks = seq(0, 60000, by = 10000)) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  xlab("Average Daily Distance Dispersed (km)") +
  theme_awesome() +
  theme(
      legend.position  = "none"
    , panel.grid.minor = element_blank()
    , axis.text.y      = element_blank()
    , axis.title.y     = element_blank()
  ) +
  facet_wrap(~ "Daily Distance Dispersed")

################################################################################
#### Arrange the Plots
################################################################################
# Arrange plots and store to file
p <- ggarrange(p1, p2, p3, nrow = 1, align = "h", widths = c(4, 1.5, 2))
ggsave("04_Manuscript/Figures/GPSSummaries.png"
  , plot   = p
  , width  = 5
  , height = 3
  , bg     = "white"
  , scale  = 2
  , device = png
)

