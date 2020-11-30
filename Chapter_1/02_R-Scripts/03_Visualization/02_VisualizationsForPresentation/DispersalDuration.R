################################################################################
#### Plot of Dispersal Durations
################################################################################
# Description: A simple plot of the Dispersal Durations

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "C:/Users/david/switchdrive/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)    # For data wrangling
library(lubridate)    # To handle dates nicely
library(ggpubr)       # For nice plots
library(fitdistrplus) # To fit a distribution

################################################################################
#### Data Preparation
################################################################################
# Load required data
cut <- read_csv("03_Data/02_CleanData/00_General_Dispersers_Popecol_CutoffDates.csv")

# Add Two hours to the timestamps
cut$StartDate <- cut$StartDate + hours(2)
cut$EndDate <- cut$EndDate + hours(2)

# Calculate dispersal duration
durations <- cut %>%
  mutate(DispersalDuration = EndDate - StartDate) %>%
  group_by(DogName) %>%
  summarize(DispersalDuration = as.numeric(sum(DispersalDuration)))

# Fit an exponential distribution to the durations
dist <- fitdist(durations$DispersalDuration, distr = "gamma")

# Create new data from the distribution
x <- seq(0, 300, by = 0.1)
y <- dgamma(x, shape = dist$estimate[["shape"]], rate = dist$estimate[["rate"]])

# Visualize
p <- ggdensity(
    durations
  , x     = "DispersalDuration"
  , fill  = "orange"
  , color = "orange"
  # , add   = "mean"
  , rug   = TRUE
  , xlab  = "Dispersal Duration (days)"
  , ylab  = "Density"
  , font.label = list(color = "black")
  ) +
  geom_line(data = data.frame(x, y), aes(x = x, y = y), lty = 2, col = "gray50") +
  labs(color = "white") +
  theme(
      axis.ticks       = element_line(color = "white")
    , axis.line        = element_line(color = "white")
    , axis.text        = element_text(color = "white", size = 10)
    , text             = element_text(color = "white")
    , panel.background = element_rect(fill = "black")
    , plot.background  = element_rect(fill = "black", color = "black")
  )
ggsave(p, device = "png", filename = "test.png", bg = "transparent", width = 3, height = 2)
