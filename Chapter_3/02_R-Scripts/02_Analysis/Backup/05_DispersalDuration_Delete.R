################################################################################
#### Dispersal Duration
################################################################################
# Description: Use observed dispersal events to determine distribution of
# dispersal durations.

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_3")

# Load required packages
library(tidyverse)      # To wrangle data
library(lubridate)      # To handle dates
library(fitdistrplus)   # To fit distributions

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load the dispersal data
dat <- "03_Data/02_CleanData/Dispersers.csv" %>%
  read_csv() %>%
  subset(State == "Disperser")

# Some individuals exhibit multiple dispersal events. Let's figure out how the
# threshold we use to separate one dispersal event from another affects the
# derived dispersal durations.
durations <- data.frame(Threshold = seq(5, 30, by = 5)) %>%
  mutate(Durations = map(Threshold, function(x) {
    dat_periods <- dat %>%
      nest(Data = -ID) %>%
      mutate(Data = map(Data, function(y) {
        y <- arrange(y, Timestamp)
        y$dt = as.numeric(difftime(lead(y$Timestamp), y$Timestamp, units = "days"))
        y$Burst <- lag(y$dt > x)
        y$Burst[1] <- 1
        y$Burst <- cumsum(y$Burst)
        return(y)
      })) %>%
      unnest(Data) %>%
      group_by(ID, Burst) %>%
      summarize(
          First   = min(Timestamp)
        , Last    = max(Timestamp)
        , .groups = "drop"
      ) %>%
      mutate(Difference = difftime(Last, First, units = "days"))
    return(as.numeric(dat_periods$Difference))
  }))

# Fit gamma distribution to the durations
durations <- durations %>%
  mutate(Fitted = map(Durations, function(x) {
    fitted <- suppressWarnings({coef(fitdist(x, distr = "gamma"))})
    fitted <- data.frame(shape = fitted["shape"], scale = 1 / fitted["rate"])
    rownames(fitted) <- NULL
    return(fitted)
  }))

# Generate data so we can visualize parametric gamma distributions
gammas <- durations %>%
  dplyr::select(Threshold, Fitted) %>%
  mutate(Data = map(Fitted, function(i) {
    x <- seq(0, 300, length.out = 1000)
    y <- dgamma(x, shape = i$shape, scale = i$scale)
    xy <- data.frame(x, y)
    return(xy)
  }))

# Unnest all
durations <- unnest(durations, Durations)
gammas    <- unnest(gammas, c(Data, Fitted))

# Let's also prepare a data frame for the fitted parameters
labels <- durations %>%
  dplyr::select(Threshold, Fitted) %>%
  unnest(Fitted) %>%
  distinct() %>%
  mutate(shape = round(shape, 2)) %>%
  mutate(scale = round(scale, 2)) %>%
  mutate(Label = paste0("Shape: ", shape, ", Scale: ", scale)) %>%
  mutate(x = 150, y = 0.02)

# Visualize
ggplot() +
  geom_histogram(data = durations, aes(x = Durations, y = ..density..), fill = "gray50", col = "white", linewidth = 0.2, bins = 30) +
  geom_line(data = gammas, aes(x = x, y = y)) +
  geom_text(data = labels, aes(x = x, y = y, label = Label), size = 2) +
  facet_wrap(~ Threshold) +
  theme_minimal() +
  theme(
      strip.background = element_rect(fill = "gray95", color = "white")
    , panel.grid.minor = element_blank()
  ) +
  xlab("Dispersal Durations (in days)") +
  ylab("Density")
