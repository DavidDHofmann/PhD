################################################################################
#### Overview
################################################################################
# Description: Overview of the Jargon used in the paper

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)     # For data wrangling
library(ggpubr)        # To combine multiple plots

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Function to identify bursts
computeBursts <- function(data, forgiveness) {

  # Nest data by id (a burst cannot expand across multiple ids)
  data_bursted <- data %>%
    group_by(ID) %>%
    nest() %>%

    # Compute bursts by id. A new burst is defined if the duration is >
    # forgiveness
    mutate(data = map(data, function(x) {
      x$duration <- lead(x$step_number) - x$step_number
      x$irregular <- x$duration > forgiveness
      x$burst <- NA
      x$burst[1] <- 1
      x$burst[2:nrow(x)] <- lag(cumsum(x$irregular) + 1)[-1]
    return(x)

  # Unnest the data again
  })) %>% unnest(data)

  # Return the "bursted" data
  return(data_bursted)
}

# Function to compute step metrics
computeMetrics <- function(data) {

  # Nest by ID and burst
  data_metrics <- data %>%
    group_by(ID, burst) %>%
    nest() %>%

    # Compute step metrics
    mutate(data = map(data, function(z) {
      metrics <- stepMet(x = z$x, y = z$y)
      metrics <- cbind(z, metrics)
      return(metrics)
    })) %>%

    # Unnest again and tidy up
    unnest(data) %>%
    dplyr::select(burst, step_number, step_id, everything()) %>%
    ungroup()

  # Return the data containing step metrics
  return(data_metrics)
}

# Simulate a track
set.seed(123)
df <- tibble(
    x           = sort(runif(6, min = 0, max = 100))
  , y           = runif(6, min = 0, max = 25)
  , ID          = 1
  , step_number = 1:length(x)
  , step_id     = 1:length(x)
  , Successful  = T
)

################################################################################
#### Scenario 1: Missingness = 0, Forgiveness = 1
################################################################################
# Identify bursts and compute step metrics
dat <- computeBursts(df, forgiveness = 1)
dat <- computeMetrics(dat)

# Compute endpoint of each step
dat$x_to = dat$x + dat$sl * sin(dat$absta)
dat$y_to = dat$y + dat$sl * cos(dat$absta)

# Steps are only valid if they have a relative turning angle
dat$valid <- !is.na(dat$relta)

# Visualize the steps
p1 <- ggplot(df, aes(x = x, y = y)) +
  geom_segment(
      data    = dat
    , mapping = aes(xend = x_to, yend = y_to, col = valid)
  ) +
  geom_point() +
  coord_equal(ylim = c(8, 26), clip = "off") +
  geom_text(
      data     = dat
    , mapping  = aes(x = (x + x_to) / 2, y = (y + y_to) / 2, label = paste0("(", duration, ")"))
    , nudge_y  = 2
    , fontface = 3
    , size     = 3
    , col      = "gray50"
  ) +
  geom_text(
      data     = df
    , mapping  = aes(label = step_number)
    , nudge_y  = c(2, -2, 2, 2, -2, 2)
    , size     = 3
    , col      = "gray20"
  ) +
  geom_bracket(
      xmin       = df$x[1]
    , xmax       = df$x[nrow(df)]
    , y.position = 25
    , label      = "Burst 1"
    , col        = "gray50"
    , fontface   = 3
    , label.size = 3.5
    , vjust      = -0.2
  ) +
  scale_color_viridis_d(begin = 0.7, end = 1, direction = -1, name = "Valid Step") +
  theme_classic() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

################################################################################
#### Scenario 2: Missingness = 1 / 6, Forgiveness = 1
################################################################################
# Assume one of the fixes is missing
df_missing <- df[-4, ]

# Identify bursts and compute step metrics
dat <- computeBursts(df_missing, forgiveness = 1)
dat <- computeMetrics(dat)

# Compute endpoint of each step
dat$x_to = dat$x + dat$sl * sin(dat$absta)
dat$y_to = dat$y + dat$sl * cos(dat$absta)

# Steps are only valid if they have a relative turning angle
dat$valid <- !is.na(dat$relta)

# Visualize the steps
p2 <- ggplot(df, aes(x = x, y = y)) +
  geom_segment(
      data    = dat
    , mapping = aes(xend = x_to, yend = y_to, col = valid)
  ) +
  geom_point() +
  geom_point(data = df[4, ], col = "red", shape = 13, size = 4) +
  coord_equal(ylim = c(8, 26), clip = "off") +
  geom_text(
      data     = dat
    , mapping  = aes(x = (x + x_to) / 2, y = (y + y_to) / 2, label = paste0("(", duration, ")"))
    , nudge_y  = 2
    , fontface = 3
    , size     = 3
    , col      = "gray50"
  ) +
  geom_text(
      data     = df
    , mapping  = aes(label = step_number)
    , nudge_y  = c(2, -2, 2, 2, -2, 2)
    , size     = 3
    , col      = "gray20"
  ) +
  geom_bracket(
      xmin       = df$x[1]
    , xmax       = df$x[3]
    , y.position = 25
    , label      = "Burst 1"
    , col        = "gray50"
    , fontface   = 3
    , label.size = 3.5
    , vjust      = -0.2
  ) +
  geom_bracket(
      xmin       = df$x[5]
    , xmax       = df$x[6]
    , y.position = 25
    , label      = "Burst 2"
    , col        = "gray50"
    , fontface   = 3
    , label.size = 3.5
    , vjust      = -0.2
  ) +
  scale_color_viridis_d(begin = 0.7, end = 1, direction = -1, name = "Valid Step") +
  theme_classic() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

################################################################################
#### Scenario 3: Missingness = 1 / 6, Forgiveness = 2
################################################################################
# Assume one of the fixes is missing
df_missing <- df[-4, ]

# Identify bursts and compute step metrics
dat <- computeBursts(df_missing, forgiveness = 2)
dat <- computeMetrics(dat)

# Compute endpoint of each step
dat$x_to = dat$x + dat$sl * sin(dat$absta)
dat$y_to = dat$y + dat$sl * cos(dat$absta)

# Steps are only valid if they have a relative turning angle
dat$valid <- !is.na(dat$relta)

# Visualize the steps
p3 <- ggplot(df, aes(x = x, y = y)) +
  geom_segment(
      data    = dat
    , mapping = aes(xend = x_to, yend = y_to, col = valid)
  ) +
  geom_point() +
  geom_point(data = df[4, ], col = "red", shape = 13, size = 4) +
  coord_equal(ylim = c(8, 26), clip = "off") +
  geom_text(
      data     = dat
    , mapping  = aes(x = (x + x_to) / 2, y = (y + y_to) / 2, label = paste0("(", duration, ")"))
    , nudge_y  = 2
    , fontface = 3
    , size     = 3
    , col      = "gray50"
  ) +
  geom_text(
      data     = df
    , mapping  = aes(label = step_number)
    , nudge_y  = c(2, -2, 2, 2, -2, 2)
    , size     = 3
    , col      = "gray20"
  ) +
  geom_bracket(
      xmin       = df$x[1]
    , xmax       = df$x[nrow(df)]
    , y.position = 25
    , label      = "Burst 1"
    , col        = "gray50"
    , fontface   = 3
    , label.size = 3.5
    , vjust      = -0.2
  ) +
  scale_color_viridis_d(begin = 0.7, end = 1, direction = -1, name = "Valid Step") +
  theme_classic() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

################################################################################
#### Combine into Plot
################################################################################
# Arrange plots
p <- ggarrange(p1, p2, p3
  , ncol          = 1
  , common.legend = T
  , legend        = "bottom"
  , labels        = c(
      "(a) Missingness = 0 / 0, Forgiveness = 1"
    , "(b) Missingness = 1 / 6, Forgiveness = 1"
    , "(c) Missignness = 1 / 6, Forgiveness = 2"
  )
)

# Store it
ggsave("04_Manuscript/99_Overview.png"
  , plot   = p
  , width  = 10
  , height = 10
  , scale  = 0.7
  , bg     = "white"
)
