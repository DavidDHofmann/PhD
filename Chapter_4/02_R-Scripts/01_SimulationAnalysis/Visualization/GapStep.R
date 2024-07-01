################################################################################
#### How the Distribution of Turning-Angles Changes
################################################################################
# Description: Visualization that shows how the relative turning angle changes
# in relation to the duration of the previous step.

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)     # For data wrangling
library(pbmcapply)     # For multicore abilities

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/Functions.R")

# Load observed movement data and covariates
# obs <- read_csv("03_Data/01_RawData/ObservedMovements.csv")
dat <- "03_Data/Simulation.rds" %>%
  read_rds() %>%
  dplyr::select(AutocorrRange, Movement, Replicate)

# Go through the replicate, rarify the data, compute metrics for forgiveness = 5
dat$Reltadur <- pbmclapply(dat$Movement, ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
  relta_durs <- lapply(1:10, function(y) {
    relta_dur <- x %>%
      rarifyData(missingness = 0.5) %>%
      computeBursts(max_duration = 5) %>%
      computeMetrics() %>%
      mutate(duration_previous = lag(duration)) %>%
      dplyr::select(relta, duration, duration_previous)
    return(relta_dur)
  }) %>% do.call(rbind, .)
  return(relta_durs)
})

# Unnest and clean up
final <- dat %>%
  dplyr::select(AutocorrRange, Reltadur) %>%
  unnest(Reltadur) %>%
  subset(duration == 1) %>%
  mutate(duration = as.factor(duration)) %>%
  mutate(duration_previous = as.factor(duration_previous)) %>%
  subset(!is.na(relta) & !is.na(duration_previous))

# Plot
p1 <- ggplot(final, aes(x = relta, col = duration_previous)) +
  stat_density(geom = "line", position = "identity") +
  facet_wrap(~ AutocorrRange, nrow = 3) +
  scale_color_viridis_d(name = "Duration of Previous Step") +
  theme_minimal() +
  scale_x_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi)
    , labels = c(expression(-pi, -pi/2, 0, pi/2, pi))
  ) +
  theme(
      legend.position  = "bottom"
    , panel.grid.minor = element_blank()
    , strip.background = element_rect(fill = "gray95", color = "white")
    , legend.text      = element_text(size = 5)
    , legend.title     = element_text(size = 10)
  ) +
  xlab("Relative Turning Angle") +
  ylab("Density")

# Store it
ggsave("04_Manuscript/99_TurningAnglePreviousDuration.png"
  , plot   = p1
  , width  = 3
  , height = 4
  , bg     = "white"
  , scale  = 1.4
  , device = png
)
