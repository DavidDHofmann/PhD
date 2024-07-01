################################################################################
#### Overview
################################################################################
# Description: Introduction to the problem of missingness and how it may be
# addressed by increasing the forgiveness

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)     # For data wrangling
library(ggpubr)        # To combine multiple plots
library(ggarchery)     # For arrows
library(ggh4x)         # For nested facet_wrap
library(latex2exp)     # For latex expressions

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/Functions.R")

################################################################################
#### Problem of Missingness
################################################################################
# Simulate a track
set.seed(123)
df <- tibble(
    x           = sort(runif(6, min = 0, max = 100))
  , y           = runif(6, min = 0, max = 25)
  , id          = 1
  , step_number = 1:length(x)
  , step_id     = 1:length(x)
  , Successful  = T
  , timestamp   = lubridate::ymd_hms("2000-01-01 00:00:00") + hours(1:length(x))
)

# Prepare a dataframe with the different scenarios we want to depict. We'll then
# prepare the associated data step by step
all <- tibble(
    Data        = list(NA, NA, NA)
  , Scenario    = c(
      "Missingness = 0%, Forgiveness = 1, Valid Steps = 4"
    , "Missingness = 16.7%, Forgiveness = 1, Valid Steps = 1"
    , "Missingness = 16.7%, Forgiveness = 2, Valid Steps = 3"
  )
)

# Identify bursts and compute step metrics
dat <- computeBursts(df, max_duration = 1)
dat <- computeMetrics(dat)

# Compute endpoint of each step
dat$x_to = dat$x + dat$sl * sin(dat$absta)
dat$y_to = dat$y + dat$sl * cos(dat$absta)

# Steps are only valid if they have a relative turning angle
dat$valid <- if_else(!is.na(dat$relta), "valid", "invalid")

# Add to list
all$Data[[1]] <- dat

# Assume one of the fixes is missing
df_missing <- df[-4, ]

# Identify bursts and compute step metrics
dat <- computeBursts(df_missing, max_duration = 1)
dat <- computeMetrics(dat)

# Compute endpoint of each step
dat$x_to = dat$x + dat$sl * sin(dat$absta)
dat$y_to = dat$y + dat$sl * cos(dat$absta)

# Steps are only valid if they have a relative turning angle
dat$valid <- if_else(!is.na(dat$relta), "valid", "invalid")

# Add to list
all$Data[[2]] <- dat

# Assume one of the fixes is missing
df_missing <- df[-4, ]

# Identify bursts and compute step metrics
dat <- computeBursts(df_missing, max_duration = 2)
dat <- computeMetrics(dat)

# Compute endpoint of each step
dat$x_to = dat$x + dat$sl * sin(dat$absta)
dat$y_to = dat$y + dat$sl * cos(dat$absta)

# Steps are only valid if they have a relative turning angle
dat$valid <- if_else(!is.na(dat$relta), "valid", "invalid")

# Add to list
all$Data[[3]] <- dat

# Put all together and to a bit of cleaning
all <- all %>%
  unnest(Data) %>%
  mutate(label_duration = TeX(paste("\\textit{$\\Delta t =", duration, "$}"), output = "character")) %>%
  mutate(valid = factor(valid, levels = c("valid", "invalid")))

# Create pseudo-fixes for the data that we're missing
missing      <- expand_grid(df[4, ], unique(all$Scenario))
missing$Type <- "Missing"

# Define colors
cols <- hcl.colors(n = 6, palette = "viridis")

# Plot
p1 <- ggplot(all
    , mapping = aes(x = x, y = y, xend = x_to, yend = y_to, col = as.factor(burst))
  ) +
  geom_point(data = missing
    , mapping     = aes(x = x, y = y, alpha = Type)
    , inherit.aes = F
    , size        = 3
    , pch         = 20
    , col         = "gray"
  ) +
  geom_segment(lineend = "round", linewidth = 7, alpha = 0.2) +
  geom_segment(aes(linetype = valid)) +
  geom_point(size = 5) +
  geom_text(
      mapping  = aes(x = (x + x_to) / 2, y = (y + y_to) / 2, label = label_duration)
    , nudge_y  = 2
    , fontface = 3
    , size     = 3
    , col      = "black"
    , parse    = T
  ) +
  geom_text(
      mapping  = aes(label = step_number)
    , size     = 3
    , col      = "white"
  ) +
  coord_cartesian(clip = "off") +
  scale_alpha_manual(values = 1) +
  scale_color_manual(values = cols[c(1:2)], name = "Burst") +
  scale_fill_manual(values = cols[c(1:2)], name = "Burst") +
  scale_linetype_manual(values = c("solid", "22"), name = "Step Validity") +
  facet_wrap(~ Scenario, ncol = 1) +
  theme_minimal() +
  theme(
      strip.background = element_rect(fill = "gray95", color = "white")
    , panel.grid.major = element_blank()
    , axis.ticks.y     = element_line(color = "gray")
    , axis.text.x      = element_blank()
    , axis.title.y     = element_text(angle = 0, vjust = 0.5)
    , legend.position  = "none"
  )

################################################################################
#### Solution of Forgiveness
################################################################################
# Load custom functions
set.seed(99)
n <- 26

# Load simulations and pick one as an example
dat <- "03_Data/Simulation.rds" %>%
  read_rds() %>%
  subset(Replicate == 2 & AutocorrRange == 10) %>%
  dplyr::select(AutocorrRange, Movement) %>%
  unnest(Movement) %>%
  slice(1:n)

# We'll simplify the track to a straight line
dat$x <- 1
dat$y <- 1:nrow(dat)

# Rarify the trajectory
dat <- rarifyData(dat, missingness = 0.4)

# Compute bursts assuming different missingness values
dat_bursted <- lapply(1:5, function(x) {
  res <- computeBursts(dat, max_duration = x)
  res <- computeMetrics(res)
  res$x_to <- res$x + sin(res$absta) * res$sl
  res$y_to <- res$y + cos(res$absta) * res$sl
  res$valid <- ifelse(is.na(res$relta), "invalid", "valid")
  res$Forgiveness <- x
  return(res)
}) %>% do.call(rbind, .) %>% mutate(valid = factor(valid, levels = c("valid", "invalid")))

# Create pseudo-fixes for the data that we're missing
missing <- expand_grid(y = 1:(n - 1), x = 1, Forgiveness = 1:5)
missing <- anti_join(missing, dat_bursted, by = c("Forgiveness", "x", "y"))
missing$Type <- "Missing"

# Identify rectangle for each burst
burst <- dat_bursted %>%
  group_by(Forgiveness, burst) %>%
  summarize(
      ymin = min(y, na.rm = T)
    , ymax = max(y, na.rm = T)
    , .groups = "drop"
  ) %>%
  mutate(xmin = 1, xmax = 1)

# Visualize
p2 <- ggplot(dat_bursted
    , mapping = aes(x = x, y = y, xend = x_to, yend = y_to, col = as.factor(burst))
  ) +
  geom_point(data = missing
    , mapping     = aes(x = x, y = y, alpha = Type)
    , inherit.aes = F
    , size        = 3
    , pch         = 20
    , col         = "gray"
  ) +
  geom_line(aes(linetype = valid)) +
  geom_segment(data = burst
    , lineend   = "round"
    , mapping   = aes(x = xmin, xend = xmax, y = ymin, yend = ymax)
    , linewidth = 7
    , alpha     = 0.2
  ) +
  geom_point(size = 3) +
  facet_nested_wrap(~ "Forgiveness" + Forgiveness, nrow = 1) +
  scale_alpha_manual(values = 1) +
  scale_color_viridis_d(name = "Burst") +
  scale_fill_viridis_d(name = "Burst") +
  scale_linetype_manual(values = c("solid", "22"), name = "Step Validity") +
  scale_x_continuous(limits = c(0, 2), breaks = 1) +
  scale_y_continuous(breaks = 1:n) +
  theme_minimal() +
  guides(
      fill     = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)
    , color    = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)
    , linetype = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)
    , alpha    = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)

  ) +
  theme(
      strip.background = element_rect(fill = "gray95", color = "white")
    , panel.grid.minor = element_blank()
    , panel.grid.major = element_blank()
    , axis.ticks.y     = element_line(color = "gray")
    , axis.text.x      = element_blank()
    , legend.position  = "bottom"
    , legend.direction = "horizontal"
  ) +
  ylab("Animal Location-Number") +
  xlab("")

# Extract the legend
legend <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")

# Arrange plots and add legend back
p <- ggarrange(p1, p2, labels = "auto", nrow = 1, widths = c(1, 0.6))
p <- ggarrange(p, legend, nrow = 2, heights = c(1, 0.1))

# Store plot to file
ggsave("04_Manuscript/99_Overview.png"
  , plot   = p
  , bg     = "white"
  , width  = 6
  , height = 4
  , device = png
  , scale  = 1.7
)
