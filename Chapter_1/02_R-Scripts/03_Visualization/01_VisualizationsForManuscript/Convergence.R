################################################################################
#### Convergence
################################################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)      # To wrangle data
library(lemon)          # For capped axes in plot
library(ggpubr)         # To arrange multiple plots

# Load data
convergence <- read_rds("03_Data/03_Results/99_Convergence.rds")

# Give the checkpoint ID a nicer name
convergence$CheckID <- paste0("Checkpoint ", convergence$CheckID)

# Summarize
summarized <- convergence %>%
  group_by(NTracks, Replicate) %>%
  summarize(
      RelativeTraversals = mean(RelativeTraversals)
    , .groups            = "drop"
  ) %>%
  subset(NTracks > 0)

# Set seed
set.seed(1234)

# Look for local convergence
p1 <- convergence %>%
  subset(CheckID %in% sample(unique(convergence$CheckID), 3)) %>%
  group_by(NTracks, CheckID) %>%
  summarize(
      MeanRelativeTraversals = mean(RelativeTraversals)
    , Upper                  = quantile(RelativeTraversals, 0.975)
    , Lower                  = quantile(RelativeTraversals, 0.025)
    , .groups                = "drop"
  ) %>%
  subset(MeanRelativeTraversals > 0) %>%
  ggplot(aes(x = NTracks, y = MeanRelativeTraversals)) +
  geom_ribbon(aes(
      ymin = Lower
    , ymax = Upper
  ), alpha = 0.5, fill = "orange", color = "orange") +
  geom_line() +
  facet_wrap("CheckID", scales = "free", nrow = 3) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  scale_x_continuous(
      labels = function(x){format(x, big.mark = "'")}
    , breaks = c(0, 25000, 50000)
  ) +
  xlab("# Simulated Trajectories") +
  ylab("Relative Traversal Frequency") +
  theme(plot.margin = margin(0, 10, 0, 5))

# Look for global convergence
p2 <- convergence %>%
  group_by(NTracks, Replicate) %>%
  summarize(
      RelativeTraversals = mean(RelativeTraversals)
    , .groups            = "drop"
  ) %>%
  group_by(NTracks) %>%
  summarize(
      MeanRelativeTraversals = mean(RelativeTraversals)
    , Upper                  = quantile(RelativeTraversals, 0.975)
    , Lower                  = quantile(RelativeTraversals, 0.025)
    , .groups                = "drop"
  ) %>%
  subset(NTracks > 0) %>%
  ggplot(aes(x = NTracks, y = MeanRelativeTraversals)) +
  # geom_point(data = summarized, aes(x = NTracks, y = RelativeTraversals), col = "orange", alpha = 0.2, size = 0.5, pch = 19) +
  geom_ribbon(aes(
      ymin = Lower
    , ymax = Upper
  ), alpha = 0.5, fill = "orange", color = "orange") +
  geom_line() +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  scale_x_continuous(
    labels = function(x){format(x, big.mark = "'")}
  ) +
  xlab("# Simulated Trajectories") +
  ylab("Mean Relative Traversal Frequency")

# Check the width of the confidence interval over time
p3 <- convergence %>%
  subset(CheckID %in% sample(unique(convergence$CheckID), 8)) %>%
  group_by(NTracks, CheckID) %>%
  summarize(
      MeanRelativeTraversals = mean(RelativeTraversals)
    , Upper                  = quantile(RelativeTraversals, 0.975)
    , Lower                  = quantile(RelativeTraversals, 0.025)
    , Width                  = Upper - Lower
    , .groups                = "drop"
  ) %>%
  subset(NTracks > 0) %>%
  ggplot(aes(x = NTracks, y = Width)) +
  geom_hline(yintercept = 0.01, lty = 2, col = "gray40") +
  geom_line() +
  facet_wrap("CheckID", nrow = 4) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  scale_x_continuous(
    labels = function(x){format(x, big.mark = "'")}
  ) +
  ylim(c(0, 0.03)) +
  xlab("# Simulated Trajectories") +
  ylab("Width of 95%-PI")

# Check global confidence interval over time
p4 <- convergence %>%
  group_by(NTracks, Replicate) %>%
  summarize(
      RelativeTraversals = mean(RelativeTraversals)
    , .groups            = "drop"
  ) %>%
  group_by(NTracks) %>%
  summarize(
      MeanRelativeTraversals = mean(RelativeTraversals)
    , Upper                  = quantile(RelativeTraversals, 0.975)
    , Lower                  = quantile(RelativeTraversals, 0.025)
    , Width                  = Upper - Lower
    , .groups                = "drop"
  ) %>%
  subset(MeanRelativeTraversals > 0) %>%
  ggplot(aes(x = NTracks, y = Width)) +
  geom_line() +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  scale_x_continuous(
    labels = function(x){format(x, big.mark = "'")}
  ) +
  xlab("# Simulated Trajectories") +
  ylab("Width of 95%-PI")

# Arrange plots
p5 <- ggarrange(p2, p4, ncol = 1, labels = c("b", "c"), heights = c(1, 0.4), label.x = 0.05)
p6 <- ggarrange(p1, p5, ncol = 2, labels = c("a"), widths = c(0.5, 1))

# Show the final plot
p6

# Store plot
ggsave("04_Manuscript/99_Convergence.png", plot = p6, width = 10, height = 6)
