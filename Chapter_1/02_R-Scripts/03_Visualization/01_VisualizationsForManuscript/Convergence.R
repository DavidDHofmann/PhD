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

# Look for local convergence
p1 <- convergence %>%
  subset(CheckID %in% sample(unique(convergence$CheckID), 8)) %>%
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
  # geom_point(size = 0.2) +
  facet_wrap("CheckID", scales = "free", nrow = 4) +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  scale_x_continuous(
    labels = function(x){format(x, big.mark = "'")}
  ) +
  xlab("# Simulated Trajectories") +
  ylab("Relative Traversal Frequency")

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
  subset(MeanRelativeTraversals > 0) %>%
  ggplot(aes(x = NTracks, y = MeanRelativeTraversals)) +
  geom_ribbon(aes(
      ymin = Lower
    , ymax = Upper
  ), alpha = 0.5, fill = "orange", color = "orange") +
  geom_line() +
  # geom_point() +
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
  # subset(MeanRelativeTraversals > 0) %>%
  ggplot(aes(x = NTracks, y = Width)) +
  geom_hline(yintercept = 0.01, lty = 2, col = "gray40") +
  geom_line() +
  # geom_point(col = "orange", size = 0.2) +
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
  ylab("Width of 95% CI")

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
  # geom_point() +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  scale_x_continuous(
    labels = function(x){format(x, big.mark = "'")}
  ) +
  xlab("# Simulated Trajectories") +
  ylab("Mean Width of 95% CI")

# Arrange plots
p5 <- ggarrange(p2, p4, ncol = 1, labels = c("auto"))

# Store plot
ggsave("04_Manuscript/99_Convergence.png", plot = p5)
