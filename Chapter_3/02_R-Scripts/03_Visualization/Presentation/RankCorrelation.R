################################################################################
#### Rank Correlation
################################################################################
# Description: Spearman rank correlation validation results

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(tidyverse)  # To wrangle data
library(ggh4x)      # For nested facettes
library(ggpubr)     # For group comparisons
library(colorspace) # Do darken and lighten colors
library(latex2exp)  # To use latex expressions
library(ggdark)     # For dark plot themes

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Reload rank-correlations (for now, we only care about true preferences)
rank <- "03_Data/03_Results/RankFrequency.rds" %>%
  read_rds() %>%
  mutate(Formula = gsub(Formula, pattern = "Full", replacement = "Complex")) %>%
  mutate(Formula = factor(Formula, levels = c("Simple", "Complex"), labels = c("Simple Formula", "Complex Formula"))) %>%
  subset(Preferences == "True")

rank_rand <- "03_Data/03_Results/RankFrequency.rds" %>%
  read_rds() %>%
  mutate(Formula = gsub(Formula, pattern = "Full", replacement = "Complex")) %>%
  mutate(Formula = factor(Formula, levels = c("Simple", "Complex"), labels = c("Simple Formula", "Complex Formula"))) %>%
  subset(Preferences == "Randomized")

# Some summary stats
rank %>%
  group_by(Formula) %>%
  summarize(Spearman = mean(Spearman), .groups = "drop")

rank %>%
  group_by(Formula, ModelCode) %>%
  summarize(Spearman = mean(Spearman), .groups = "drop")

# Compute averages
rank_summarized <- rank %>%
  group_by(ModelCode, Formula) %>%
  summarize(Spearman = mean(Spearman), .groups = "drop")

# Compare means
mod <- aov(Spearman ~ ModelCode, data = subset(rank, Formula == "Simple Formula"))
summary(mod)

mod <- aov(Spearman ~ ModelCode, data = subset(rank, Formula == "Complex Formula"))
summary(mod)

################################################################################
#### Absolute Values
################################################################################
# Visualize results from simple model
p1 <- ggplot(subset(rank, Formula == "Simple Formula"), aes(x = ModelCode, y = Spearman, col = ModelCode, fill = ModelCode)) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) +
  # geom_boxplot(data = subset(rank_rand, Formula == "Simple Formula"), inherit.aes = F, aes(x = ModelCode, y = Spearman, group = ModelCode), width = 0.1, color = "gray90", fill = "gray90", position = position_nudge(x = 0.2), linewidth = 0.25, outlier.size = 0.2) +
  geom_boxplot(width = 0.4, alpha = 0.75) +
  theme_awesome() +
  dark_theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  scale_y_reverse(name = "Predictive Performance", limits = c(0, -1), breaks = c(0, -0.5, -1), labels = c("Bad", "Good", "Perfect"), sec.axis = sec_axis(~ ., breaks = c(0, -0.5, -1), name = "Spearman's Rank Correlation")) +
  scale_x_discrete(position = "bottom") +
  theme(
      legend.position    = "none"
    , axis.title.y       = element_text(angle = 90)
    , panel.grid.minor   = element_blank()
    , panel.grid.major   = element_line(color = adjustcolor("white", alpha.f = 0.25))
    , axis.text.y.left  = element_text(angle = 90, hjust = 0.5)
    , axis.title.y.right = element_text(angle = 90)
    , panel.background   = element_rect(fill = adjustcolor("white", alpha.f = 0.1))
    , plot.background    = element_blank()
  ) +
  xlab("") +
  ylab("Spearman's Rank Correlation")

# Visualize results from complex model
p2 <- ggplot(subset(rank, Formula == "Complex Formula"), aes(x = ModelCode, y = Spearman, col = ModelCode, fill = ModelCode)) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) +
  # geom_boxplot(data = subset(rank_rand, Formula == "Complex Formula"), inherit.aes = F, aes(x = ModelCode, y = Spearman, group = ModelCode), width = 0.1, color = "gray90", fill = "gray90", position = position_nudge(x = 0.2), linewidth = 0.25, outlier.size = 0.2) +
  geom_boxplot(width = 0.4, alpha = 0.75) +
  theme_awesome() +
  dark_theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  scale_y_reverse(name = "Predictive Performance", limits = c(0, -1), breaks = c(0, -0.5, -1), labels = c("Bad", "Good", "Perfect"), sec.axis = sec_axis(~ ., breaks = c(0, -0.5, -1), name = "Spearman's Rank Correlation")) +
  scale_x_discrete(position = "bottom") +
  theme(
      legend.position    = "none"
    , axis.title.y       = element_text(angle = 90)
    , panel.grid.minor   = element_blank()
    , panel.grid.major   = element_line(color = adjustcolor("white", alpha.f = 0.25))
    , axis.text.y.left   = element_text(angle = 90, hjust = 0.5)
    , axis.title.y.right = element_text(angle = 90)
    , panel.background   = element_rect(fill = adjustcolor("white", alpha.f = 0.1))
    , plot.background    = element_blank()
  ) +
  xlab("") +
  ylab("Spearman's Rank Correlation")

# Store the plot to file
ggsave("/home/david/ownCloud/University/15. PhD/Chapter_3/05_Presentation/RankCorrelationSimple.png"
  , plot   = p1
  , device = png
  , bg     = "transparent"
  , width  = 7.5
  , height = 2.5
  , scale  = 1.2
)

# Store the plot to file
ggsave("/home/david/ownCloud/University/15. PhD/Chapter_3/05_Presentation/RankCorrelationComplex.png"
  , plot   = p2
  , device = png
  , bg     = "transparent"
  , width  = 7.5
  , height = 2.5
  , scale  = 1.2
)
