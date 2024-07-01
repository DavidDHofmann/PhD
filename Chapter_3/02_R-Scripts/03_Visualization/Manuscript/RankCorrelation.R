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
# Visualize
p1 <- ggplot(rank, aes(x = ModelCode, y = Spearman, col = ModelCode, fill = ModelCode)) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.1) +
  geom_boxplot(data = rank_rand, inherit.aes = F, aes(x = ModelCode, y = Spearman, group = ModelCode), width = 0.1, color = "gray90", fill = "gray90", position = position_nudge(x = 0.2), linewidth = 0.25, outlier.size = 0.2) +
  geom_boxplot(width = 0.4, alpha = 0.2) +
  stat_compare_means(
      method  = "anova"
    , label.y = 1.05
    , label.x = 4.5
    , size    = 3
  ) +
  stat_compare_means(
      aes(label = after_stat(p.signif))
    , comparisons   = list(c("SSS", "SMS"), c("SSS", "DSS"), c("SSS", "DSD"), c("DSS", "DSD"), c("DMS", "DMD"), c("DSD", "DMD"), c("SSS", "DMD"))
    , method        = "t.test"
    , tip.length    = 0
    , step.increase = 0.05
    , size          = 1.5
  ) +
  theme_awesome() +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  scale_y_reverse(breaks = c(0, -0.5, -1), sec.axis = sec_axis(~ ., breaks = c(-1, -0.5, 0), labels = c("Perfect", "Good", "Bad"), name = "Predictive Performance")) +
  scale_x_discrete(position = "bottom") +
  facet_grid(~ Formula) +
  theme(
      legend.position    = "none"
    , axis.title.y       = element_text(angle = 90)
    , panel.grid.minor   = element_blank()
    , axis.text.y.right  = element_text(angle = 90, hjust = 0.5)
    , axis.title.y.right = element_text(angle = 90)
    # , panel.spacing.x = unit(0, "line")
    # , panel.spacing.y  = unit(1.5, "line")
  ) +
  # xlab("Configuration (Dynamism)") +
  xlab("") +
  ylab("Spearman's Rank Correlation")

# Store the plot to file
ggsave("04_Manuscript/Figures/RankCorrelation.png"
  , plot   = p1
  , device = png
  , bg     = "white"
  , width  = 5
  , height = 3
  , scale  = 1.4
)

################################################################################
#### Differences
################################################################################
# Compute differences
diffs <- expand_grid(
    From    = unique(rank$ModelCode)
  , To      = unique(rank$ModelCode)
  , Formula = unique(rank$Formula)
  , Difference = NA
)
for (i in 1:nrow(diffs)) {
  value1 <- subset(rank_summarized, ModelCode == diffs$From[i] & Formula == diffs$Formula[i])$Spearman
  value2 <- subset(rank_summarized, ModelCode == diffs$To[i] & Formula == diffs$Formula[i])$Spearman
  diffs$Difference[i] <- value2 - value1
}

# Plot for simple model
p2 <- ggplot(subset(diffs, Formula == "Simple Formula"), aes(x = From, y = To, fill = Difference, color = Difference)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = round(Difference, 2)), size = 2) +
  scale_fill_gradient2(
      low    = "cornflowerblue"
    , mid    = "white"
    , high   = "orange"
    , breaks = c(-0.1, 0, 0.1)
    , name   = TeX("$\\Delta r_s = r_s^{DMD} - r_s^{SSS}$")
    , guide  = guide_colorbar(
      , title.position = "bottom"
      , title.hjust    = 0.5
      , ticks          = F
      , barheight      = unit(0.2, "cm")
    )
  ) +
  scale_color_gradient2(
      low    = darken("cornflowerblue", amount = 0.4)
    , mid    = darken("white", amount = 0.4)
    , high   = darken("orange", amount = 0.4)
  ) +
  guides(color = "none") +
  theme_awesome() +
  theme(
    , axis.title.y    = element_text(angle = 90)
  ) +
  facet_wrap(~ "Simple Formula")

# Plot for full model
p3 <- ggplot(subset(diffs, Formula == "Complex Formula"), aes(x = From, y = To, fill = Difference, color = Difference)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = round(Difference, 2)), size = 2) +
  scale_fill_gradient2(
      low    = "cornflowerblue"
    , mid    = "white"
    , high   = "orange"
    , breaks = c(-0.02, 0, 0.02)
    , name   = TeX("$\\Delta r_s = r_s^{DMD} - r_s^{SSS}$")
    , guide  = guide_colorbar(
      , title.position = "bottom"
      , title.hjust    = 0.5
      , ticks          = F
      , barheight      = unit(0.2, "cm")
    )
  ) +
  scale_color_gradient2(
      low    = darken("cornflowerblue", amount = 0.4)
    , mid    = darken("white", amount = 0.4)
    , high   = darken("orange", amount = 0.4)
  ) +
  guides(color = "none") +
  theme_awesome() +
  theme(
      axis.text.y  = element_blank()
    , axis.title.y = element_blank()
    , axis.ticks.y = element_blank()
  ) +
  facet_wrap(~ "Complex Formula")

# Combine them
p4 <- ggarrange(p2, p3, align = "hv")

# Store the plot to file
ggsave("04_Manuscript/Figures/RankCorrelationDifferences.png"
  , plot   = p4
  , device = png
  , bg     = "white"
  , width  = 5
  , height = 3
  , scale  = 1.4
)
