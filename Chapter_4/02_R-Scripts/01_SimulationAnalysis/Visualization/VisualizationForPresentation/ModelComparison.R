################################################################################
#### Comparing Model Estimates
################################################################################
# Description: Comparison of estimates from different models

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)   # For plotting
library(ggpubr)      # To combine multiple plots
library(ggh4x)       # For nested facet wraps
library(ggdark)      # For dark themes
library(latex2exp)   # For latex expressions

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/Functions.R")

# Reload the model results
dat <- "03_Data/SimulationResultsConsolidated.rds" %>%
  read_rds() %>%
  select(-Results) %>%
  unnest(ResultsCleaned) %>%
  rename(Estimate = Value) %>%
  subset(Variable != "mu") %>%
  mutate(Kernel = ifelse(Kernel == "Habitat", "Habitat-Selection Function", "Movement Kernel")) %>%
  mutate(Variable = gsub(Variable, pattern = "kappa", replacement = "concentration")) %>%
  mutate(Approach = factor(Approach, levels = c(
      "uncorrected"   = "uncorrected"
    , "naive"         = "naive"
    , "dynamic+model" = "dynamic+model"
    , "multistep"     = "multistep"
    , "imputed"       = "imputed"
  )))

# Let's also prepare a dataframe containing the "truth", i.e. the simulation
# parameters
truth <- data.frame(
    Variable = c("concentration", "shape", "scale", "dist", "elev", "forest")
  , Truth    = c(0.5, 3, 1, -20, 0.5, -1)
  , Kernel   = rep(c("Movement Kernel", "Habitat-Selection Function"), each = 3)
)

# Join the Estimates with the true values and compute the root mean squared
# error
dat_rmse <- dat %>%
  left_join(., truth[, -3], by = "Variable") %>%
  mutate(SquaredError = (Truth - Estimate) ** 2) %>%
  group_by(Variable, Kernel, AutocorrRange, Missingness, Forgiveness, Approach) %>%
  summarize(RMSE = sqrt(mean(SquaredError)), .groups = "drop")

# Join the Estimates with the true values and compute the bias
dat <- dat %>%
  left_join(., truth[, -3], by = "Variable") %>%
  mutate(Bias = Truth - Estimate)

# Function to make the factor columns fancy
fancyFactors <- function(x) {
  x <- mutate(x, Variable = factor(Variable
    , levels = c(
        "concentration" = "concentration"
      , "shape"         = "shape"
      , "scale"         = "scale"
      , "dist"          = "dist"
      , "elev"          = "elev"
      , "forest"        = "forest"
    )
    , labels = c(
        "concentration" = TeX("Concentration ($\\hat{\\kappa}$)")
      , "shape"         = TeX("Shape ($\\hat{\\k}$)")
      , "scale"         = TeX("Scale ($\\hat{\\theta}$)")
      , "dist"          = TeX("$\\hat{\\beta}_{Dist}")
      , "elev"          = TeX("$\\hat{\\beta}_{Elev}$")
      , "forest"        = TeX("$\\hat{\\beta}_{Forest}$")
    )
  ))
  x <- mutate(x, Kernel = factor(Kernel
    , levels = c(
        "Habitat-Selection Function" = "Habitat-Selection Function"
      , "Movement Kernel"            = "Movement Kernel"
    )
    , labels = c(
        "Habitat-Selection Function" = TeX("Habitat-Selection Function", output = "character")
      , "Movement Kernel"            = TeX("Movement Kernel", output = "character")
    )
  ))
  return(x)
}

# Apply it
dat       <- fancyFactors(dat)
dat_rmse  <- fancyFactors(dat_rmse)
truth     <- fancyFactors(truth)

# Visualize the bias for different forgiveness values
p1 <- dat %>%
  subset(Missingness == 0.5 & AutocorrRange == 20 & Approach %in% c("uncorrected", "dynamic+model")) %>%
  ggplot(aes(x = as.factor(Forgiveness), y = Estimate, col = Approach, fill = Approach)) +
    geom_hline(data = truth, aes(yintercept = Truth), lty = 1, lwd = 0.5, color = "gray30") +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.5, outlier.size = 0.1, width = 0.5, linewidth = 0.5) +
    facet_nested_wrap(Kernel ~ Variable, scales = "free", labeller = label_parsed) +
    dark_theme_minimal() +
    scale_color_viridis_d(begin = 0.7) +
    scale_fill_viridis_d(begin = 0.7) +
    labs(caption = "Scenario: Missingness = 50%, Autocorrelation-Range = 20") +
    xlab("Forgiveness") +
    ylab("Estimate") +
    theme(
        legend.position  = "bottom"
      , panel.background = element_blank()
      , plot.background  = element_blank()
      , panel.grid       = element_blank()
      , strip.background = element_rect(fill = "gray10", color = NA)
      , plot.caption     = element_text(face = 3, color = "gray30", hjust = 0.5)
      , axis.text        = element_text(size = 7)
    )

p1_b <- dat %>%
  subset(Missingness == 0.5 & AutocorrRange == 20 & Approach %in% c("uncorrected", "dynamic+model")) %>%
  ggplot(aes(x = as.factor(Forgiveness), y = Estimate, col = Approach, fill = Approach)) +
    geom_hline(data = truth, aes(yintercept = Truth), lty = 1, lwd = 0.5, color = "gray30") +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0, outlier.size = 0.1, width = 0.5, linewidth = 0) +
    facet_nested_wrap(Kernel ~ Variable, scales = "free", labeller = label_parsed) +
    dark_theme_minimal() +
    scale_color_viridis_d(begin = 0.7) +
    scale_fill_viridis_d(begin = 0.7) +
    labs(caption = "Scenario: Missingness = 50%, Autocorrelation-Range = 20") +
    xlab("Forgiveness") +
    ylab("Estimate") +
    theme(
        legend.position  = "bottom"
      , panel.background = element_blank()
      , plot.background  = element_blank()
      , panel.grid       = element_blank()
      , strip.background = element_rect(fill = "gray10", color = NA)
      , plot.caption     = element_text(face = 3, color = "gray30", hjust = 0.5)
      , axis.text        = element_text(size = 7)
    )

# Visualize the bias for different missingness values
p2 <- dat %>%
  subset(Forgiveness == 5 & AutocorrRange == 20 & Approach %in% c("uncorrected", "dynamic+model")) %>%
  ggplot(aes(x = as.factor(Missingness), y = Estimate, col = Approach, fill = Approach)) +
    geom_hline(data = truth, aes(yintercept = Truth), lty = 1, lwd = 0.5, color = "gray30") +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.5, outlier.size = 0.1, width = 0.5, linewidth = 0.5) +
    facet_nested_wrap(Kernel ~ Variable, scales = "free", labeller = label_parsed) +
    dark_theme_minimal() +
    scale_color_viridis_d(begin = 0.7) +
    scale_fill_viridis_d(begin = 0.7) +
    scale_x_discrete(labels = function(x) {paste0(as.numeric(x) * 100, "%")}) +
    labs(caption = "Scenario: Forgiveness = 5, Autocorrelation-Range = 20") +
    xlab("Missingness") +
    ylab("Estimate") +
    theme(
        legend.position  = "bottom"
      , panel.background = element_blank()
      , plot.background  = element_blank()
      , panel.grid       = element_blank()
      , strip.background = element_rect(fill = "gray10", color = NA)
      , plot.caption     = element_text(face = 3, color = "gray30", hjust = 0.5)
      , axis.text        = element_text(size = 7)
    )

# Visualize root mean squared error for different forgiveness values
p3 <- dat_rmse %>%
  subset(Missingness == 0.5 & AutocorrRange == 20 & Approach %in% c("uncorrected", "dynamic+model")) %>%
  ggplot(aes(x = Forgiveness, y = RMSE, color = Approach, shape = Approach)) +
    geom_line() +
    geom_point(size = 2.5) +
    facet_nested_wrap(Kernel ~ Variable, scales = "free", labeller = label_parsed) +
    scale_color_viridis_d(begin = 0.7) +
    geom_hline(yintercept = 0, color = "transparent") +
    scale_shape_manual(values = c("\u25B2", "\u25BC")) +
    # scale_shape_manual(values = c("\u25FE", "\u25B2", "\u25BC", "\u25CF", "\u25C6")) +
    dark_theme_minimal() +
    labs(caption = "Scenario: Missingness = 50%, Autocorrelation-Range = 20") +
    xlab("Forgiveness") +
    ylab("RMSE") +
    theme(
        legend.position  = "bottom"
      , panel.background = element_blank()
      , plot.background  = element_blank()
      , panel.grid       = element_blank()
      , strip.background = element_rect(fill = "gray10", color = NA)
      , plot.caption     = element_text(face = 3, color = "gray30", hjust = 0.5)
      , axis.text        = element_text(size = 7)
    )

# Visualize root mean squared error for different missingness values
p4 <- dat_rmse %>%
  subset(Forgiveness == 5 & AutocorrRange == 20 & Approach %in% c("uncorrected", "dynamic+model")) %>%
  ggplot(aes(x = Missingness, y = RMSE, color = Approach, shape = Approach)) +
    geom_line() +
    geom_point(size = 2.5) +
    facet_nested_wrap(Kernel ~ Variable, scales = "free", labeller = label_parsed) +
    scale_color_viridis_d(begin = 0.7) +
    scale_x_continuous(labels = function(x) {paste0(as.numeric(x) * 100, "%")}) +
    geom_hline(yintercept = 0, color = "transparent") +
    scale_shape_manual(values = c("\u25B2", "\u25BC")) +
    # scale_shape_manual(values = c("\u25FE", "\u25B2", "\u25BC", "\u25CF", "\u25C6")) +
    dark_theme_minimal() +
    labs(caption = "Scenario: Forgiveness = 5, Autocorrelation-Range = 20") +
    xlab("Missingness") +
    ylab("RMSE") +
    theme(
        legend.position  = "bottom"
      , panel.background = element_blank()
      , plot.background  = element_blank()
      , panel.grid       = element_blank()
      , strip.background = element_rect(fill = "gray10", color = NA)
      , plot.caption     = element_text(face = 3, color = "gray30", hjust = 0.5)
      , axis.text        = element_text(size = 7)
    )

# Arrange the figures for the RMSE on a single plot
p5 <- ggarrange(p1, p3, nrow = 2, align = "hv")
p5_b <- ggarrange(p1_b, p3, nrow = 2, align = "hv")
p6 <- ggarrange(p2, p4, nrow = 2, align = "hv")

# Store them
ggsave("05_Presentation/99_ModelComparisonForgiveness.png"
  , plot   = p5
  , width  = 10
  , height = 13
  , scale  = 0.7
  , bg     = "transparent"
  , device = png
)
ggsave("05_Presentation/99_ModelComparisonForgivenessEmpty.png"
  , plot   = p5_b
  , width  = 10
  , height = 13
  , scale  = 0.7
  , bg     = "transparent"
  , device = png
)
ggsave("05_Presentation/99_ModelComparisonMissingness.png"
  , plot   = p6
  , width  = 10
  , height = 13
  , scale  = 0.7
  , bg     = "transparent"
  , device = png
)
