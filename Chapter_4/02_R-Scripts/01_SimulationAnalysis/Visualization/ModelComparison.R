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

# The imputation approach is unaffected by forgiveness. However, since we still
# replicated the analysis for different forgiveness values, there are slight
# differences in model estimates. To avoid confusion, we'll loop through all
# cases where "approach = imputed" and "forgiveness > 1" and replace the results
# with those from a forgiveness of 1.
dat <- dat %>% nest(Results = -c(Missingness, Forgiveness, Approach))
for (i in 1:nrow(dat)) {
  appr <- as.character(dat$Approach[i])
  forg <- dat$Forgiveness[i]
  miss <- dat$Missingness[i]
  if (appr == "imputed" & forg > 1) {
    dat$Results[[i]] <- dat$Results[[which(
      as.character(dat$Approach) == appr &
      dat$Forgiveness == 1 &
      dat$Missingness == miss
    )]]
  }
}
dat <- unnest(dat, Results)

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
  mutate(SquaredError = (Estimate - Truth) ** 2) %>%
  group_by(Variable, Kernel, AutocorrRange, Missingness, Forgiveness, Approach) %>%
  summarize(RMSE = sqrt(mean(SquaredError)), .groups = "drop")

# Join the Estimates with the true values and compute the bias
dat <- dat %>%
  left_join(., truth[, -3], by = "Variable") %>%
  mutate(Bias = Estimate - Truth)

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

# Visualize estimates for different forgiveness values
p1 <- dat %>%
  subset(Missingness == 0.2 & AutocorrRange == 20) %>%
  ggplot(aes(x = as.factor(Forgiveness), y = Estimate, col = Approach, fill = Approach)) +
    geom_hline(data = truth, aes(yintercept = Truth), lty = 1, lwd = 0.5, color = "gray30") +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.5, outlier.size = 0.1, width = 0.5, linewidth = 0.5) +
    facet_nested_wrap(Kernel ~ Variable, scales = "free", labeller = label_parsed) +
    theme_minimal() +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    # labs(caption = "Scenario: Missingness = 20%, Autocorrelation-Range = 20") +
    xlab("Forgiveness") +
    ylab("Estimate") +
    theme(
        legend.position  = "bottom"
      , panel.grid.minor = element_blank()
      , strip.background = element_rect(fill = "gray95", color = "white")
      , plot.caption     = element_text(face = 3, color = "gray30")
      , axis.text        = element_text(size = 7)
    )

# # If one wants to visualize the bias (estimate - truth) instead
# p1 <- dat %>%
#   subset(Missingness == 0.2 & AutocorrRange == 20) %>%
#   ggplot(aes(x = as.factor(Forgiveness), y = Bias, col = Approach, fill = Approach)) +
#     geom_hline(yintercept = 0, lty = 1, lwd = 0.5, color = "gray30") +
#     geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.5, outlier.size = 0.1, width = 0.5, linewidth = 0.5) +
#     facet_nested_wrap(Kernel ~ Variable, scales = "free", labeller = label_parsed) +
#     theme_minimal() +
#     scale_color_viridis_d() +
#     scale_fill_viridis_d() +
#     # labs(caption = "Scenario: Missingness = 20%, Autocorrelation-Range = 20") +
#     xlab("Forgiveness") +
#     ylab("Bias (Estimate - Truth)") +
#     theme(
#         legend.position  = "bottom"
#       , panel.grid.minor = element_blank()
#       , strip.background = element_rect(fill = "gray95", color = "white")
#       , plot.caption     = element_text(face = 3, color = "gray30")
#       , axis.text        = element_text(size = 7)
#     )

# Visualize the estimates for different missingness values
p2 <- dat %>%
  subset(Forgiveness == 2 & AutocorrRange == 20) %>%
  ggplot(aes(x = as.factor(Missingness), y = Estimate, col = Approach, fill = Approach)) +
    geom_hline(data = truth, aes(yintercept = Truth), lty = 1, lwd = 0.5, color = "gray30") +
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.5, outlier.size = 0.1, width = 0.5, linewidth = 0.5) +
    facet_nested_wrap(Kernel ~ Variable, scales = "free", labeller = label_parsed) +
    theme_minimal() +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    scale_x_discrete(labels = function(x) {paste0(as.numeric(x) * 100, "%")}) +
    # labs(caption = "Scenario: Forgiveness = 2, Autocorrelation-Range = 20") +
    xlab("Missingness") +
    ylab("Estimate") +
    theme(
        legend.position  = "bottom"
      , panel.grid.minor = element_blank()
      , strip.background = element_rect(fill = "gray95", color = "white")
      , plot.caption     = element_text(face = 3, color = "gray30")
      , axis.text        = element_text(size = 7)
    )

# # If one wants to visualize the bias (estimate - truth) instead
# p2 <- dat %>%
#   subset(Forgiveness == 2 & AutocorrRange == 20) %>%
#   ggplot(aes(x = as.factor(Missingness), y = Bias, col = Approach, fill = Approach)) +
#     geom_hline(yintercept = 0, lty = 1, lwd = 0.5, color = "gray30") +
#     geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.5, outlier.size = 0.1, width = 0.5, linewidth = 0.5) +
#     facet_nested_wrap(Kernel ~ Variable, scales = "free", labeller = label_parsed) +
#     theme_minimal() +
#     scale_color_viridis_d() +
#     scale_fill_viridis_d() +
#     scale_x_discrete(labels = function(x) {paste0(as.numeric(x) * 100, "%")}) +
#     # labs(caption = "Scenario: Forgiveness = 2, Autocorrelation-Range = 20") +
#     xlab("Missingness") +
#     ylab("Bias (Estimate - Truth)") +
#     theme(
#         legend.position  = "bottom"
#       , panel.grid.minor = element_blank()
#       , strip.background = element_rect(fill = "gray95", color = "white")
#       , plot.caption     = element_text(face = 3, color = "gray30")
#       , axis.text        = element_text(size = 7)
#     )

# Visualize root mean squared error for different forgiveness values
p3 <- dat_rmse %>%
  subset(Missingness == 0.2 & AutocorrRange == 20) %>%
  ggplot(aes(x = Forgiveness, y = RMSE, color = Approach, shape = Approach)) +
    geom_line() +
    geom_point(size = 2.5) +
    facet_nested_wrap(Kernel ~ Variable, scales = "free", labeller = label_parsed) +
    scale_color_viridis_d() +
    geom_hline(yintercept = 0, color = "transparent") +
    scale_shape_manual(values = c("\u25FE", "\u25B2", "\u25BC", "\u25CF", "\u25C6")) +
    theme_minimal() +
    labs(caption = "Scenario: Missingness = 20%, Autocorrelation-Range = 20") +
    xlab("Forgiveness") +
    ylab("RMSE") +
    theme(
        legend.position  = "bottom"
      , panel.grid.minor = element_blank()
      , strip.background = element_rect(fill = "gray95", color = "white")
      , plot.caption     = element_text(face = 3, color = "gray30", hjust = 0.5)
      , axis.text        = element_text(size = 7)
    )

# Visualize root mean squared error for different missingness values
p4 <- dat_rmse %>%
  subset(Forgiveness == 2 & AutocorrRange == 20) %>%
  ggplot(aes(x = Missingness, y = RMSE, color = Approach, shape = Approach)) +
    geom_line() +
    geom_point(size = 2.5) +
    facet_nested_wrap(Kernel ~ Variable, scales = "free", labeller = label_parsed) +
    scale_color_viridis_d() +
    scale_x_continuous(labels = function(x) {paste0(as.numeric(x) * 100, "%")}) +
    geom_hline(yintercept = 0, color = "transparent") +
    scale_shape_manual(values = c("\u25FE", "\u25B2", "\u25BC", "\u25CF", "\u25C6")) +
    theme_minimal() +
    labs(caption = "Scenario: Forgiveness = 2, Autocorrelation-Range = 20") +
    xlab("Missingness") +
    ylab("RMSE") +
    theme(
        legend.position  = "bottom"
      , panel.grid.minor = element_blank()
      , strip.background = element_rect(fill = "gray95", color = "white")
      , plot.caption     = element_text(face = 3, color = "gray30", hjust = 0.5)
      , axis.text        = element_text(size = 7)
    )

# Arrange the figures for the RMSE on a single plot
p5 <- ggarrange(p1, p3, nrow = 2, align = "hv", labels = "auto")
p6 <- ggarrange(p2, p4, nrow = 2, align = "hv", labels = "auto")

# Store them
ggsave("04_Manuscript/99_ModelComparisonForgiveness.png"
  , plot   = p5
  , width  = 10
  , height = 13
  , scale  = 0.7
  , bg     = "white"
  , device = png
)
ggsave("04_Manuscript/99_ModelComparisonMissingness.png"
  , plot   = p6
  , width  = 10
  , height = 13
  , scale  = 0.7
  , bg     = "white"
  , device = png
)

################################################################################
#### Appendix Plots
################################################################################
# Aggregate the data
dat_aggr <- dat %>%
  group_by(Variable, Kernel, AutocorrRange, Missingness, Forgiveness, Approach) %>%
  summarize(
    , meanEstimate = mean(Estimate)
    , LCI_BT       = quantile(Estimate, 0.025, na.rm = T)
    , UCI_BT       = quantile(Estimate, 0.975, na.rm = T)
    , .groups      = "drop"
  ) %>%
  rename(Estimate = meanEstimate)

# Shorten labels
truth <- mutate(truth, Variable = factor(Variable, labels = c(
      TeX("$\\hat{\\kappa}$")
    , TeX("$\\hat{k}$")
    , TeX("$\\hat{\\theta}$")
    , TeX("$\\hat{\\beta}_{Dist}")
    , TeX("$\\hat{\\beta}_{Elev}$")
    , TeX("$\\hat{\\beta}_{Forest}$")
  )
))
dat_aggr <- mutate(dat_aggr, Variable = factor(Variable, labels = c(
      TeX("$\\hat{\\kappa}$")
    , TeX("$\\hat{k}$")
    , TeX("$\\hat{\\theta}$")
    , TeX("$\\hat{\\beta}_{Dist}")
    , TeX("$\\hat{\\beta}_{Elev}$")
    , TeX("$\\hat{\\beta}_{Forest}$")
  )
))

# Visualize
p <- lapply(unique(dat_aggr$AutocorrRange), function(x) {

  # Prepare the plot
  p <- dat_aggr %>%
    subset(AutocorrRange == x) %>%
    mutate(Missingness = paste0(Missingness* 100, "%")) %>%
    ggplot(aes(x = Forgiveness, y = Estimate, col = Approach, shape = Approach)) +
      geom_hline(data = truth, aes(yintercept = Truth), lty = 1, lwd = 0.4) +
      # geom_hline(yintercept = 0, lty = 2, lwd = 0.25) +
      geom_point(position = position_dodge(width = 0.8), size = 1.5) +
      geom_errorbar(aes(ymin = LCI_BT, ymax = UCI_BT), width = 0, position = position_dodge(width = 0.8), linewidth = 0.5) +
      facet_grid(Variable ~ TeX(Missingness, output = "character"), scales = "free", labeller = label_parsed) +
      theme_minimal() +
      scale_color_viridis_d() +
      scale_shape_manual(values = c("\u25FE", "\u25B2", "\u25BC", "\u25CF", "\u25C6", "\u2605")) +
      theme(
          legend.position  = "bottom"
        , panel.grid.minor = element_blank()
        , strip.background = element_rect(fill = "gray95", color = "white")
      ) +
      theme() +
      # ylab(expression(beta * "-Estimate"))
      ylab("Estimate")

  # Return it
  return(p)
})

# Put them together
all <- ggarrange(p[[1]], p[[2]], p[[3]]
  , ncol          = 1
  , labels        = "auto"
  , align         = "hv"
  , common.legend = T
  , legend        = "bottom"
)

# Store it
ggsave("04_Manuscript/99_ModelComparisonAppendix.png"
  , plot   = all
  , width  = 10
  , height = 12
  , scale  = 1.3
  , bg     = "white"
  , device = png
)
