################################################################################
#### Comparing Model Estimates
################################################################################
# Description: Comparison of estimates from different models

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)   # For plotting
library(ggpubr)      # To combine multiple plots

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Laod all data
dat <- "03_Data/SimulationDesign.rds" %>%
  read_rds() %>%
  mutate(Coefs = map(Filename, read_rds)) %>%
  unnest(Coefs)

# Separate Movement from Habitat Kernel
dat$Kernel <- with(dat, ifelse(
    Coefficient %in% c("forest", "dist", "elev")
  , yes = "Habitat"
  , no  = "Movement"
))

# For now, ignore the "model" approach
dat <- subset(dat, Approach != "model")

# Let's also prepare a dataframe containing the "truth", i.e. the simulation
# parameters
truth <- data.frame(
    Coefficient = c("sl", "log_sl", "cos_ta", "forest", "elev", "dist")
  , Estimate    = c(0, 0, 0, -1, 0.5, -20)
  , Kernel      = c("Movement", "Movement", "Movement", "Habitat", "Habitat", "Habitat")
)

# Coefficients for a missingness of 0.5
p1 <- dat %>%
  subset(Missingness == 0.5) %>%
  subset(Kernel == "Habitat") %>%
  group_by(Coefficient, Forgiveness, Approach) %>%
  summarize(
    , SD       = sd(Estimate)
    , Estimate = mean(Estimate)
    , LCI      = Estimate - 2 * SD
    , UCI      = Estimate + 2 * SD
    , .groups  = "drop"
  ) %>%
  ggplot(aes(x = Forgiveness, y = Estimate, col = Approach)) +
    geom_hline(data = subset(truth, Kernel == "Habitat"), aes(yintercept = Estimate), lty = 2, lwd = 0.5) +
    geom_point(position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0, position = position_dodge(width = 0.6)) +
    # geom_point(data = truth, col = "black") +
    facet_wrap(~ Coefficient, scales = "free") +
    theme_minimal() +
    scale_color_viridis_d() +
    theme(
        legend.position  = "bottom"
      , strip.background = element_rect(fill = "gray95", color = NA)
    ) +
    theme() +
    ylab(expression(beta * "-Estimate"))

# Coefficients for a missingness of 0.5
p2 <- dat %>%
  subset(Missingness == 0.5) %>%
  subset(Kernel == "Movement") %>%
  group_by(Coefficient, Forgiveness, Approach) %>%
  summarize(
    , SD       = sd(Estimate)
    , Estimate = mean(Estimate)
    , LCI      = Estimate - 2 * SD
    , UCI      = Estimate + 2 * SD
    , .groups  = "drop"
  ) %>%
  ggplot(aes(x = Forgiveness, y = Estimate, col = Approach)) +
    geom_hline(data = subset(truth, Kernel == "Movement"), aes(yintercept = Estimate), lty = 2, lwd = 0.5) +
    geom_point(position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0, position = position_dodge(width = 0.6)) +
    # geom_point(data = truth, col = "black") +
    facet_wrap(~ Coefficient, scales = "free") +
    theme_minimal() +
    scale_color_viridis_d() +
    theme(
        legend.position  = "bottom"
      , strip.background = element_rect(fill = "gray95", color = NA)
    ) +
    theme() +
    ylab(expression(beta * "-Estimate"))

# Store plots to file
ggsave("04_Manuscript/99_ResultsHabitatKernel.png"
  , plot   = p1
  , bg     = "white"
  , height = 5
  , width  = 10
  , scale  = 0.8
)
ggsave("04_Manuscript/99_ResultsMovementKernel.png"
  , plot   = p2
  , bg     = "white"
  , height = 5
  , width  = 10
  , scale  = 0.8
)

# ################################################################################
# #### LEGACY
# ################################################################################
# # Let's start with the "base" model and look at its results
# dat %>%
#   subset(Missingness == 0.4 & Forgiveness == 5) %>%
#   group_by(Coefficient, Approach) %>%
#   summarize(
#     , SD       = sd(Estimate)
#     , Estimate = mean(Estimate)
#     , LCI      = Estimate - 2 * SD
#     , UCI      = Estimate + 2 * SD
#     , .groups  = "drop"
#   ) %>%
#   ggplot(aes(x = Approach, y = Estimate)) +
#     geom_hline(data = truth, aes(yintercept = Estimate), lty = 2, lwd = 0.5) +
#     geom_point(col = "red") +
#     geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) +
#     # geom_point(data = truth, col = "black") +
#     facet_wrap(~ Coefficient, scales = "free") +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45))
#
# # Visualize Results
# dat %>%
#   subset(Kernel == "Habitat") %>%
#   group_by(Missingness, Forgiveness, Approach, Coefficient) %>%
#   summarize(
#     , SD       = sd(Estimate)
#     , Estimate = mean(Estimate)
#     , LCI      = Estimate - 2 * SD
#     , UCI      = Estimate + 2 * SD
#     , .groups  = "drop"
#   ) %>%
#   ggplot(aes(x = Missingness, y = Estimate, col = as.factor(Forgiveness), fill = as.factor(Forgiveness))) +
#     geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5, lwd = 0.5) +
#     geom_hline(data = subset(truth, Kernel == "Habitat"), aes(yintercept = Estimate), lty = 2, lwd = 0.5) +
#     # geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.05, alpha = 0.5) +
#     # geom_point() +
#     # geom_line() +
#     facet_wrap(~ Approach + Coefficient, scales = "free", nrow = 3) +
#     theme_minimal() +
#     theme(legend.position = "bottom") +
#     scale_fill_viridis_d(name = "Forgiveness") +
#     scale_color_viridis_d(name = "Forgiveness") +
#     scale_x_continuous(breaks = seq(0, 0.5, by = 0.1)) +
#     theme(axis.text.x = element_text(angle = 45))
#     # scale_fill_manual(values = c("orange", "cornflowerblue")) +
#     # scale_color_manual(values = c("orange", "cornflowerblue"))
