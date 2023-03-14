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

# Prepare model files
dat <- tibble(
    Filepath = dir(path = "03_Data/ModelResults", pattern = ".rds$", full.names = T)
  , Autocorr_Range = sapply(Filepath, function(x) {
      arange <- str_split(basename(x), pattern = "\\_")[[1]][[1]]
      arange <- as.numeric(gsub(arange, pattern = "A", replacement = ""))
      return(arange)
    })
  , Replicate = sapply(Filepath, function(x) {
      rep <- str_split(basename(x), pattern = "\\_")[[1]][[2]]
      rep <- as.numeric(substr(rep, 2, 4))
      return(rep)
  })
)

# Load results
dat <- dat %>%
  mutate(Data = map(Filepath, read_rds)) %>%
  unnest(Data)

# Separate Movement from Habitat Kernel
dat$Kernel <- with(dat, ifelse(
    Coefficient %in% c("forest", "dist", "elev")
  , yes = "Habitat"
  , no  = "Movement"
))

# # For now, ignore the "model" approach
# dat <- subset(dat, Approach != "model")

# Let's also prepare a dataframe containing the "truth", i.e. the simulation
# parameters
truth <- data.frame(
    Coefficient = c("sl", "log_sl", "cos_ta", "forest", "elev", "dist")
  , Estimate    = c(0, 0, 0, -1, 0.5, -20)
  , Kernel      = c("Movement", "Movement", "Movement", "Habitat", "Habitat", "Habitat")
)

# Summarize Results
dat_aggr <- dat %>%
  group_by(Coefficient, Kernel, Autocorr_Range, Missingness, Forgiveness, Approach) %>%
  summarize(
    , meanEstimate = mean(Estimate)
    , SD           = sd(Estimate)
    # , SE           = SD / sqrt(n())
    # , LCI_SE       = meanEstimate - 1.96 * SE
    # , UCI_SE       = meanEstimate + 1.96 * SE
    # , LCI_BT       = quantile(Estimate, 0.025, na.rm = T)
    # , UCI_BT       = quantile(Estimate, 0.975, na.rm = T)
    , .groups      = "drop"
  ) %>%
  rename(Estimate = meanEstimate)

# Visualize
dat %>%
  subset(Kernel == "Habitat" & Missingness == 0.5) %>%
  ggplot(aes(x = as.factor(Forgiveness), y = Estimate, col = Approach, fill = Approach)) +
    geom_hline(data = subset(truth, Kernel == "Habitat"), aes(yintercept = Estimate), lty = 2, lwd = 0.5) +
    geom_boxplot(outlier.size = 0.1, size = 0.4) +
    facet_wrap(Autocorr_Range ~ Coefficient, scales = "free") +
    theme_minimal() +
    scale_color_viridis_d() +
    scale_fill_viridis_d(alpha = 0.2) +
    theme(
        legend.position  = "bottom"
      , strip.background = element_rect(fill = "gray95", color = NA)
    ) +
    theme() +
    ylab(expression(beta * "-Estimate"))

dat_aggr %>%
  subset(Kernel == "Habitat" & Missingness == 0.5) %>%
  ggplot(aes(x = Forgiveness, y = Estimate, col = Approach)) +
    geom_hline(data = subset(truth, Kernel == "Habitat"), aes(yintercept = Estimate), lty = 2, lwd = 0.5) +
    geom_point(position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = LCI_SE, ymax = UCI_SE), width = 0, position = position_dodge(width = 0.6)) +
    facet_wrap(Autocorr_Range ~ Coefficient, scales = "free") +
    theme_minimal() +
    scale_color_viridis_d() +
    theme(
        legend.position  = "bottom"
      , strip.background = element_rect(fill = "gray95", color = NA)
    ) +
    theme() +
    ylab(expression(beta * "-Estimate"))
