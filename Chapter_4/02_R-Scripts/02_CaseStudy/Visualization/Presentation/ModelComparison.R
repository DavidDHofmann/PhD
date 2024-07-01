################################################################################
#### Comparing Model Estimates
################################################################################
# Description: Comparison of estimates from different models. We will first show
# the model relating to a duration of 2 hours (the reference) and then also the
# interaction terms for habitat covariates.

# Clear R's brain
rm(list = ls())

# Suppress scientific notation
options(scipen = 999)

# Load required packages
library(kableExtra)  # To generate a nice talbe
library(survival)    # To work with conditional logistic regression model
library(tidyverse)   # For plotting
library(ggpubr)      # To combine multiple plots
library(ggh4x)       # For nested facet wraps
library(latex2exp)   # For latex expressions
library(ggdark)      # For dark themes

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/Functions.R")

# Reload the model results
dat <- read_rds("03_Data/CaseStudyResults.rds")

# Give a model name to each treatment
dat <- mutate(dat, Name = case_when(
    Forgiveness == 3 & !HabitatDurationInteractions ~ "F3-S"
  , Forgiveness == 3 & HabitatDurationInteractions ~ "F3-SH"
  , .default = "F1"
))

################################################################################
#### Results Table
################################################################################
# Compile a table that compares the three models
dat$Table <- lapply(dat$Model, function(x) {

  # Extract model coefficients
  coefs <- summary(x)$coefficients
  coefs <- tibble(
      Coefficient  = rownames(coefs)
    , Estimate     = coefs[, "coef"]
    , SE           = coefs[, "se(coef)"]
    , Pvalue       = coefs[, "Pr(>|z|)"]
  )

  # Create significance stars
  coefs <- mutate(coefs, Significance = case_when(
      Pvalue <= 0.01 ~ "***"
    , Pvalue <= 0.05 ~ "**"
    , Pvalue <= 0.10 ~ "*"
    , .default = ""
  ))

  # Round values
  coefs <- mutate(coefs
    , Estimate = paste0(round(Estimate, 5), Significance)
    , SE       = paste0("(", round(SE, 5), ")")
  )

  # Pivot
  coefs <- coefs %>%
    dplyr::select(Coefficient, Estimate, SE) %>%
    pivot_longer(Estimate:SE, names_to = "Variable", values_to = "Value")
})

################################################################################
#### Figure for Regular Step-Duration
################################################################################
# Reformat results a bit
dat <- dat %>%
  select(-Model) %>%
  unnest(Results) %>%
  rename(Estimate = Value) %>%
  mutate(Kernel = case_when(
      Kernel   == "Habitat" ~ "Habitat-Selection Function"
    , Kernel   == "Movement" ~ "Movement Kernel"
  )) %>%
  mutate(Variable = gsub(Variable, pattern = "kappa", replacement = "concentration")) %>%
  pivot_longer(LCI_90:UCI_99) %>%
  separate(name, into = c("Limit", "Level"), sep = "_", convert = T) %>%
  spread(key = Limit, value = value)

# Function to make the factor columns fancy
fancyFactors <- function(x) {
  x <- mutate(x, Variable = factor(Variable
    , levels = c(
        "concentration"   = "concentration"
      , "shape"           = "shape"
      , "scale"           = "scale"
      , "Water"           = "Water"
      , "DistanceToWater" = "DistanceToWater"
      , "Trees"           = "Trees"
    )
    , labels = c(
        "concentration"   = TeX("Concentration ($\\hat{\\kappa}$)")
      , "shape"           = TeX("Shape ($\\hat{\\k}$)")
      , "scale"           = TeX("Scale ($\\hat{\\theta}$)")
      , "Water"           = TeX("$\\hat{\\beta}_{Water}")
      , "DistanceToWater" = TeX("$\\hat{\\beta}_{DistanceToWater}$")
      , "Trees"           = TeX("$\\hat{\\beta}_{Trees}$")
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
dat <- fancyFactors(dat)

# Visualize results for a duration of 2 (the reference)
p <- ggplot(dat, aes(x = Name, y = Estimate, col = Name)) +
  geom_errorbar(
      mapping  = aes(ymin = LCI, ymax = UCI, linewidth = factor(Level))
    , width    = 0
    , alpha    = 0.5
  ) +
  geom_point(size = 2) +
  facet_nested_wrap(~ Kernel + Variable, scales = "free", ncol = 3, labeller = label_parsed) +
  dark_theme_minimal() +
  scale_color_viridis_d(begin = 0.25, end = 0.75, name = "Model") +
  scale_fill_viridis_d(begin = 0.25, end = 0.75, name = "Model") +
  scale_linewidth_manual(
      name   = "Confidence Level"
    , values = c(2, 1, 0.6)
    , guide  = "none"
  ) +
  theme(
      legend.position  = "bottom"
    , panel.grid.minor = element_blank()
    , panel.grid.major = element_line(color = adjustcolor("black", alpha.f = 0.25))
    , strip.background = element_rect(fill = "gray30", color = NA)
    , plot.caption     = element_text(face = 3, color = "gray30")
    , axis.text        = element_text(size = 7)
    , plot.background  = element_blank()
    , panel.background = element_blank()
  )

# Store it
ggsave("05_Presentation/99_CaseStudy.png"
  , plot   = p
  , width  = 10
  , height = 6.5
  , scale  = 0.8
  , bg     = "transparent"
  , device = png
)
