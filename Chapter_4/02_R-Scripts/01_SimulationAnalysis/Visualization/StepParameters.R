################################################################################
#### Step Parameters
################################################################################
# Description: Plot of step distribution parameters

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)  # For data wrangling
library(viridis)    # To get viridis colors
library(latex2exp)  # For latex expressions

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load model results on fitted distributions
params <- "03_Data/Analysis.rds" %>%
  read_rds() %>%
  mutate(Data = map(ResultsFilename, read_rds)) %>%
  select(-ResultsFilename) %>%
  unnest(Data) %>%
  subset(Approach == "dynamic+model") %>%
  mutate(Results = map(Results, function(x) {
    sl <- x$Dists$dynamic$sl
    ta <- x$Dists$dynamic$ta
    all <- tibble(
        duration = sl$duration
      , shape    = sl$shape
      , scale    = sl$scale
      , kappa    = ta$kappa
      # , mu       = ta$mu
    )
    return(all)
  })) %>%
  unnest(Results) %>%
  pivot_longer(shape:kappa, names_to = "Parameter", values_to = "Value") %>%
  mutate(Parameter = gsub(Parameter, pattern = "kappa", replacement = "concentration")) %>%
  mutate(Parameter = factor(Parameter
    , levels = c(
        "concentration" = "concentration"
      , "shape"         = "shape"
      , "scale"         = "scale"
    )
    , labels = c(
        "concentration" = TeX("Concentration ($\\kappa_0$)")
      , "shape"         = TeX("Shape ($k_0$)")
      , "scale"         = TeX("Scale ($\\theta_0$)")
    )
  ))

# Compute means
params_means <- params %>%
  group_by(Parameter, AutocorrRange, duration) %>%
  summarize(Value = mean(Value), .groups = "drop")

# Visualize everything
p <- params %>%
  subset(Value < 9) %>%
  ggplot(aes(x = duration, y = Value)) +
    geom_jitter(width = 0.25, alpha = 0.2, size = 0.5) +
    geom_point(data = params_means, col = viridis(10)[10], size = 3) +
    geom_line(data = params_means, col = viridis(10)[10]) +
    facet_grid(Parameter ~ AutocorrRange, scales = "free", labeller = label_parsed) +
    theme_minimal() +
    xlab("Step Duration") +
    ylab("Parameter Estimate") +
    theme(
        strip.background = element_rect(fill = "gray95", color = NA)
      , panel.grid.minor = element_blank()
    )

# Store the plot
ggsave("04_Manuscript/99_DistributionParameters.png"
  , plot   = p
  , width  = 5
  , height = 5
  , scale  = 1.4
  , bg     = "white"
  , device = png
)
