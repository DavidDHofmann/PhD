################################################################################
#### Valid Steps
################################################################################
# Description: Plot of number of valid steps

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)     # For data wrangling
library(ggdark)        # For dark themes

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load simulated data
dat <- read_rds("03_Data/ValidSteps.rds")

# Compute summary stats from the simulations
valid_steps <- dat %>%
  group_by(Missingness, Forgiveness) %>%
  summarize(
      LCI     = quantile(NumberSteps, 0.025)
    , UCI     = quantile(NumberSteps, 0.975)
    , Steps   = mean(NumberSteps)
    , .groups = "drop"
  ) %>%
  mutate(Metric = "Number of valid steps")

# What's the decrease in valid steps when moving from 0% missingness to 25%
# missingness?
drop <- valid_steps %>%
  subset(Missingness %in% c(0, 0.25)) %>%
  select(-c(LCI, UCI, Metric)) %>%
  mutate(Missingness = paste0("Missing_", as.character(Missingness * 100))) %>%
  pivot_wider(names_from = Missingness, values_from = Steps)
drop %>%
  mutate(DropPercent = 100 - round(Missing_25 / Missing_0 * 100)) %>%
  slice(1) %>%
  pull(DropPercent)
drop %>%
  mutate(GainPercent = round((Missing_25 / first(Missing_25) - 1) * 100)) %>%
  slice(2) %>%
  pull(GainPercent)

# Let's also compute the number of steps gained when increasing the forgiveness
gained_steps <- dat %>%
  group_by(Replicate, Missingness) %>%
  mutate(StepsGained = NumberSteps - first(NumberSteps)) %>%
  group_by(Missingness, Forgiveness) %>%
  summarize(
      LCI         = quantile(StepsGained, 0.025)
    , UCI         = quantile(StepsGained, 0.975)
    , Steps = mean(StepsGained)
    , .groups     = "drop"
  ) %>%
  mutate(Metric = "Number of steps to be gained")

# Bind all data together
all <- rbind(valid_steps, gained_steps) %>%
  mutate(Metric = factor(Metric, levels = c("Number of valid steps", "Number of steps to be gained")))

# Visualize
p <- ggplot(all, aes(x = Missingness, y = Steps, col = as.factor(Forgiveness), fill = as.factor(Forgiveness), ymin = LCI, ymax = UCI)) +
  geom_ribbon(alpha = 0.4, linewidth = 0) +
  geom_line(linewidth = 0.4) +
  geom_vline(xintercept = 1 - 0.78, lty = 2, color = "white") +
  scale_fill_viridis_d(name = "Forgiveness") +
  scale_color_viridis_d(name = "Forgiveness") +
  scale_linetype_discrete(name = "Forgiveness") +
  annotate("text"
      , label    = "Average missingness in terrestrial\nGPS studies (0.22; Hofman, 2019)"
      , x        = 0.3
      , y        = 900
      , hjust    = 0
      , size     = 3
      , fontface = 3
    ) +
  dark_theme_minimal() +
  # ylab("Number of Steps to Gain") +
  theme(
      legend.position   = "bottom", panel.grid.minor = element_blank()
    , strip.background  = element_rect(fill = "gray10", color = NA)
    , legend.key.height = unit(0.1, "cm")
    , panel.grid        = element_blank()
    # , panel.grid.major  = element_line(color = "gray10")
    , panel.background  = element_blank()
    , plot.background   = element_blank()
  ) +
  facet_wrap(~ Metric)

# Store the plot
ggsave("05_Presentation/99_NumberOfSteps.png"
  , plot   = p
  , width  = 5
  , height = 2.5
  , bg     = "transparent"
  , scale  = 1.5
  , device = png
)
