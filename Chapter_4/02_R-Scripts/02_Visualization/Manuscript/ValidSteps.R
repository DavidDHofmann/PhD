################################################################################
#### Valid Steps
################################################################################
# Description: Plot of number of valid steps

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)     # For data wrangling

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load simulated data
dat <- read_csv("03_Data/03_Results/ValidSteps.csv")

# Compute summary stats from the simulations
summaries <- dat %>%
  group_by(Missingness, Forgiveness) %>%
  summarize(
      LCI         = quantile(NumberSteps, 0.025)
    , UCI         = quantile(NumberSteps, 0.975)
    , NumberSteps = mean(NumberSteps)
    , .groups     = "drop"
  )

# Visualize
p <- ggplot(summaries, aes(x = Missingness, y = NumberSteps, col = as.factor(Forgiveness), fill = as.factor(Forgiveness), ymin = LCI, ymax = UCI)) +
  geom_ribbon(alpha = 0.2, lwd = 0.4) +
  geom_line(size = 0.4) +
  scale_fill_viridis_d(name = "Forgiveness") +
  scale_color_viridis_d(name = "Forgiveness") +
  theme_minimal() +
  ylab("Number of Valid Steps") +
  theme(legend.position = "bottom")
  # geom_point() +
  # geom_line() +
  # geom_errorbar(width = 0.05)

# Store the plot
ggsave("04_Manuscript/99_NumberOfSteps.png"
  , plot   = p
  , width  = 5
  , height = 3
  , bg     = "white"
  , scale  = 1.2
)
