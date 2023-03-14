################################################################################
#### Threshold Used to Define Activity
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/Schreibtisch")

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling

# Load data from threshold simulations
dat <- read_csv("03_Data/03_Results/99_ActiveInactiveSeparation.csv")

# Identify the minimum
min <- dat %>% subset(Minimize == min(Minimize))

# Labels to add
lab <- data.frame(
    x = c(0 - 5, min$activity)
  , y = c(min$Minimize, 0 - 0.02)
  , labels = round(c(min$Minimize, min$activity), 2)
)

# Create a plot
ggplot(dat, aes(x = activity, y = Minimize, col = as.factor(period))) +
  geom_line() +
  geom_segment(aes(x = min$activity, xend = min$activity, y = 0, yend = min$Minimize), lty = 2, col = "black") +
  geom_segment(aes(x = 0, xend = min$activity, y = min$Minimize, yend = min$Minimize), lty = 2, col = "black") +
  geom_point(data = min, aes(x = activity, y = Minimize), col = "red", inherit.aes = F) +
  geom_text(data = lab, aes(x = x, y = y, label = labels), inherit.aes = F, size = 2) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = hcl.colors(6, palette = "Spectral"), name = "Consecutive activity values considered") +
  ylab("Percentage of individuals active\nbefore 14:00 or never active")

# Store the plot
ggsave("04_Manuscript/99_ActivitySeparation.png", plot = last_plot(), bg = "white", width = 6, height = 4)
